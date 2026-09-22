using System.Collections.Concurrent;
using System.Diagnostics;
using System.Globalization;
using System.Numerics;
using System.Text;
using OpenTabletDriver.Plugin;
using OpenTabletDriver.Plugin.Attributes;
using OpenTabletDriver.Plugin.Output;
using OpenTabletDriver.Plugin.Tablet;
using TrajectoryEigenFilterCore;

namespace TrajectoryEigenFilter.OTD;

[PluginName("A Spectral Reconstructor")]
public sealed class TrajectoryEigenFilterPlugin : AsyncPositionedPipelineElement<IDeviceReport>
{
    private const float DefaultOutputFrequencyHz = 1000f;
    private readonly object gate = new();
    private TrajectoryFilter? filter;
    private FilterSettings? activeSettings;
    private readonly Queue<Vector2> rawHistory = new();

    private long inputSequence, outputSequence;
    private long lastArrivalTicks, lastOutputTicks, rateAnchorTicks, rateAnchorSequence;
    private double estimatedRate = 700.0, lastOutputIntervalMs = double.NaN;
    private bool rateInitialized;
    private int stableRateWindows;
    private double rateErrorEma;

    // The filter operates in uniformly-spaced sample-index space at the estimated
    // average input rate. Wall-clock arrival time is used only for rate estimation,
    // pen-stream discontinuity detection, and selecting the 1000 Hz output position.
    private FilteredWindow? latestWindow;
    private double latestInputWallMs;
    private Vector2 latestRaw, latestOutput;
    private bool haveRaw, haveOutput;
    private double lastEvalIndex = double.NaN;
    private double lastInputStalenessMs = double.NaN;
    private long streamResets;

    private Task<ModelBuildResult>? modelBuildTask;
    private double requestedModelRate;
    private double lastModelBuildMs = double.NaN;
    private int modelSwaps;

    private readonly ConcurrentQueue<LogRow> logRows = new();
    private CancellationTokenSource? logCts;
    private Task? logTask;
    private volatile bool loggingFailed;
    private int droppedLogRows;
    private float resamplingFrequencyHz = DefaultOutputFrequencyHz;

    public override PipelinePosition Position => PipelinePosition.PreTransform;

    [Property("Resampling frequency"), Unit("Hz"), DefaultPropertyValue(1000.0f), ToolTip("Output frequency used in async mode.\nDefault: 1000 Hz. Recommended: 125, 250, 500, or 1000 Hz to use OTD's native timer. Other values use OTD's fallback timer.")]
    public new float Frequency
    {
        get => resamplingFrequencyHz;
        set
        {
            resamplingFrequencyHz = float.IsFinite(value)
                ? Math.Clamp(value, 1f, 4000f)
                : DefaultOutputFrequencyHz;
            base.Frequency = resamplingFrequencyHz;
        }
    }

    [BooleanProperty("Enable async mode", "Resamples output on OTD's scheduler instead of emitting only on tablet reports."), DefaultPropertyValue(true), ToolTip("Emits reconstructed positions at the configured Resampling frequency. When disabled, one reconstructed output is emitted for each tablet report.\nDefault: On. Recommended: On for resampling.")]
    public bool EnableAsyncMode { get; set; } = true;

    [Property("Window duration"), Unit("ms"), DefaultPropertyValue(100.0), ToolTip("Trajectory history used by the filter.\nDefault: 100 ms. Recommended: 70-150 ms.\n100 ms was used for current tuning; larger windows cost more model-build CPU and memory.")]
    public double WindowMilliseconds { get; set; } = 100.0;

    [Property("Filter strength"), DefaultPropertyValue(1.0), ToolTip("Inverse retained-mode control. Strength 1 uses 12 modes at zero latency; 0 retains all modes. Higher values retain fewer modes.\nDefault: 1. Recommended: 1 for high-rate tablets; 1.5 for 133 Hz tablets.")]
    public double FilterStrength { get; set; } = 1.0;

    [BooleanProperty("Use adaptive modes", "Varies the retained mode count with reconstruction look-ahead."), DefaultPropertyValue(true), ToolTip("Uses the sample-look-ahead schedule 12, 9, 8, 7, 6, then 5, scaled by Filter strength.\nDefault: On. Recommended: On.")]
    public bool UseAdaptiveModes { get; set; } = true;

    [Property("Outlier threshold"), DefaultPropertyValue(6.0), ToolTip("Robust residual threshold for isolated bad samples.\nDefault: 6. Recommended: 5-8.\nLower values reject more aggressively; higher values reserve replacement for more extreme excursions.")]
    public double OutlierThreshold { get; set; } = 6.0;

    [Property("Latency"), Unit("samples"), DefaultPropertyValue(4), ToolTip("Input samples buffered before output.\nDefault: 4 samples. Recommended: 1-10 samples.\nZero latency is fully functional but not recommended: noise rejection and input-jitter tolerance are reduced. Linear extrapolation is used if the output timer runs ahead of the latest input report.")]
    public int LatencySamples { get; set; } = 4;

    [Property("Staleness timeout"), Unit("ms"), DefaultPropertyValue(25.0), ToolTip("Stops the 1000 Hz output scheduler when no physical tablet report has arrived for this long.\nDefault: 25 ms.\nThis prevents a cached in-range report from being extrapolated after the pen leaves proximity. Output resumes on the next physical report.")]
    public double StalenessTimeoutMs { get; set; } = 25.0;


    [BooleanProperty("Enable CSV logging", "Writes diagnostic data to CSV."), DefaultPropertyValue(false), ToolTip("Writes diagnostic data to CSV.\nDefault: Off. Recommended: enable only for testing and troubleshooting.")]
    public bool EnableCsvLogging { get; set; } = false;

    private const int MaximumDifferenceOrder = 8;
    private const int MaximumOutlierBlock = 3;
    private const int MaximumOutlierPasses = 2;
    private const double RateTransitionFraction = 0.10;
    private const int RateEstimateMinimumSamples = 256;
    private const double RateEstimateMinimumSeconds = 0.35;
    private const int StableRateWindowsBeforeRebuild = 2;

    protected override void ConsumeState()
    {
        if (State is not ITabletReport report)
        {
            OnEmit();
            return;
        }
        lock (gate)
        {
            long now = Stopwatch.GetTimestamp();
            if (!rateInitialized)
            {
                estimatedRate = 700.0;
                rateInitialized = true;
            }

            bool discontinuity = IsStreamDiscontinuity(now);
            inputSequence++;
            if (discontinuity) ResetStream(report.Position, now);
            else UpdateRateEstimate(now);

            latestInputWallMs = TicksToMilliseconds(now);
            lastArrivalTicks = now;
            latestRaw = report.Position;
            haveRaw = true;
            if (!discontinuity) AddRawHistory(latestRaw);

            EnsureFilterNonBlocking();
            TryFinishModelBuild();

            if (filter is not null)
            {
                int lag = Math.Min(EffectiveLatencySamples(), Math.Max(0, filter.WindowSamples - 2));
                bool ready = activeSettings?.UseAdaptiveModes == true
                    ? filter.PushAdaptiveDelayAverage(report.Position, lag, out var result)
                    : filter.Push(report.Position, out result);
                if (ready) latestWindow = result;
            }

            UpdateLoggingState();
            if (latestWindow is not null) EnqueueInputLogRow(now, report.Position, latestWindow);
        }

        if (!EnableAsyncMode)
            EmitCurrentState();
    }

    protected override void UpdateState()
    {
        if (!EnableAsyncMode) return;
        EmitCurrentState();
    }

    private void EmitCurrentState()
    {
        if (State is not ITabletReport report) return;
        Vector2 output;
        lock (gate)
        {
            TryFinishModelBuild();
            long now = Stopwatch.GetTimestamp();

            // Some tablets leave the last ITabletReport cached when the pen exits
            // proximity. Do not depend on an out-of-range report to stop this timer:
            // otherwise the final reconstructed velocity is extrapolated forever.
            lastInputStalenessMs = lastArrivalTicks == 0 ? double.NaN : (now - lastArrivalTicks) * 1000.0 / Stopwatch.Frequency;
            if (lastArrivalTicks != 0 && lastInputStalenessMs >= EffectiveStalenessTimeoutMs())
            {
                lastEvalIndex = double.NaN;
                return;
            }

            if (lastOutputTicks != 0) lastOutputIntervalMs = (now - lastOutputTicks) * 1000.0 / Stopwatch.Frequency;
            lastOutputTicks = now;

            if (latestWindow is null || filter is null)
            {
                if (!haveRaw) return;
                output = latestRaw;
                lastEvalIndex = double.NaN;
            }
            else
            {
                // No PLL or virtual playback clock. The output timer directly asks for
                // the reconstructed trajectory at wall-clock (now - configured latency).
                double requestedWallMs = TicksToMilliseconds(now) - EffectiveLatencySamples() * 1000.0 / filter.SampleRateHz;
                double deltaFromNewestMs = requestedWallMs - latestInputWallMs;
                lastEvalIndex = (filter.WindowSamples - 1) + deltaFromNewestMs * filter.SampleRateHz / 1000.0;

                // Clamp only the historical side. If the playback head runs ahead
                // of the newest report, FilteredWindow performs linear extrapolation.
                lastEvalIndex = Math.Max(0.0, lastEvalIndex);
                output = latestWindow.At(lastEvalIndex);
            }

            latestOutput = output;
            haveOutput = true;
            outputSequence++;
            EnqueueOutputLogRow(now, output);
        }
        report.Position = output;
        OnEmit();
    }

    private bool IsStreamDiscontinuity(long now)
    {
        if (lastArrivalTicks == 0) return true;
        double gapMs = (now - lastArrivalTicks) * 1000.0 / Stopwatch.Frequency;
        double thresholdMs = Math.Max(25.0, 5000.0 / Math.Max(30.0, estimatedRate));
        return gapMs > thresholdMs;
    }

    private void ResetStream(Vector2 firstPosition, long now)
    {
        // A new in-range epoch has no meaningful trajectory connection to the old one.
        // Initialize history as a stationary pen at the first new coordinate.
        rawHistory.Clear();
        for (int i = 0; i < 256; i++)
            rawHistory.Enqueue(firstPosition);
        if (filter is not null) filter.ResetInput(Enumerable.Repeat(firstPosition, filter.WindowSamples));
        latestWindow = null;
        rateAnchorTicks = now;
        rateAnchorSequence = inputSequence;
        stableRateWindows = 0;
        rateErrorEma = 0;
        streamResets++;
    }

    private FilterSettings MakeSettings(double rate)
    {
        var settings = new FilterSettings
        {
            NominalSampleRateHz = rate,
            WindowMilliseconds = Math.Clamp(WindowMilliseconds, 10, 500),
            MaximumDifferenceOrder = MaximumDifferenceOrder,
            ModesRetained = 12,
            UseAdaptiveModes = UseAdaptiveModes,
            AdaptiveStrength = Math.Max(0, FilterStrength),
            OutlierThreshold = Math.Max(.5, OutlierThreshold),
            MaxOutlierBlock = MaximumOutlierBlock,
            MaxOutlierPasses = MaximumOutlierPasses
        };
        int n = settings.WindowSamples(rate);
        int lag = Math.Min(EffectiveLatencySamples(), Math.Max(0, n - 2));
        int modes = UseAdaptiveModes
            ? ComplexityModel.AdaptiveModesForLookaheadSamples(settings, lag, n)
            : ComplexityModel.ScaleModesForStrength(12, settings.AdaptiveStrength, n);
        return settings with { ModesRetained = modes };
    }

    private void EnsureFilterNonBlocking()
    {
        var settings = MakeSettings(estimatedRate);
        if (filter is null)
        {
            if (modelBuildTask is null) StartModelBuild(settings, estimatedRate);
            else requestedModelRate = estimatedRate;
            return;
        }
        bool settingsChanged = activeSettings is null || !EquivalentSettings(settings, activeSettings);
        double rel = Math.Abs(filter.SampleRateHz - estimatedRate) / Math.Max(1, filter.SampleRateHz);
        if (!settingsChanged && (rel <= RateTransitionFraction || stableRateWindows < StableRateWindowsBeforeRebuild)) return;
        requestedModelRate = estimatedRate;
        if (modelBuildTask is null) StartModelBuild(settings, estimatedRate);
    }

    private void StartModelBuild(FilterSettings settings, double rate)
    {
        requestedModelRate = rate;
        modelBuildTask = Task.Factory.StartNew(() =>
        {
            Thread.CurrentThread.Priority = ThreadPriority.BelowNormal;
            var sw = Stopwatch.StartNew();
            var builtFilter = new TrajectoryFilter(settings, rate);
            sw.Stop();
            return new ModelBuildResult(builtFilter, settings, rate, sw.Elapsed.TotalMilliseconds);
        }, CancellationToken.None, TaskCreationOptions.LongRunning, TaskScheduler.Default);
    }

    private void TryFinishModelBuild()
    {
        if (modelBuildTask is null || !modelBuildTask.IsCompleted) return;
        if (!modelBuildTask.IsCompletedSuccessfully)
        {
            if (modelBuildTask.Exception is not null)
                Log.Exception(modelBuildTask.Exception.GetBaseException());
            modelBuildTask = null;
            return;
        }

        var built = modelBuildTask.Result;
        modelBuildTask = null;
        lastModelBuildMs = built.Milliseconds;
        // A stale build is safe to use temporarily, but if rate moved substantially
        // schedule the newest model immediately after the atomic swap.
        built.Filter.ResetInput(rawHistory);
        filter = built.Filter;
        activeSettings = built.Settings;
        modelSwaps++;
        latestWindow = null; // next input computes a full window under the new basis.
        if (Math.Abs(requestedModelRate - built.Rate) / Math.Max(1, built.Rate) > RateTransitionFraction && stableRateWindows >= StableRateWindowsBeforeRebuild)
            StartModelBuild(MakeSettings(requestedModelRate), requestedModelRate);
    }

    private static bool EquivalentSettings(FilterSettings a, FilterSettings b) =>
        a.WindowMilliseconds == b.WindowMilliseconds && a.MaximumDifferenceOrder == b.MaximumDifferenceOrder &&
        a.ModesRetained == b.ModesRetained && a.UseAdaptiveModes == b.UseAdaptiveModes && a.AdaptiveStrength == b.AdaptiveStrength &&
        a.OutlierThreshold == b.OutlierThreshold &&
        a.MaxOutlierBlock == b.MaxOutlierBlock && a.MaxOutlierPasses == b.MaxOutlierPasses;

    private void UpdateRateEstimate(long now)
    {
        if (rateAnchorTicks == 0)
        {
            rateAnchorTicks = now;
            rateAnchorSequence = inputSequence;
            return;
        }
        long count = inputSequence - rateAnchorSequence;
        double elapsed = (now - rateAnchorTicks) / (double)Stopwatch.Frequency;
        if (count < RateEstimateMinimumSamples || elapsed < RateEstimateMinimumSeconds) return;
        double observed = count / elapsed;
        if (observed is < 30 or > 4000)
        {
            rateAnchorTicks = now;
            rateAnchorSequence = inputSequence;
            return;
        }
        double rel = (observed - estimatedRate) / estimatedRate;
        rateErrorEma = .8 * rateErrorEma + .2 * rel;
        // Long observation windows reject scheduler jitter. Only sustained ~10% changes
        // acquire quickly; smaller fluctuations have little value for model sizing.
        double alpha = Math.Abs(rateErrorEma) > RateTransitionFraction ? .50 : .05;
        estimatedRate = Math.Clamp(estimatedRate + alpha * (observed - estimatedRate), 30, 4000);
        stableRateWindows = Math.Abs(rel) < RateTransitionFraction ? stableRateWindows + 1 : 0;
        rateAnchorTicks = now;
        rateAnchorSequence = inputSequence;
    }

    private static double TicksToMilliseconds(long ticks) => ticks * 1000.0 / Stopwatch.Frequency;
    private int EffectiveLatencySamples() => Math.Clamp(LatencySamples, 0, 20);
    private double EffectiveStalenessTimeoutMs() => Math.Clamp(StalenessTimeoutMs, 1.0, 1000.0);

    private void AddRawHistory(Vector2 point)
    {
        rawHistory.Enqueue(point);
        while (rawHistory.Count > 256)
            rawHistory.Dequeue();
    }

    private void UpdateLoggingState()
    {
        if (!EnableCsvLogging)
        {
            logCts?.Cancel();
            logCts = null;
            logTask = null;
            loggingFailed = false;
            return;
        }

        if (logCts is not null || loggingFailed) return;

        try
        {
            string root = Path.Combine(GetOpenTabletDriverDataDirectory(), "SpectralReconstructorLogs");
            Directory.CreateDirectory(root);
            string path = Path.Combine(root, $"spectral-reconstructor-{DateTime.Now:yyyyMMdd-HHmmss}.csv");
            logCts = new CancellationTokenSource();
            var token = logCts.Token;
            logTask = Task.Factory.StartNew(
                () => LogWriterLoop(path, token),
                token,
                TaskCreationOptions.LongRunning,
                TaskScheduler.Default);
        }
        catch (Exception ex)
        {
            loggingFailed = true;
            Log.Exception(ex);
        }
    }

    private static string GetOpenTabletDriverDataDirectory()
    {
        if (!OperatingSystem.IsLinux())
            return Path.Combine(
                Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData),
                "OpenTabletDriver");

        string? configRoot = Environment.GetEnvironmentVariable("XDG_CONFIG_HOME");
        if (string.IsNullOrWhiteSpace(configRoot) || !Path.IsPathRooted(configRoot))
        {
            configRoot = Path.Combine(
                Environment.GetFolderPath(Environment.SpecialFolder.UserProfile),
                ".config");
        }

        return Path.Combine(configRoot, "OpenTabletDriver");
    }

    private void LogWriterLoop(string path, CancellationToken token)
    {
        try
        {
            Thread.CurrentThread.Priority = ThreadPriority.BelowNormal;
            using var w = new StreamWriter(path, false, new UTF8Encoding(false), 1 << 16);
            w.WriteLine("event_type,event_sequence,time_ms,raw_x,raw_y,output_x,output_y,output_minus_raw_x,output_minus_raw_y,window_had_replacement,max_outlier_score,complexity,estimated_rate_hz,model_rate_hz,window_samples,input_sequence,output_sequence,latest_input_wall_ms,evaluation_index,input_staleness_ms,stream_resets,last_output_interval_ms,model_build_active,model_build_ms,model_swaps,dropped_log_rows");
            while (!token.IsCancellationRequested)
            {
                while (logRows.TryDequeue(out var row))
                    w.WriteLine(FormatLogRow(row));
                w.Flush();
                Thread.Sleep(100);
            }
            while (logRows.TryDequeue(out var row))
                w.WriteLine(FormatLogRow(row));
            w.Flush();
        }
        catch (Exception ex)
        {
            loggingFailed = true;
            Log.Exception(ex);
        }
    }

    private void EnqueueInputLogRow(long now, Vector2 raw, FilteredWindow result)
    {
        if (logCts is null) return;
        if (logRows.Count > 4096)
        {
            droppedLogRows++;
            return;
        }
        double maxScore = result.Replacements.Count == 0 ? 0 : result.Replacements.Max(r => r.RelativeResidualScore);
        Vector2 o = haveOutput ? latestOutput : raw;
        logRows.Enqueue(CreateLogRow("input", inputSequence, now, raw, o,
            result.Replacements.Count != 0, maxScore, result.Complexity));
    }

    private void EnqueueOutputLogRow(long now, Vector2 output)
    {
        if (logCts is null) return;
        if (logRows.Count > 4096)
        {
            droppedLogRows++;
            return;
        }
        logRows.Enqueue(CreateLogRow("output", outputSequence, now, latestRaw, output,
            null, double.NaN, latestWindow?.Complexity ?? double.NaN));
    }

    private LogRow CreateLogRow(
        string eventType,
        long eventSequence,
        long now,
        Vector2 raw,
        Vector2 output,
        bool? windowHadReplacement,
        double maxOutlierScore,
        double complexity) =>
        new(
            eventType,
            eventSequence,
            now * 1000.0 / Stopwatch.Frequency,
            raw,
            output,
            windowHadReplacement,
            maxOutlierScore,
            complexity,
            estimatedRate,
            filter?.SampleRateHz ?? double.NaN,
            filter?.WindowSamples ?? 0,
            inputSequence,
            outputSequence,
            latestInputWallMs,
            lastEvalIndex,
            lastInputStalenessMs,
            streamResets,
            lastOutputIntervalMs,
            modelBuildTask is not null && !modelBuildTask.IsCompleted,
            lastModelBuildMs,
            modelSwaps,
            droppedLogRows);

    private static string FormatLogRow(LogRow row)
    {
        var c = CultureInfo.InvariantCulture;
        string F(double value) => double.IsNaN(value) ? "NaN" : value.ToString("R", c);
        string B(bool? value) => value.HasValue ? (value.Value ? "1" : "0") : "";
        var fields = new[]
        {
            row.EventType,
            row.EventSequence.ToString(c),
            F(row.TimeMs),
            F(row.Raw.X),
            F(row.Raw.Y),
            F(row.Output.X),
            F(row.Output.Y),
            F(row.Output.X - row.Raw.X),
            F(row.Output.Y - row.Raw.Y),
            B(row.WindowHadReplacement),
            F(row.MaxOutlierScore),
            F(row.Complexity),
            F(row.EstimatedRate),
            F(row.ModelRate),
            row.WindowSamples.ToString(c),
            row.InputSequence.ToString(c),
            row.OutputSequence.ToString(c),
            F(row.LatestInputWallMs),
            F(row.EvaluationIndex),
            F(row.InputStalenessMs),
            row.StreamResets.ToString(c),
            F(row.OutputIntervalMs),
            row.ModelBuildActive ? "1" : "0",
            F(row.ModelBuildMs),
            row.ModelSwaps.ToString(c),
            row.DroppedLogRows.ToString(c)
        };
        return string.Join(',', fields);
    }

    private readonly record struct LogRow(
        string EventType,
        long EventSequence,
        double TimeMs,
        Vector2 Raw,
        Vector2 Output,
        bool? WindowHadReplacement,
        double MaxOutlierScore,
        double Complexity,
        double EstimatedRate,
        double ModelRate,
        int WindowSamples,
        long InputSequence,
        long OutputSequence,
        double LatestInputWallMs,
        double EvaluationIndex,
        double InputStalenessMs,
        long StreamResets,
        double OutputIntervalMs,
        bool ModelBuildActive,
        double ModelBuildMs,
        int ModelSwaps,
        int DroppedLogRows);

    private sealed record ModelBuildResult(TrajectoryFilter Filter, FilterSettings Settings, double Rate, double Milliseconds);
}
