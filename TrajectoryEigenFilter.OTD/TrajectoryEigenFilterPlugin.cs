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

[PluginName("AA The Perfect Filter")]
public sealed class TrajectoryEigenFilterPlugin : AsyncPositionedPipelineElement<IDeviceReport>
{
    private const float OutputFrequencyHz = 1000f; // OTD native timer path.
    private readonly object gate = new();
    private TrajectoryFilter? filter;
    private FilterSettings? activeSettings;
    private readonly Queue<Vector2> rawHistory = new();

    private long inputSequence, outputSequence;
    private long lastArrivalTicks, lastOutputTicks, rateAnchorTicks, rateAnchorSequence;
    private double estimatedRate = 700.0, lastOutputIntervalMs = double.NaN;
    private bool rateInitialized, schedulerInitialized;
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

    private readonly ConcurrentQueue<string> logLines = new();
    private CancellationTokenSource? logCts;
    private Task? logTask;
    private volatile bool loggingFailed;
    private int droppedLogRows;

    public override PipelinePosition Position => PipelinePosition.PreTransform;

    // OTD scheduler is always 1000 Hz. This remains an implementation detail;
    // EnsureScheduler1000Hz applies it after OTD initializes the base timer.
    [DefaultPropertyValue(1000.0)]
    public new float Frequency
    {
        get => base.Frequency;
        set => base.Frequency = OutputFrequencyHz;
    }

    [Property("Window duration"), Unit("ms"), DefaultPropertyValue(100.0), ToolTip("Trajectory history used by the filter.\nDefault: 100 ms. Recommended: 70-150 ms.\n100 ms was used for current tuning; larger windows cost more model-build CPU and memory.")]
    public double WindowMilliseconds { get; set; } = 100.0;

    [Property("Filter strength"), DefaultPropertyValue(1.0), ToolTip("Inverse retained-mode control. Strength 1 uses 12 modes at zero latency; 0 retains all modes. Higher values retain fewer modes.\nDefault: 1. Recommended: 1.")]
    public double FilterStrength { get; set; } = 1.0;

    [BooleanProperty("Use adaptive modes", "Varies the retained mode count with reconstruction look-ahead."), DefaultPropertyValue(true), ToolTip("Uses the sample-look-ahead schedule 12, 9, 8, 7, 6, then 5, scaled by Filter strength.\nDefault: On. Recommended: On.")]
    public bool UseAdaptiveModes { get; set; } = true;

    [Property("Outlier threshold"), DefaultPropertyValue(6.0), ToolTip("Robust residual threshold for isolated bad samples.\nDefault: 6. Recommended: 5-8.\nLower values reject more aggressively; higher values reserve replacement for more extreme excursions.")]
    public double OutlierThreshold { get; set; } = 6.0;

    [Property("Latency"), Unit("ms"), DefaultPropertyValue(5.0), ToolTip("Look-ahead/buffer before output.\nDefault: 5 ms. Recommended: 2-20 ms.\nZero latency is fully functional but not recommended: noise rejection and input-jitter tolerance are reduced. Linear extrapolation is used if the output timer runs ahead of the latest input report.")]
    public double LatencyMs { get; set; } = 5.0;

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
            EnsureScheduler1000Hz();
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

            if (filter is not null && filter.Push(report.Position, out var result))
                latestWindow = result;

            UpdateLoggingState();
            if (latestWindow is not null) EnqueueInputLogRow(now, report.Position, latestWindow);
        }
        // Positioned async elements emit from UpdateState at the fixed OTD timer rate.
    }

    protected override void UpdateState()
    {
        if (State is not ITabletReport report) return;
        Vector2 output;
        lock (gate)
        {
            EnsureScheduler1000Hz();
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
                double requestedWallMs = TicksToMilliseconds(now) - EffectiveLatencyMs();
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
        }
        report.Position = output;
        OnEmit();
    }

    private void EnsureScheduler1000Hz()
    {
        if (schedulerInitialized) return;
        try
        {
            base.Frequency = OutputFrequencyHz;
            schedulerInitialized = true;
        }
        catch (NullReferenceException) { /* OTD initializes the timer after construction. */ }
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
        int lag = (int)Math.Round(EffectiveLatencyMs() * rate / 1000.0);
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
    private double EffectiveLatencyMs() => Math.Clamp(LatencyMs, 0.0, 20.0);
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

        string root = Path.Combine(
            Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData),
            "OpenTabletDriver",
            "TrajectoryEigenFilterLogs");
        Directory.CreateDirectory(root);
        string path = Path.Combine(root, $"trajectory-eigen-{DateTime.Now:yyyyMMdd-HHmmss}.csv");
        logCts = new CancellationTokenSource();
        var token = logCts.Token;
        logTask = Task.Factory.StartNew(
            () => LogWriterLoop(path, token),
            token,
            TaskCreationOptions.LongRunning,
            TaskScheduler.Default);
    }

    private void LogWriterLoop(string path, CancellationToken token)
    {
        try
        {
            Thread.CurrentThread.Priority = ThreadPriority.BelowNormal;
            using var w = new StreamWriter(path, false, new UTF8Encoding(false), 1 << 16);
            w.WriteLine("sample,time_ms,raw_x,raw_y,output_x,output_y,output_minus_raw_x,output_minus_raw_y,window_had_replacement,max_outlier_score,complexity,estimated_rate_hz,model_rate_hz,window_samples,input_sequence,output_sequence,latest_input_wall_ms,evaluation_index,input_staleness_ms,stream_resets,last_output_interval_ms,model_build_active,model_build_ms,model_swaps,dropped_log_rows");
            while (!token.IsCancellationRequested)
            {
                while (logLines.TryDequeue(out var line))
                    w.WriteLine(line);
                w.Flush();
                Thread.Sleep(100);
            }
            while (logLines.TryDequeue(out var line))
                w.WriteLine(line);
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
        if (logLines.Count > 4096)
        {
            droppedLogRows++;
            return;
        }
        double maxScore = result.Replacements.Count == 0 ? 0 : result.Replacements.Max(r => r.RelativeResidualScore);
        var c = CultureInfo.InvariantCulture;
        string F(double v) => double.IsNaN(v) ? "NaN" : v.ToString("R", c);
        double timeMs = now * 1000.0 / Stopwatch.Frequency;
        Vector2 o = haveOutput ? latestOutput : raw;
        var fields = new[] { inputSequence.ToString(c), F(timeMs), F(raw.X), F(raw.Y), F(o.X), F(o.Y), F(o.X - raw.X), F(o.Y - raw.Y), result.Replacements.Count != 0 ? "1" : "0", F(maxScore), F(result.Complexity), F(estimatedRate), F(filter?.SampleRateHz ?? double.NaN), (filter?.WindowSamples ?? 0).ToString(c), inputSequence.ToString(c), outputSequence.ToString(c), F(latestInputWallMs), F(lastEvalIndex), F(lastInputStalenessMs), streamResets.ToString(c), F(lastOutputIntervalMs), modelBuildTask is not null && !modelBuildTask.IsCompleted ? "1" : "0", F(lastModelBuildMs), modelSwaps.ToString(c), droppedLogRows.ToString(c) };
        logLines.Enqueue(string.Join(',', fields));
    }

    private sealed record ModelBuildResult(TrajectoryFilter Filter, FilterSettings Settings, double Rate, double Milliseconds);
}
