namespace TrajectoryEigenFilterCore;

public sealed record FilterSettings
{
    public double NominalSampleRateHz { get; init; } = 700.0;
    public double WindowMilliseconds { get; init; } = 100.0;
    public int MinimumDifferenceOrder { get; init; } = 0;
    public int MaximumDifferenceOrder { get; init; } = 0;

    // Exact hard cutoff in eigenvalue space. Modes at or below the cutoff pass.
    public double LambdaCutoff { get; init; } = 1.5;

    public double LocalDifferenceStrength { get; init; } = 1.0;

    // Experimental latency-dependent hard cutoffs used when averaging repeated
    // estimates of the same physical sample.  Anchors are interpolated in ms.
    public double AdaptiveLambdaAt0Ms { get; init; } = 1.50;
    public double AdaptiveLambdaAt2Ms { get; init; } = 1.20;
    public bool UseAdaptiveLambda { get; init; } = true;
    public double AdaptiveLambdaAt5Ms { get; init; } = 1.05;
    public double AdaptiveLambdaAt10Ms { get; init; } = 1.04;
    public double AdaptiveLambdaAt20Ms { get; init; } = 1.03;

    // Relative robust residual score. A candidate is anomalous when its
    // model-consistent residual magnitude exceeds median + threshold*MADsigma.
    public double OutlierThreshold { get; init; } = 6.0;
    public int MaxOutlierBlock { get; init; } = 3;
    public int MaxOutlierPasses { get; init; } = 2;

    public int WindowSamples(double sampleRateHz)=>Math.Max(MinimumDifferenceOrder+2,(int)Math.Round(sampleRateHz*WindowMilliseconds/1000.0));
}
