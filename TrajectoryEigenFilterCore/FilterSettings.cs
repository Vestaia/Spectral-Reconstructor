namespace TrajectoryEigenFilterCore;

public sealed record FilterSettings
{
    public double NominalSampleRateHz { get; init; } = 700.0;
    public double WindowMilliseconds { get; init; } = 100.0;
    public int MaximumDifferenceOrder { get; init; } = 0;

    // Provisional hard spectral ranks. Eigenvalues define the ordering only.
    public int ModesRetained { get; init; } = 12;
    public bool UseAdaptiveModes { get; init; } = true;
    public double AdaptiveStrength { get; init; } = 1.0;

    // Relative robust residual score. A candidate is anomalous when its
    // model-consistent residual magnitude exceeds median + threshold*MADsigma.
    public double OutlierThreshold { get; init; } = 6.0;
    public int MaxOutlierBlock { get; init; } = 3;
    public int MaxOutlierPasses { get; init; } = 2;

    public int WindowSamples(double sampleRateHz) => Math.Max(2, (int)Math.Round(sampleRateHz * WindowMilliseconds / 1000.0));
}
