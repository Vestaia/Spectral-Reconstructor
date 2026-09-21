namespace TrajectoryEigenFilterCore;

public sealed class ComplexityModel
{
    public int N { get; }
    public double SampleRateHz { get; }
    public double[,] L { get; }
    public double[] Eigenvalues { get; }
    public double[,] Q { get; }
    readonly FilterSettings settings;
    readonly double[,] reconstruction;

    public ComplexityModel(FilterSettings settings, double sampleRateHz)
    {
        this.settings = settings; SampleRateHz = sampleRateHz; N = settings.WindowSamples(sampleRateHz);
        L = BuildOperator(settings, sampleRateHz, N);

        var (ev, q) = SymmetricEigen.Decompose(L);
        Eigenvalues = ev;
        Q = q;

        reconstruction = BuildReconstruction(Q, ClampModes(settings.ModesRetained));
    }

    static double[,] BuildReconstruction(double[,] q, int modesRetained)
    {
        int n = q.GetLength(0); var r = new double[n, n];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
        {
            double sum = 0;
            for (int k = 0; k < modesRetained; k++) sum += q[i, k] * q[j, k];
            r[i, j] = sum;
        }
        return r;
    }

    static double[,] BuildOperator(FilterSettings s, double fs, int n)
    {
        int hi = s.MaximumDifferenceOrder <= 0 ? n - 1 : Math.Min(n - 1, s.MaximumDifferenceOrder);
        var l = new double[n, n];

        // Every constituent operator is built from unit-L2 rows and then
        // normalized to unit spectral radius (largest eigenvalue = 1).
        // Order 1 is deliberately absent, so every term above identity
        // annihilates affine trajectories.
        for (int order = 0; order <= hi; order++)
        {
            if (order == 1) continue;
            int width = Math.Min(n, order + 2);
            if (width <= order) continue;
            var orderL = new double[n, n];
            for (int eval = 0; eval < n; eval++)
            {
                int stencilStart = StencilStart(eval, width, n);
                var nodes = ConsecutiveNodes(stencilStart, width);
                var row = FiniteDifferenceWeights(nodes, eval, order);
                AddNormalizedRowPenalty(orderL, row, stencilStart);
            }
            NormalizeOperatorSpectral(orderL);
            AddOperator(l, orderL);
        }

        return l;
    }

    static void NormalizeOperatorSpectral(double[,] matrix)
    {
        var (values, _) = SymmetricEigen.Decompose(matrix);
        if (values.Length == 0) return;
        double maxEigenvalue = values[^1];
        if (!(maxEigenvalue > 1e-14) || !double.IsFinite(maxEigenvalue)) return;
        double inv = 1.0 / maxEigenvalue;
        int n = matrix.GetLength(0);
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) matrix[i, j] *= inv;
    }

    static void AddOperator(double[,] destination, double[,] source)
    {
        int n = destination.GetLength(0);
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) destination[i, j] += source[i, j];
    }

    static void AddNormalizedRowPenalty(double[,] matrix, double[] row, int start)
    {
        double norm = Math.Sqrt(row.Sum(v => v * v));
        if (!(norm > 1e-14) || !double.IsFinite(norm)) return;
        for (int a = 0; a < row.Length; a++)
        {
            double va = row[a] / norm;
            for (int b = 0; b < row.Length; b++) matrix[start + a, start + b] += va * (row[b] / norm);
        }
    }

    static int StencilStart(int eval, int width, int n)
    {
        // Nearest available support around the derivative evaluation point.
        // Clamp at the observation boundary rather than assuming continuation.
        int left = (width - 1) / 2;
        return Math.Clamp(eval - left, 0, n - width);
    }

    static double[] ConsecutiveNodes(int start, int width)
    { var x = new double[width]; for (int i = 0; i < width; i++) x[i] = start + i; return x; }

    // Fornberg finite-difference weights, determined by polynomial exactness on
    // the available sample locations.
    static double[] FiniteDifferenceWeights(ReadOnlySpan<double> nodes, double x0, int derivative)
    {
        int n = nodes.Length;
        if (derivative < 0 || derivative >= n) throw new ArgumentOutOfRangeException(nameof(derivative));
        var c = new double[n, derivative + 1];
        c[0, 0] = 1.0;
        double c1 = 1.0, c4 = nodes[0] - x0;
        for (int i = 1; i < n; i++)
        {
            int mn = Math.Min(i, derivative);
            double c2 = 1.0, c5 = c4; c4 = nodes[i] - x0;
            for (int j = 0; j < i; j++)
            {
                double c3 = nodes[i] - nodes[j];
                if (c3 == 0) throw new ArgumentException("Finite-difference nodes must be distinct.", nameof(nodes));
                c2 *= c3;
                if (j == i - 1)
                    for (int k = mn; k >= 1; k--)
                        c[i, k] = c1 * (k * c[i - 1, k - 1] - c5 * c[i - 1, k]) / c2;
                if (j == i - 1) c[i, 0] = -c1 * c5 * c[i - 1, 0] / c2;
                for (int k = mn; k >= 1; k--)
                    c[j, k] = (c4 * c[j, k] - k * c[j, k - 1]) / c3;
                c[j, 0] = c4 * c[j, 0] / c3;
            }
            c1 = c2;
        }
        var w = new double[n]; for (int i = 0; i < n; i++) w[i] = c[i, derivative]; return w;
    }

    public void MultiplyL(ReadOnlySpan<double> x, Span<double> y)
    { for (int i = 0; i < N; i++) { double s = 0; for (int j = 0; j < N; j++) s += L[i, j] * x[j]; y[i] = s; } }

    public double Complexity(ReadOnlySpan<double> x)
    { double c = 0; for (int i = 0; i < N; i++) { double r = 0; for (int j = 0; j < N; j++) r += L[i, j] * x[j]; c += x[i] * r; } return c; }

    public double Gain(int k)
    {
        if ((uint)k >= (uint)N) throw new ArgumentOutOfRangeException(nameof(k));
        return k < ClampModes(settings.ModesRetained) ? 1.0 : 0.0;
    }

    // Hard-rank kernel used by the latency-dependent ensemble. Building it
    // is model/configuration work; runtime is only one dot product per axis.
    public (double[] Kernel, double NoiseTransmission) BuildModeCountKernel(int index, int modesRetained)
    {
        if ((uint)index >= (uint)N) throw new ArgumentOutOfRangeException(nameof(index));
        modesRetained = ClampModes(modesRetained);
        var h = new double[N]; double noise = 0;
        for (int j = 0; j < N; j++)
        {
            double v = 0;
            for (int k = 0; k < modesRetained; k++) v += Q[index, k] * Q[j, k];
            h[j] = v; noise += v * v;
        }
        return (h, noise);
    }

    public int AdaptiveModesForLookaheadSamples(int lookaheadSamples)
    {
        if (!settings.UseAdaptiveModes) return ClampModes(settings.ModesRetained);
        return AdaptiveModesForLookaheadSamples(settings, lookaheadSamples, N);
    }

    public static int AdaptiveModesForLookaheadSamples(FilterSettings settings, int lookaheadSamples, int windowSamples)
    {
        if (!settings.UseAdaptiveModes) return ClampModes(settings.ModesRetained, windowSamples);
        int baseModes = Math.Max(0, lookaheadSamples) switch { 0 => 12, 1 => 9, 2 => 8, 3 => 7, 4 => 6, _ => 5 };
        return ScaleModesForStrength(baseModes, settings.AdaptiveStrength, windowSamples);
    }

    public static int ScaleModesForStrength(int baseModes, double strength, int windowSamples)
    {
        strength = Math.Max(0, strength);
        if (strength == 0) return windowSamples;
        double scaled = 2.0 + (baseModes - 2.0) / (strength * strength);
        return ClampModes((int)Math.Round(scaled, MidpointRounding.AwayFromZero), windowSamples);
    }

    int ClampModes(int modesRetained) => ClampModes(modesRetained, N);
    static int ClampModes(int modesRetained, int n) => Math.Clamp(modesRetained, Math.Min(2, n), n);

    public double SmoothAt(ReadOnlySpan<double> x, int index)
    {
        if ((uint)index >= (uint)N) throw new ArgumentOutOfRangeException(nameof(index));
        double y = 0; for (int j = 0; j < N; j++) y += reconstruction[index, j] * x[j];
        return y;
    }

    public double[] Smooth(ReadOnlySpan<double> x)
    {
        var y = new double[N];
        for (int i = 0; i < N; i++)
        {
            double sum = 0;
            for (int j = 0; j < N; j++) sum += reconstruction[i, j] * x[j];
            y[i] = sum;
        }
        return y;
    }

    public NumericalDiagnostics GetNumericalDiagnostics()
    {
        var entries = new List<double>(N * N);
        double maxL = double.NegativeInfinity;
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++)
        {
            double v = L[i, j];
            maxL = Math.Max(maxL, v);
            if (v != 0.0) entries.Add(Math.Abs(v));
        }
        entries.Sort();
        double medianNonZeroAbsL = entries.Count == 0 ? 0.0 :
            (entries.Count % 2 != 0 ? entries[entries.Count / 2] : 0.5 * (entries[entries.Count / 2 - 1] + entries[entries.Count / 2]));
        return new(maxL, medianNonZeroAbsL);
    }

    public AffineEigenspaceDiagnostics GetAffineEigenspaceDiagnostics()
    {
        var constant = new double[N];
        var linear = new double[N];
        double center = (N - 1) * 0.5;
        for (int i = 0; i < N; i++) { constant[i] = 1.0; linear[i] = i - center; }
        Normalize(constant); Normalize(linear);
        return new(Eigenvalues[0], Eigenvalues[1],
            ProjectionResidual(constant), ProjectionResidual(linear),
            OperatorEigenResidual(constant, Eigenvalues[0]), OperatorEigenResidual(linear, Eigenvalues[0]));
    }

    double ProjectionResidual(double[] vector)
    {
        var residual = (double[])vector.Clone();
        for (int k = 0; k < 2; k++)
        {
            double coefficient = 0;
            for (int i = 0; i < N; i++) coefficient += Q[i, k] * vector[i];
            for (int i = 0; i < N; i++) residual[i] -= coefficient * Q[i, k];
        }
        return Math.Sqrt(residual.Sum(v => v * v));
    }

    double OperatorEigenResidual(double[] vector, double eigenvalue)
    {
        var product = new double[N]; MultiplyL(vector, product);
        double sum = 0; for (int i = 0; i < N; i++) { double r = product[i] - eigenvalue * vector[i]; sum += r * r; }
        return Math.Sqrt(sum);
    }

    static void Normalize(double[] vector)
    {
        double norm = Math.Sqrt(vector.Sum(v => v * v));
        for (int i = 0; i < vector.Length; i++) vector[i] /= norm;
    }

}

public sealed record NumericalDiagnostics(double MaxOperatorEntry, double MedianNonZeroAbsOperatorEntry);
public sealed record AffineEigenspaceDiagnostics(double FirstEigenvalue, double SecondEigenvalue, double ConstantProjectionResidual, double LinearProjectionResidual, double ConstantOperatorResidual, double LinearOperatorResidual);
