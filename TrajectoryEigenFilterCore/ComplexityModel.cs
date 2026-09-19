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
        this.settings=settings; SampleRateHz=sampleRateHz; N=settings.WindowSamples(sampleRateHz);
        L=BuildOperator(settings,sampleRateHz,N);

        var (ev,q)=SymmetricEigen.Decompose(L);
        Eigenvalues=ev;
        Q=q;

        reconstruction=BuildReconstruction(Q,Eigenvalues,settings);
    }

    static double[,] BuildReconstruction(double[,] q,double[] eigenvalues,FilterSettings settings)
    {
        int n=eigenvalues.Length; var r=new double[n,n];
        for(int i=0;i<n;i++) for(int j=0;j<n;j++)
        {
            double sum=0;
            for(int k=0;k<n;k++) sum+=q[i,k]*GainForLambda(eigenvalues[k],settings)*q[j,k];
            r[i,j]=sum;
        }
        return r;
    }

    static double GainForLambda(double lambda,FilterSettings settings)
        => lambda<=settings.LambdaCutoff ? 1.0 : 0.0;

    static double[,] BuildOperator(FilterSettings s,double fs,int n)
    {
        int hi=s.MaximumDifferenceOrder<=0 ? n-1 : Math.Min(n-1,s.MaximumDifferenceOrder);
        int lo=Math.Clamp(s.MinimumDifferenceOrder,0,hi);
        var l=new double[n,n];

        // Every constituent operator is built from unit-L2 rows and then
        // normalized to unit spectral radius (largest eigenvalue = 1).
        // This removes all reference-sinusoid a_n weighting and gives local
        // derivative orders and non-local scales a common eigenvalue scale.
        if(s.LocalDifferenceStrength>0)
        {
            for(int order=lo;order<=hi;order++)
            {
                int width=Math.Min(n,order+2);
                if(width<=order) continue;
                var orderL=new double[n,n];
                for(int eval=0;eval<n;eval++)
                {
                    int stencilStart=StencilStart(eval,width,n);
                    var nodes=ConsecutiveNodes(stencilStart,width);
                    var row=FiniteDifferenceWeights(nodes,eval,order);
                    AddNormalizedRowPenalty(orderL,row,stencilStart);
                }
                NormalizeOperatorSpectral(orderL);
                AddOperator(l,orderL,s.LocalDifferenceStrength);
            }
        }

        return l;
    }

    static void NormalizeOperatorSpectral(double[,] matrix)
    {
        var (values,_)=SymmetricEigen.Decompose(matrix);
        if(values.Length==0) return;
        double maxEigenvalue=values[^1];
        if(!(maxEigenvalue>1e-14) || !double.IsFinite(maxEigenvalue)) return;
        double inv=1.0/maxEigenvalue;
        int n=matrix.GetLength(0);
        for(int i=0;i<n;i++) for(int j=0;j<n;j++) matrix[i,j]*=inv;
    }

    static void AddOperator(double[,] destination,double[,] source,double strength)
    {
        int n=destination.GetLength(0);
        for(int i=0;i<n;i++) for(int j=0;j<n;j++) destination[i,j]+=strength*source[i,j];
    }

    static void AddNormalizedRowPenalty(double[,] matrix,double[] row,int start)
    {
        double norm=Math.Sqrt(row.Sum(v=>v*v));
        if(!(norm>1e-14) || !double.IsFinite(norm)) return;
        for(int a=0;a<row.Length;a++)
        {
            double va=row[a]/norm;
            for(int b=0;b<row.Length;b++) matrix[start+a,start+b]+=va*(row[b]/norm);
        }
    }

    static int StencilStart(int eval,int width,int n)
    {
        // Nearest available support around the derivative evaluation point.
        // Clamp at the observation boundary rather than assuming continuation.
        int left=(width-1)/2;
        return Math.Clamp(eval-left,0,n-width);
    }

    static double[] ConsecutiveNodes(int start,int width)
    { var x=new double[width];for(int i=0;i<width;i++)x[i]=start+i;return x; }

    // Fornberg finite-difference weights.  This is Option A: coefficients are
    // determined only by polynomial exactness on the available sample locations;
    // there is no noise minimization or smoothing inside the derivative operator.
    static double[] FiniteDifferenceWeights(ReadOnlySpan<double> nodes,double x0,int derivative)
    {
        int n=nodes.Length;
        if(derivative<0 || derivative>=n) throw new ArgumentOutOfRangeException(nameof(derivative));
        var c=new double[n,derivative+1];
        c[0,0]=1.0;
        double c1=1.0,c4=nodes[0]-x0;
        for(int i=1;i<n;i++)
        {
            int mn=Math.Min(i,derivative);
            double c2=1.0,c5=c4;c4=nodes[i]-x0;
            for(int j=0;j<i;j++)
            {
                double c3=nodes[i]-nodes[j];
                if(c3==0) throw new ArgumentException("Finite-difference nodes must be distinct.",nameof(nodes));
                c2*=c3;
                if(j==i-1)
                    for(int k=mn;k>=1;k--)
                        c[i,k]=c1*(k*c[i-1,k-1]-c5*c[i-1,k])/c2;
                if(j==i-1)c[i,0]=-c1*c5*c[i-1,0]/c2;
                for(int k=mn;k>=1;k--)
                    c[j,k]=(c4*c[j,k]-k*c[j,k-1])/c3;
                c[j,0]=c4*c[j,0]/c3;
            }
            c1=c2;
        }
        var w=new double[n];for(int i=0;i<n;i++)w[i]=c[i,derivative];return w;
    }

    public void MultiplyL(ReadOnlySpan<double> x,Span<double> y)
    { for(int i=0;i<N;i++){double s=0;for(int j=0;j<N;j++)s+=L[i,j]*x[j];y[i]=s;} }

    public double Complexity(ReadOnlySpan<double> x)
    { double c=0;for(int i=0;i<N;i++){double r=0;for(int j=0;j<N;j++)r+=L[i,j]*x[j];c+=x[i]*r;}return c; }

    public double Gain(int k)
    {
        if ((uint)k >= (uint)N) throw new ArgumentOutOfRangeException(nameof(k));
        return GainForLambda(Eigenvalues[k],settings);
    }

    // Hard-cutoff kernel used by the latency-dependent ensemble.  Building it
    // is model/configuration work; runtime is only one dot product per axis.
    public (double[] Kernel,double NoiseTransmission) BuildHardCutoffKernel(int index,double lambdaCutoff)
    {
        if((uint)index>=(uint)N)throw new ArgumentOutOfRangeException(nameof(index));
        var h=new double[N];double noise=0;
        for(int j=0;j<N;j++)
        {
            double v=0;
            for(int k=0;k<N;k++) if(Eigenvalues[k]<=lambdaCutoff) v+=Q[index,k]*Q[j,k];
            h[j]=v;noise+=v*v;
        }
        return (h,noise);
    }

    public double AdaptiveLambdaForLookaheadSamples(int lookaheadSamples)
    {
        double ms=Math.Max(0,lookaheadSamples)*1000.0/SampleRateHz;
        if(!settings.UseAdaptiveLambda) return settings.LambdaCutoff;
        if(ms<=2.0)return Lerp(settings.AdaptiveLambdaAt0Ms,settings.AdaptiveLambdaAt2Ms,ms/2.0);
        if(ms<=5.0)return Lerp(settings.AdaptiveLambdaAt2Ms,settings.AdaptiveLambdaAt5Ms,(ms-2.0)/3.0);
        if(ms<=10.0)return Lerp(settings.AdaptiveLambdaAt5Ms,settings.AdaptiveLambdaAt10Ms,(ms-5.0)/5.0);
        if(ms<=20.0)return Lerp(settings.AdaptiveLambdaAt10Ms,settings.AdaptiveLambdaAt20Ms,(ms-10.0)/10.0);
        return settings.AdaptiveLambdaAt20Ms;
    }

    static double Lerp(double a,double b,double t)=>a+(b-a)*t;

    public double SmoothAt(ReadOnlySpan<double> x,int index)
    {
        if((uint)index>=(uint)N)throw new ArgumentOutOfRangeException(nameof(index));
        double y=0; for(int j=0;j<N;j++) y+=reconstruction[index,j]*x[j];
        return y;
    }

    // Deliberately slow reference path used by regression tests/diagnostics.
    public double SmoothAtReference(ReadOnlySpan<double> x,int index)
    {
        double y=0;
        for(int k=0;k<N;k++){double ck=0;for(int i=0;i<N;i++)ck+=Q[i,k]*x[i];y+=Q[index,k]*Gain(k)*ck;}
        return y;
    }

    public ReconstructionKernelDiagnostics GetKernelDiagnostics(int index)
    {
        if((uint)index>=(uint)N)throw new ArgumentOutOfRangeException(nameof(index));
        double sum=0,m1=0,m2=0,l1=0,l2=0,max=0,min=double.PositiveInfinity,maxCoeff=double.NegativeInfinity;
        for(int j=0;j<N;j++){double w=reconstruction[index,j];sum+=w;m1+=w*j;m2+=w*j*j;l1+=Math.Abs(w);l2+=w*w;max=Math.Max(max,Math.Abs(w));min=Math.Min(min,w);maxCoeff=Math.Max(maxCoeff,w);}
        return new(index,sum,m1-index,m2-index*(double)index,l1,Math.Sqrt(l2),max,min,maxCoeff);
    }
    public double[] Smooth(ReadOnlySpan<double> x)
    {
        var y=new double[N];
        for(int i=0;i<N;i++)
        {
            double sum=0;
            for(int j=0;j<N;j++)sum+=reconstruction[i,j]*x[j];
            y[i]=sum;
        }
        return y;
    }

    public double[] MinimizeBlock(ReadOnlySpan<double> x,int start,int count)
    {
        var a=new double[count,count];var b=new double[count];
        for(int i=0;i<count;i++){int gi=start+i;for(int j=0;j<count;j++)a[i,j]=L[gi,start+j];double rhs=0;for(int j=0;j<N;j++)if(j<start||j>=start+count)rhs+=L[gi,j]*x[j];b[i]=-rhs;}
        return Solve(a,b);
    }
    static double[] Solve(double[,] a,double[] b)
    { int n=b.Length;var m=(double[,])a.Clone();var x=(double[])b.Clone();double ridge=1e-12*Math.Max(1.0,Enumerable.Range(0,n).Select(i=>Math.Abs(m[i,i])).DefaultIfEmpty(1).Max());for(int i=0;i<n;i++)m[i,i]+=ridge;for(int k=0;k<n;k++){int p=k;for(int i=k+1;i<n;i++)if(Math.Abs(m[i,k])>Math.Abs(m[p,k]))p=i;if(p!=k){for(int j=k;j<n;j++)(m[k,j],m[p,j])=(m[p,j],m[k,j]);(x[k],x[p])=(x[p],x[k]);}double d=m[k,k];if(Math.Abs(d)<1e-20)continue;for(int i=k+1;i<n;i++){double f=m[i,k]/d;for(int j=k;j<n;j++)m[i,j]-=f*m[k,j];x[i]-=f*x[k];}}var z=new double[n];for(int i=n-1;i>=0;i--){double s=x[i];for(int j=i+1;j<n;j++)s-=m[i,j]*z[j];z[i]=Math.Abs(m[i,i])<1e-20?0:s/m[i,i];}return z; }

    public NumericalDiagnostics GetNumericalDiagnostics()
    {
        var entries=new List<double>(N*N);
        double maxL=double.NegativeInfinity;
        for(int i=0;i<N;i++) for(int j=0;j<N;j++)
        {
            double v=L[i,j];
            maxL=Math.Max(maxL,v);
            if(v!=0.0) entries.Add(Math.Abs(v));
        }
        entries.Sort();
        double medianNonZeroAbsL=entries.Count==0?0.0:
            (entries.Count%2!=0?entries[entries.Count/2]:0.5*(entries[entries.Count/2-1]+entries[entries.Count/2]));
        return new(maxL,medianNonZeroAbsL);
    }

}

public sealed record ReconstructionKernelDiagnostics(int Index,double DcGain,double LinearMomentError,double QuadraticMomentError,double L1Norm,double L2Norm,double MaxAbsCoefficient,double MinCoefficient,double MaxCoefficient);

public sealed record NumericalDiagnostics(double MaxOperatorEntry,double MedianNonZeroAbsOperatorEntry);
