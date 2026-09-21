using System.Diagnostics;
using System.Numerics;

namespace TrajectoryEigenFilterCore;

public sealed class TrajectoryFilter
{
    public FilterSettings Settings { get; }
    public ComplexityModel Model { get; private set; }
    readonly Queue<Vector2> input=new();
    readonly Dictionary<long,AdaptiveEstimateAccumulator> adaptiveAccumulators=new();
    long pushedSamples; int adaptiveLag=-1; AdaptiveKernel[] adaptiveKernels=Array.Empty<AdaptiveKernel>();
    public double SampleRateHz=>Model.SampleRateHz; public int WindowSamples=>Model.N;
    public ProcessingTimings Timings { get; }=new();

    public TrajectoryFilter(FilterSettings settings,double? sampleRateHz=null){Settings=settings;Model=new(settings,sampleRateHz??settings.NominalSampleRateHz);}
    public void ReconfigureRate(double hz){if(hz<=0)throw new ArgumentOutOfRangeException(nameof(hz));var old=input.ToArray();Model=new(Settings,hz);input.Clear();foreach(var p in old.TakeLast(Model.N))input.Enqueue(p);ResetAdaptiveState();}

    public Vector2[] SnapshotInput() => input.ToArray();
    public void ResetInput(IEnumerable<Vector2> points)
    {
        input.Clear();
        foreach (var p in points.TakeLast(Model.N)) input.Enqueue(p);
        ResetAdaptiveState();
    }

    public bool Push(Vector2 point,out FilteredWindow result)=>Push(point,null,out result);
    public bool Push(Vector2 point,int? targetIndex,out FilteredWindow result)
    {
        input.Enqueue(point);while(input.Count>Model.N)input.Dequeue();if(input.Count<Model.N){result=default!;return false;}
        long t0=Stopwatch.GetTimestamp();
        var src=input.ToArray();var x=new double[Model.N];var y=new double[Model.N];for(int i=0;i<Model.N;i++){x[i]=src[i].X;y[i]=src[i].Y;}
        Timings.WindowPreparationTicks+=Stopwatch.GetTimestamp()-t0;

        var replaced=RejectOutliers(x,y);

        t0=Stopwatch.GetTimestamp();
        double[] sx,sy;
        if(targetIndex is int ti){sx=new[]{Model.SmoothAt(x,ti)};sy=new[]{Model.SmoothAt(y,ti)};}
        else {sx=Model.Smooth(x);sy=Model.Smooth(y);}
        Timings.ReconstructionTicks+=Stopwatch.GetTimestamp()-t0;

        t0=Stopwatch.GetTimestamp();double complexity=Model.Complexity(x)+Model.Complexity(y);Timings.ComplexityTicks+=Stopwatch.GetTimestamp()-t0;
        result=new FilteredWindow(sx,sy,replaced,complexity,targetIndex);Timings.Windows++;return true;
    }

    void ResetAdaptiveState(){adaptiveAccumulators.Clear();pushedSamples=0;adaptiveLag=-1;adaptiveKernels=Array.Empty<AdaptiveKernel>();}

    // Reconstruct each physical sample repeatedly as it moves inward from the
    // newest boundary.  Early/low-lookahead estimates use a higher hard cutoff
    // to avoid the observed low-rank boundary delay.  Estimates are weighted by
    // inverse unit-white-noise transmission.  R=0 is intentionally just the
    // zero-lookahead high-rank estimator; there are no negative-delay terms.
    public bool PushAdaptiveDelayAverage(Vector2 point,int lag,out FilteredWindow result)
    {
        if(lag<0||lag>=Model.N)throw new ArgumentOutOfRangeException(nameof(lag));
        if(adaptiveLag!=lag)
        {
            adaptiveAccumulators.Clear();pushedSamples=0;adaptiveLag=lag;
            adaptiveKernels=new AdaptiveKernel[lag+1];
            for(int r=0;r<=lag;r++)
            {
                double cutoff=Model.AdaptiveLambdaForLookaheadSamples(r);
                var built=Model.BuildHardCutoffKernel(Model.N-1-r,cutoff);
                adaptiveKernels[r]=new AdaptiveKernel(built.Kernel,1.0/Math.Max(built.NoiseTransmission,1e-12));
            }
        }
        long sequence=pushedSamples++;
        input.Enqueue(point);while(input.Count>Model.N)input.Dequeue();
        if(input.Count<Model.N){result=default!;return false;}
        long t0=Stopwatch.GetTimestamp();
        var src=input.ToArray();var x=new double[Model.N];var y=new double[Model.N];
        for(int i=0;i<Model.N;i++){x[i]=src[i].X;y[i]=src[i].Y;}
        Timings.WindowPreparationTicks+=Stopwatch.GetTimestamp()-t0;
        var replaced=RejectOutliers(x,y);
        t0=Stopwatch.GetTimestamp();
        for(int r=0;r<=lag;r++)
        {
            var ak=adaptiveKernels[r]; double ex=0,ey=0;
            for(int j=0;j<Model.N;j++){double h=ak.Kernel[j];ex+=h*x[j];ey+=h*y[j];}
            double weight=ak.Weight;
            long target=sequence-r;
            if(!adaptiveAccumulators.TryGetValue(target,out var a))adaptiveAccumulators[target]=a=new AdaptiveEstimateAccumulator();
            a.X+=weight*ex;a.Y+=weight*ey;a.Weight+=weight;a.Count++;
        }
        Timings.ReconstructionTicks+=Stopwatch.GetTimestamp()-t0;
        long outputSequence=sequence-lag;
        if(!adaptiveAccumulators.TryGetValue(outputSequence,out var output)||output.Count<lag+1){result=default!;return false;}
        adaptiveAccumulators.Remove(outputSequence);
        foreach(var key in adaptiveAccumulators.Keys.Where(k=>k<outputSequence).ToArray())adaptiveAccumulators.Remove(key);
        t0=Stopwatch.GetTimestamp();double complexity=Model.Complexity(x)+Model.Complexity(y);Timings.ComplexityTicks+=Stopwatch.GetTimestamp()-t0;
        int targetIndex=Model.N-1-lag;
        result=new FilteredWindow(new[]{output.X/output.Weight},new[]{output.Y/output.Weight},replaced,complexity,targetIndex);Timings.Windows++;return true;
    }
    sealed class AdaptiveEstimateAccumulator{public double X,Y,Weight;public int Count;}
    readonly record struct AdaptiveKernel(double[] Kernel,double Weight);

    List<OutlierReplacement> RejectOutliers(double[] x,double[] y)
    {
        var accepted=new List<OutlierReplacement>();int n=Model.N;
        for(int pass=0;pass<Settings.MaxOutlierPasses;pass++)
        {
            // Localize anomalies using the spectral reconstruction residual itself.
            // Scores are computed before any replacement in this pass, so an
            // anomalous sample remains localized at its own raw coordinate rather
            // than through the oscillatory stencil of L*x/L_ii.
            long t0=Stopwatch.GetTimestamp();
            var smoothX=Model.Smooth(x);
            var smoothY=Model.Smooth(y);
            var residual=new double[n];
            for(int i=0;i<n;i++)
            {
                double rx=x[i]-smoothX[i], ry=y[i]-smoothY[i];
                residual[i]=Math.Sqrt(rx*rx+ry*ry);
            }
            Timings.ResidualTicks+=Stopwatch.GetTimestamp()-t0;

            t0=Stopwatch.GetTimestamp();
            var interior=residual.Skip(1).Take(n-2).ToArray();double med=Median(interior);var dev=interior.Select(v=>Math.Abs(v-med)).ToArray();double sigma=1.4826*Median(dev);
            // Keep a small scale floor so a nearly noiseless synthetic window remains well-defined.
            sigma=Math.Max(sigma,Math.Max(1e-9,med*1e-6));
            // The newest sample must be eligible for rejection. At zero output latency
            // it is exactly the sample exposed by the no-extrapolation boundary path.
            // Keep the robust baseline interior-only, but score every sample against it.
            int best=0;double bestZ=double.NegativeInfinity;for(int i=0;i<n;i++){double z=(residual[i]-med)/sigma;if(z>bestZ){bestZ=z;best=i;}}
            Timings.RobustScoringTicks+=Stopwatch.GetTimestamp()-t0;
            if(bestZ<=Settings.OutlierThreshold)break;

            // Only after an anomaly exists, expand to adjacent samples that are also
            // substantially abnormal. This keeps 2/3-point block solving off the normal path.
            int start=best,end=best;double adjacentThreshold=Math.Max(2.5,Settings.OutlierThreshold*0.5);
            while(start>0 && end-start+1<Settings.MaxOutlierBlock && (residual[start-1]-med)/sigma>=adjacentThreshold)start--;
            while(end<n-1 && end-start+1<Settings.MaxOutlierBlock && (residual[end+1]-med)/sigma>=adjacentThreshold)end++;
            int count=end-start+1;

            // Detection and replacement are deliberately separate.  The old
            // minimum-complexity block solve can be extremely ill-conditioned and
            // generate coordinates far outside the observed trajectory.  Replace
            // anomalous samples with a local quadratic trajectory estimate instead;
            // this is consistent with the filter's analytic nullspace and remains
            // bounded by a robust local sanity guard.
            t0=Stopwatch.GetTimestamp();
            var rxv=LocalQuadraticReplacement(x,start,count);
            var ryv=LocalQuadraticReplacement(y,start,count);
            var before=new Vector2[count];var after=new Vector2[count];
            for(int i=0;i<count;i++)
            {
                before[i]=new((float)x[start+i],(float)y[start+i]);
                x[start+i]=rxv[i];y[start+i]=ryv[i];
                after[i]=new((float)rxv[i],(float)ryv[i]);
            }
            accepted.Add(new(start,before,after,bestZ));Timings.ReplacementTicks+=Stopwatch.GetTimestamp()-t0;
        }
        return accepted;
    }

    static double[] LocalQuadraticReplacement(double[] x,int start,int count)
    {
        int n=x.Length, radius=8, lo=Math.Max(0,start-radius), hi=Math.Min(n-1,start+count-1+radius);
        var idx=new List<int>();
        for(int j=lo;j<=hi;j++) if(j<start||j>=start+count) idx.Add(j);
        // Need at least three support points.  With normal interior operation this
        // branch is never used, but keep a stable linear fallback for boundaries.
        if(idx.Count<3)
        {
            var z=new double[count];
            int l=Math.Max(0,start-1), r=Math.Min(n-1,start+count);
            for(int k=0;k<count;k++) z[k]=l==r?x[l]:x[l]+(x[r]-x[l])*(start+k-l)/(double)(r-l);
            return z;
        }

        // Fit a + b*t + c*t^2 in coordinates centered on the replacement block.
        // Centering keeps the tiny 3x3 least-squares system well conditioned.
        double center=start+(count-1)*0.5;
        double s0=0,s1=0,s2=0,s3=0,s4=0,y0=0,y1=0,y2=0;
        foreach(int j in idx)
        {
            double t=j-center, v=x[j], t2=t*t;
            s0+=1;s1+=t;s2+=t2;s3+=t2*t;s4+=t2*t2;
            y0+=v;y1+=t*v;y2+=t2*v;
        }
        var a=new double[,]{{s0,s1,s2},{s1,s2,s3},{s2,s3,s4}};
        var b=new[]{y0,y1,y2};
        var coef=Solve3(a,b);
        var raw=new double[count];
        for(int k=0;k<count;k++){double t=start+k-center;raw[k]=coef[0]+coef[1]*t+coef[2]*t*t;}

        // Robust sanity envelope from the support samples.  This is not the
        // detector threshold: it only prevents a replacement estimator from
        // inventing a remote coordinate if its fit ever becomes pathological.
        var vals=idx.Select(j=>x[j]).ToArray();double med=Median((double[])vals.Clone());
        var dev=vals.Select(v=>Math.Abs(v-med)).ToArray();double mad=1.4826*Median(dev);
        double min=vals.Min(),max=vals.Max(), span=Math.Max(max-min,Math.Max(1.0,6*mad));
        double guardLo=min-span,guardHi=max+span;
        for(int k=0;k<count;k++) raw[k]=Math.Clamp(raw[k],guardLo,guardHi);
        return raw;
    }

    static double[] Solve3(double[,] a,double[] b)
    {
        var m=(double[,])a.Clone();var v=(double[])b.Clone();
        for(int k=0;k<3;k++)
        {
            int p=k;for(int i=k+1;i<3;i++)if(Math.Abs(m[i,k])>Math.Abs(m[p,k]))p=i;
            if(p!=k){for(int j=k;j<3;j++)(m[k,j],m[p,j])=(m[p,j],m[k,j]);(v[k],v[p])=(v[p],v[k]);}
            double d=m[k,k];if(Math.Abs(d)<1e-14)return new[]{b[0]/Math.Max(1.0,a[0,0]),0.0,0.0};
            for(int i=k+1;i<3;i++){double f=m[i,k]/d;for(int j=k;j<3;j++)m[i,j]-=f*m[k,j];v[i]-=f*v[k];}
        }
        var z=new double[3];for(int i=2;i>=0;i--){double q=v[i];for(int j=i+1;j<3;j++)q-=m[i,j]*z[j];z[i]=q/m[i,i];}return z;
    }

    static double Median(double[] a){if(a.Length==0)return 0;Array.Sort(a);int m=a.Length/2;return (a.Length&1)==1?a[m]:.5*(a[m-1]+a[m]);}
}

public sealed class ProcessingTimings
{
    public long Windows,WindowPreparationTicks,ResidualTicks,RobustScoringTicks,ReplacementTicks,ReconstructionTicks,ComplexityTicks;
    public double Ms(long ticks)=>ticks*1000.0/Stopwatch.Frequency;
}
public sealed record OutlierReplacement(int Start,Vector2[] Before,Vector2[] After,double RelativeResidualScore);
public sealed record FilteredWindow(double[] X,double[] Y,IReadOnlyList<OutlierReplacement> Replacements,double Complexity,int? TargetIndex=null)
{
    public Vector2 At(int index){if(TargetIndex is int ti){if(index!=ti)throw new ArgumentOutOfRangeException(nameof(index));return new((float)X[0],(float)Y[0]);}return new((float)X[index],(float)Y[index]);}
    public Vector2 At(double index)
    {
        if(TargetIndex is not null)return At(TargetIndex.Value);
        int n=X.Length;if(n==0)return default;if(n==1)return At(0);
        if(index<=0)return At(0);
        if(index>=n-1)
        {
            double t=index-(n-1),dx=X[n-1]-X[n-2],dy=Y[n-1]-Y[n-2];
            // Respect the DCT-II reflection boundary at N-1/2. A quadratic
            // continuation reaches zero slope there; only later scheduler
            // requests use the straight-line fallback for delayed input.
            double bx=X[n-1]+dx/8.0,by=Y[n-1]+dy/8.0;
            if(t<=0.5)
            {
                double q=t-t*t;
                return new((float)(X[n-1]+0.5*dx*q),(float)(Y[n-1]+0.5*dy*q));
            }
            return new((float)(bx+dx*(t-0.5)),(float)(by+dy*(t-0.5)));
        }
        int i=(int)Math.Floor(index);double interpolationFraction=index-i;
        return new((float)(X[i]+(X[i+1]-X[i])*interpolationFraction),(float)(Y[i]+(Y[i+1]-Y[i])*interpolationFraction));
    }
}
