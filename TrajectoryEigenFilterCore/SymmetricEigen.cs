namespace TrajectoryEigenFilterCore;

// Dependency-free Jacobi eigensolver for small symmetric matrices. Construction
// happens only when the sample rate/window changes, never per tablet report.
internal static class SymmetricEigen
{
    public static (double[] values, double[,] vectors) Decompose(double[,] input)
    {
        int n = input.GetLength(0);
        var a = (double[,])input.Clone();
        var v = new double[n,n];
        for (int i=0;i<n;i++) v[i,i]=1;
        int maxIter = Math.Max(64, 30*n*n);
        for (int iter=0; iter<maxIter; iter++)
        {
            int p=0,q=1; double max=0;
            for(int i=0;i<n;i++) for(int j=i+1;j<n;j++)
            { var z=Math.Abs(a[i,j]); if(z>max){max=z;p=i;q=j;} }
            // The spectrum can span many orders of magnitude.  The previous 1e-13
            // relative cutoff stopped while off-diagonal terms were still comparable
            // to the smallest useful eigenvalues.  Continue close to double precision.
            if(max < 2e-15 * Math.Max(1.0, MaxDiagonal(a))) break;
            double app=a[p,p], aqq=a[q,q], apq=a[p,q];
            double phi=0.5*Math.Atan2(2*apq, aqq-app);
            double c=Math.Cos(phi), s=Math.Sin(phi);
            for(int k=0;k<n;k++) if(k!=p && k!=q)
            {
                double akp=a[k,p], akq=a[k,q];
                a[k,p]=a[p,k]=c*akp-s*akq;
                a[k,q]=a[q,k]=s*akp+c*akq;
            }
            a[p,p]=c*c*app-2*s*c*apq+s*s*aqq;
            a[q,q]=s*s*app+2*s*c*apq+c*c*aqq;
            a[p,q]=a[q,p]=0;
            for(int k=0;k<n;k++)
            {
                double vkp=v[k,p], vkq=v[k,q];
                v[k,p]=c*vkp-s*vkq; v[k,q]=s*vkp+c*vkq;
            }
        }
        var vals=new double[n]; for(int i=0;i<n;i++) vals[i]=a[i,i];
        var idx=Enumerable.Range(0,n).OrderBy(i=>vals[i]).ToArray();
        var sv=new double[n]; var qv=new double[n,n];
        for(int k=0;k<n;k++){sv[k]=vals[idx[k]]; for(int i=0;i<n;i++) qv[i,k]=v[i,idx[k]];}
        return (sv,qv);
    }
    static double MaxDiagonal(double[,] a){double m=0;for(int i=0;i<a.GetLength(0);i++)m=Math.Max(m,Math.Abs(a[i,i]));return m;}
}
