using System.Globalization;
using System.Numerics;
using System.Text;
using System.Text.Json;
using System.Text.RegularExpressions;
using System.Diagnostics;
using TrajectoryEigenFilterCore;

CultureInfo.CurrentCulture = CultureInfo.InvariantCulture;

if (args.Length == 0 || args.Contains("--help"))
{
    Help();
    return;
}

try
{
    var opt = Options.Parse(args);
    var samples = RecordingParser.Read(opt.Input).ToList();
    if (samples.Count == 0) throw new InvalidOperationException("No IntuosTabletReport Position:[x,y] samples were found.");

    var settings = new FilterSettings
    {
        NominalSampleRateHz = opt.Rate,
        WindowMilliseconds = opt.WindowMs,
        MinimumDifferenceOrder = opt.MinOrder,
        MaximumDifferenceOrder = opt.MaxOrder,
        LambdaCutoff = opt.LambdaCutoff,
        UseAdaptiveLambda = opt.UseAdaptiveLambda,
        LocalDifferenceStrength = opt.LocalStrength,
        AdaptiveLambdaAt0Ms = opt.AdaptiveLambda0,
        AdaptiveLambdaAt2Ms = opt.AdaptiveLambda2,
        AdaptiveLambdaAt5Ms = opt.AdaptiveLambda5,
        AdaptiveLambdaAt10Ms = opt.AdaptiveLambda10,
        AdaptiveLambdaAt20Ms = opt.AdaptiveLambda20,
        OutlierThreshold = opt.OutlierThreshold,
        MaxOutlierBlock = opt.MaxOutlierBlock,
        MaxOutlierPasses = opt.MaxOutlierPasses
    };

    var modelTimer=Stopwatch.StartNew();
    var filter = new TrajectoryFilter(settings, opt.Rate);
    modelTimer.Stop();
    var diag=filter.Model.GetNumericalDiagnostics();
    int requestedLag = (int)Math.Round(opt.LookaheadMs * opt.Rate / 1000.0);
    int lag = Math.Min(requestedLag, Math.Max(0, filter.WindowSamples - 2));
    int evalIndex = filter.WindowSamples - 1 - lag;
    var rows = new List<ResultRow>();

    var processingTimer=Stopwatch.StartNew();
    for (int i = 0; i < samples.Count; i++)
    {
        Vector2 inputPoint = new((float)samples[i].X, (float)samples[i].Y);
        FilteredWindow w;
        bool ready;
        Vector2 filtered;
        if (opt.UseAdaptiveLambda)
        {
            ready = filter.PushAdaptiveDelayAverage(inputPoint, lag, out w);
            if (!ready) continue;
            filtered = w.At(evalIndex);
        }
        else
        {
            ready = filter.Push(inputPoint, out w);
            if (!ready) continue;
            filtered = w.At(evalIndex);
        }
        int global = i - lag;
        if (global < 0 || global >= samples.Count) continue;
        var raw = samples[global];

        bool replaced = false;
        bool windowHadReplacement = w.Replacements.Count > 0;
        double outlierScore = 0;
        int replacementStart = -1;
        double replacementOriginalX = double.NaN, replacementOriginalY = double.NaN;
        double replacementNewX = double.NaN, replacementNewY = double.NaN;
        foreach (var r in w.Replacements)
        {
            outlierScore = Math.Max(outlierScore, r.RelativeResidualScore);
            if (replacementStart < 0)
            {
                replacementStart = r.Start;
                replacementOriginalX = r.Before[0].X; replacementOriginalY = r.Before[0].Y;
                replacementNewX = r.After[0].X; replacementNewY = r.After[0].Y;
            }
            if (evalIndex >= r.Start && evalIndex < r.Start + r.Before.Length) replaced = true;
        }

        rows.Add(new ResultRow(
            global,
            global * 1000.0 / opt.Rate,
            raw.X, raw.Y,
            filtered.X, filtered.Y,
            filtered.X - raw.X, filtered.Y - raw.Y,
            replaced, windowHadReplacement, outlierScore, replacementStart,
            replacementOriginalX, replacementOriginalY, replacementNewX, replacementNewY, w.Complexity));
    }

    processingTimer.Stop();
    string csv = Path.GetFullPath(opt.OutputCsv);
    string html = Path.GetFullPath(opt.OutputHtml ?? Path.ChangeExtension(csv, ".html"));
    Directory.CreateDirectory(Path.GetDirectoryName(csv)!);
    CsvWriter.Write(csv, rows);
    HtmlWriter.Write(html, rows, opt, filter.WindowSamples, lag);

    Console.WriteLine($"Parsed reports : {samples.Count}");
    Console.WriteLine($"Model build     : {modelTimer.Elapsed.TotalMilliseconds:F1} ms");
    Console.WriteLine($"Filter processing: {processingTimer.Elapsed.TotalMilliseconds:F1} ms ({(samples.Count>0?processingTimer.Elapsed.TotalMilliseconds/samples.Count:0):F3} ms/input)");
    var pt=filter.Timings;
    Console.WriteLine("Processing breakdown:");
    Console.WriteLine($"  window preparation : {pt.Ms(pt.WindowPreparationTicks):F1} ms");
    Console.WriteLine($"  spectral residuals : {pt.Ms(pt.ResidualTicks):F1} ms");
    Console.WriteLine($"  robust scoring     : {pt.Ms(pt.RobustScoringTicks):F1} ms");
    Console.WriteLine($"  block replacement  : {pt.Ms(pt.ReplacementTicks):F1} ms");
    Console.WriteLine($"  reconstruction     : {pt.Ms(pt.ReconstructionTicks):F1} ms");
    Console.WriteLine($"  complexity output  : {pt.Ms(pt.ComplexityTicks):F1} ms");
    Console.WriteLine("Numerical diagnostics:");
    Console.WriteLine($"  Max L_ij                  : {diag.MaxOperatorEntry:E6}");
    Console.WriteLine($"  median non-zero |L_ij|    : {diag.MedianNonZeroAbsOperatorEntry:E6}");
    Console.WriteLine("  eigenvalues and filtering ratios:");
    for (int k = 0; k < filter.Model.Eigenvalues.Length; k++)
        Console.WriteLine($"    lambda[{k,2}] = {filter.Model.Eigenvalues[k]:E17}    g = {filter.Model.Gain(k):E17}");
    if (opt.PrintEigen >= 0)
        PrintEigenMatrix(filter.Model, opt.PrintEigen);
    else if (opt.PrintEigen == -1)
        PrintEigenMatrix(filter.Model, -1);
    Console.WriteLine($"Window         : {filter.WindowSamples} samples ({opt.WindowMs:F2} ms)");
    Console.WriteLine($"Analysis lag   : {lag} samples ({lag * 1000.0 / opt.Rate:F3} ms)");
    Console.WriteLine($"Output rows    : {rows.Count}");
    Console.WriteLine($"Replaced target samples: {rows.Count(r => r.Replaced)}");
    if (rows.Count > 0)
    {
        var absX=rows.Select(r=>Math.Abs(r.ErrorX)).OrderBy(v=>v).ToArray();
        var absY=rows.Select(r=>Math.Abs(r.ErrorY)).OrderBy(v=>v).ToArray();
        Console.WriteLine("Absolute filtering |filtered-raw|:");
        Console.WriteLine($"  X median / top 25% / top 5%: {Percentile(absX,0.50):F4} / {Percentile(absX,0.75):F4} / {Percentile(absX,0.95):F4}");
        Console.WriteLine($"  Y median / top 25% / top 5%: {Percentile(absY,0.50):F4} / {Percentile(absY,0.75):F4} / {Percentile(absY,0.95):F4}");
    }
    Console.WriteLine($"CSV            : {csv}");
    Console.WriteLine($"Interactive HTML: {html}");
}
catch (Exception ex)
{
    Console.Error.WriteLine("Error: " + ex.Message);
    Environment.ExitCode = 1;
}

static void PrintEigenMatrix(ComplexityModel model, int selection)
{
    if (selection >= model.N)
        throw new ArgumentOutOfRangeException(nameof(selection), $"Eigenvector index {selection} is outside 0..{model.N-1}.");
    int first=selection<0 ? 0 : selection;
    int last=selection<0 ? model.N-1 : selection;
    Console.WriteLine("Eigenvectors (rows are eigenvector indices; columns are sample indices):");
    for(int k=first;k<=last;k++)
    {
        var values=new string[model.N];
        for(int i=0;i<model.N;i++) values[i]=model.Q[i,k].ToString("G3",CultureInfo.InvariantCulture);
        Console.WriteLine($"  q[{k,2}] = [ {string.Join(" ",values)} ]");
    }
}

static double Percentile(double[] sorted,double p)
{
    if(sorted.Length==0)return double.NaN;
    double pos=(sorted.Length-1)*Math.Clamp(p,0.0,1.0);
    int lo=(int)Math.Floor(pos),hi=(int)Math.Ceiling(pos);
    if(lo==hi)return sorted[lo];
    double f=pos-lo;
    return sorted[lo]*(1.0-f)+sorted[hi]*f;
}

static void Help()
{
    Console.WriteLine("""
TrajectoryEigenFilter.Analysis

Usage:
  dotnet run --project TrajectoryEigenFilter.Analysis -c Release -- <recording.txt> [options]

Options:
  --rate <Hz>                 Nominal sample rate (default 700)
  --window-ms <ms>            Physical filter window (default 100)
  --lookahead-ms <ms>         Offline evaluation look-ahead/buffer (default 5; range 0-20)
  --min-order <n>             Minimum finite-difference order (default 0)
  --max-order <n>             Maximum finite-difference order (default 8)
  --lambda-cutoff <x>         Hard eigenvalue cutoff when adaptive lambda is off (default 1.5)
  --print-eigen <n>           Print eigenvector matrix at 3 significant digits; -1 = all
  --local-strength <x>        Bulk finite-difference family weight (default 1)
  --adaptive-lambda <bool>     Enable latency-dependent cutoff/estimate averaging (default true)
  --adaptive-lambda-0 <x>     Hard cutoff for zero-lookahead estimates (default 1.5)
  --adaptive-lambda-2 <x>     Hard cutoff at 2 ms lookahead (default 1.2)
  --adaptive-lambda-5 <x>     Hard cutoff at 5 ms lookahead (default 1.05)
  --adaptive-lambda-10 <x>    Hard cutoff at 10 ms lookahead (default 1.04)
  --adaptive-lambda-20 <x>    Hard cutoff at >=20 ms lookahead (default 1.03)
  --outlier-threshold <x>     Relative residual MAD-z threshold (default 6)
  --max-outlier-block <n>     Maximum contiguous replacement (default 3)
  --max-outlier-passes <n>    Replacement passes/window (default 2)
  --output <file.csv>         CSV output (default analysis.csv)
  --html <file.html>          Interactive plot (default: CSV name with .html)

The recording Delta fields are parsed only implicitly as text and are NOT used as
trajectory time. Sample index / nominal rate defines physical time.
""");
}

sealed record Options(string Input, double Rate, double WindowMs, double LookaheadMs,
    int MinOrder, int MaxOrder, double LambdaCutoff, double LocalStrength, bool UseAdaptiveLambda, double AdaptiveLambda0, double AdaptiveLambda2, double AdaptiveLambda5, double AdaptiveLambda10, double AdaptiveLambda20, double OutlierThreshold,
    int MaxOutlierBlock, int MaxOutlierPasses, int PrintEigen, string OutputCsv, string? OutputHtml)
{
    public static Options Parse(string[] a)
    {
        string input = a[0];
        double rate=700, window=100, lookahead=5, lambdaCutoff=1.5, localStrength=1, adaptiveLambda0=1.5, adaptiveLambda2=1.2, adaptiveLambda5=1.05, adaptiveLambda10=1.04, adaptiveLambda20=1.03, threshold=6;
        bool useAdaptiveLambda=true;
        int min=0, max=8, block=3, passes=2, printEigen=-2;
        string output="analysis.csv"; string? html=null;
        for(int i=1;i<a.Length;i++)
        {
            string Next() => ++i<a.Length ? a[i] : throw new ArgumentException($"Missing value after {a[i-1]}");
            switch(a[i])
            {
                case "--rate": rate=double.Parse(Next(),CultureInfo.InvariantCulture); break;
                case "--window-ms": window=double.Parse(Next(),CultureInfo.InvariantCulture); break;
                case "--lookahead-ms": lookahead=double.Parse(Next(),CultureInfo.InvariantCulture); break;
                case "--min-order": min=int.Parse(Next()); break;
                case "--max-order": max=int.Parse(Next()); break;
                case "--lambda-cutoff": lambdaCutoff=double.Parse(Next(),CultureInfo.InvariantCulture); break;
                case "--adaptive-lambda": useAdaptiveLambda=bool.Parse(Next()); break;
                case "--print-eigen": printEigen=int.Parse(Next(),CultureInfo.InvariantCulture); break;
                case "--local-strength": localStrength=double.Parse(Next(),CultureInfo.InvariantCulture); break;
                case "--adaptive-lambda-0": adaptiveLambda0=double.Parse(Next(),CultureInfo.InvariantCulture); break;
                case "--adaptive-lambda-2": adaptiveLambda2=double.Parse(Next(),CultureInfo.InvariantCulture); break;
                case "--adaptive-lambda-5": adaptiveLambda5=double.Parse(Next(),CultureInfo.InvariantCulture); break;
                case "--adaptive-lambda-10": adaptiveLambda10=double.Parse(Next(),CultureInfo.InvariantCulture); break;
                case "--adaptive-lambda-20": adaptiveLambda20=double.Parse(Next(),CultureInfo.InvariantCulture); break;
                case "--outlier-threshold": threshold=double.Parse(Next(),CultureInfo.InvariantCulture); break;
                case "--max-outlier-block": block=int.Parse(Next()); break;
                case "--max-outlier-passes": passes=int.Parse(Next()); break;
                case "--output": output=Next(); break;
                case "--html": html=Next(); break;
                default: throw new ArgumentException($"Unknown option: {a[i]}");
            }
        }
        if(rate<=0 || window<=0 || !double.IsFinite(lookahead) || lookahead<0 || lookahead>20 || !double.IsFinite(lambdaCutoff) || lambdaCutoff<0 || localStrength<=0 || adaptiveLambda0<=0 || adaptiveLambda2<=0 || adaptiveLambda5<=0 || adaptiveLambda10<=0 || adaptiveLambda20<=0 || printEigen < -2)
            throw new ArgumentException("Rate/window and local strength must be positive, lookahead must be 0-20 ms, lambda cutoffs must be nonnegative/positive, and --print-eigen must be -1 or a nonnegative index.");
        return new(input,rate,window,lookahead,min,max,lambdaCutoff,localStrength,useAdaptiveLambda,adaptiveLambda0,adaptiveLambda2,adaptiveLambda5,adaptiveLambda10,adaptiveLambda20,threshold,block,passes,printEigen,output,html);
    }
}

readonly record struct TabletSample(double X,double Y);
static class RecordingParser
{
    static readonly Regex Position = new(@"Position:\[\s*(-?\d+(?:\.\d+)?),\s*(-?\d+(?:\.\d+)?)\].*ReportType:IntuosTabletReport", RegexOptions.Compiled);
    public static IEnumerable<TabletSample> Read(string path)
    {
        foreach(var line in File.ReadLines(path))
        {
            var m=Position.Match(line); if(!m.Success) continue;
            yield return new(double.Parse(m.Groups[1].Value,CultureInfo.InvariantCulture), double.Parse(m.Groups[2].Value,CultureInfo.InvariantCulture));
        }
    }
}

readonly record struct ResultRow(int Sample,double TimeMs,double RawX,double RawY,double FilteredX,double FilteredY,
    double ErrorX,double ErrorY,bool Replaced,bool WindowHadReplacement,double OutlierScore,int ReplacementStart,
    double ReplacementOriginalX,double ReplacementOriginalY,double ReplacementNewX,double ReplacementNewY,double Complexity);

static class CsvWriter
{
    public static void Write(string path,List<ResultRow> rows)
    {
        using var w=new StreamWriter(path,false,new UTF8Encoding(false));
        w.WriteLine("sample,time_ms,raw_x,raw_y,filtered_x,filtered_y,filtered_minus_raw_x,filtered_minus_raw_y,replaced,window_had_replacement,outlier_score,replacement_start,replacement_original_x,replacement_original_y,replacement_new_x,replacement_new_y,complexity");
        foreach(var r in rows)
            w.WriteLine(FormattableString.Invariant($"{r.Sample},{r.TimeMs:R},{r.RawX:R},{r.RawY:R},{r.FilteredX:R},{r.FilteredY:R},{r.ErrorX:R},{r.ErrorY:R},{(r.Replaced?1:0)},{(r.WindowHadReplacement?1:0)},{r.OutlierScore:R},{r.ReplacementStart},{r.ReplacementOriginalX:R},{r.ReplacementOriginalY:R},{r.ReplacementNewX:R},{r.ReplacementNewY:R},{r.Complexity:R}"));
    }
}

static class HtmlWriter
{
    public static void Write(string path,List<ResultRow> rows,Options o,int windowSamples,int lag)
    {
        var payload = rows.Select(r => new object[]{r.TimeMs,r.RawX,r.RawY,r.FilteredX,r.FilteredY,r.ErrorX,r.ErrorY,r.Replaced?1:0,r.OutlierScore}).ToArray();
        string data=JsonSerializer.Serialize(payload);
        string meta=JsonSerializer.Serialize(new {rate=o.Rate,windowMs=o.WindowMs,lookaheadMs=lag*1000.0/o.Rate,windowSamples});
        string html = """
<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Trajectory Eigen Filter Analysis</title>
<style>
:root{font-family:system-ui,sans-serif;color-scheme:light dark}body{margin:18px}#bar{display:flex;gap:12px;flex-wrap:wrap;align-items:center;margin-bottom:10px}button{padding:6px 10px}canvas{width:100%;height:340px;border:1px solid #888;touch-action:none;display:block}.small{opacity:.75;font-size:13px}.legend{display:flex;gap:18px;margin:6px 0 14px}.sw{display:inline-block;width:18px;border-top:3px solid;margin-right:5px}.raw{border-color:#888}.flt{border-color:#2878d0}.err{border-color:#d06b28}</style></head><body>
<h2>Trajectory Eigen Filter Analysis</h2><div id="bar"><button id="reset">Reset zoom</button><span id="range"></span><span class="small">Wheel: zoom · drag: pan · double-click: reset</span></div>
<div id="meta" class="small"></div>
<h3>X position</h3><div class="legend"><span><i class="sw raw"></i>Raw</span><span><i class="sw flt"></i>Filtered</span></div><canvas id="cx"></canvas>
<h3>Y position</h3><div class="legend"><span><i class="sw raw"></i>Raw</span><span><i class="sw flt"></i>Filtered</span></div><canvas id="cy"></canvas>
<h3>Filtered − raw residual</h3><div class="legend"><span><i class="sw err"></i>X residual</span><span><i class="sw flt"></i>Y residual</span></div><canvas id="ce"></canvas>
<script>
const D=__DATA__, M=__META__; const root=document.documentElement; let lo=0,hi=Math.max(1,D.length-1),drag=null;
document.getElementById('meta').textContent=`${M.rate} Hz · ${M.windowMs} ms window (${M.windowSamples} samples) · ${M.lookaheadMs.toFixed(3)} ms analysis lag`;
const specs=[['cx',[[1,'#888'],[3,'#2878d0']]],['cy',[[2,'#888'],[4,'#2878d0']]],['ce',[[5,'#d06b28'],[6,'#2878d0']]]];
function setup(c){const dpr=devicePixelRatio||1, r=c.getBoundingClientRect();c.width=Math.max(1,Math.floor(r.width*dpr));c.height=Math.max(1,Math.floor(r.height*dpr));}
function bounds(series){let mn=Infinity,mx=-Infinity;for(let i=Math.floor(lo);i<=Math.ceil(hi)&&i<D.length;i++)for(const [k] of series){const v=D[i][k];if(Number.isFinite(v)){mn=Math.min(mn,v);mx=Math.max(mx,v)}}if(!isFinite(mn)){mn=0;mx=1}if(mx===mn){mx+=1;mn-=1}const p=(mx-mn)*.06;return[mn-p,mx+p]}
function draw(){for(const [id,series] of specs){const c=document.getElementById(id);setup(c);const g=c.getContext('2d'),W=c.width,H=c.height,dpr=devicePixelRatio||1,p=36*dpr;g.clearRect(0,0,W,H);const [mn,mx]=bounds(series);g.strokeStyle='#777';g.lineWidth=1*dpr;g.beginPath();g.moveTo(p,8*dpr);g.lineTo(p,H-p);g.lineTo(W-8*dpr,H-p);g.stroke();for(const [k,col] of series){g.strokeStyle=col;g.lineWidth=1.35*dpr;g.beginPath();let first=true;const a=Math.max(0,Math.floor(lo)),b=Math.min(D.length-1,Math.ceil(hi));for(let i=a;i<=b;i++){const x=p+(i-lo)/(hi-lo)*(W-p-8*dpr),y=(H-p)-(D[i][k]-mn)/(mx-mn)*(H-p-8*dpr);if(first){g.moveTo(x,y);first=false}else g.lineTo(x,y)}g.stroke()}g.fillStyle=getComputedStyle(document.body).color;g.font=`${11*dpr}px system-ui`;g.fillText(mx.toFixed(1),2,14*dpr);g.fillText(mn.toFixed(1),2,H-p);}
 const t0=D[Math.max(0,Math.floor(lo))]?.[0]??0,t1=D[Math.min(D.length-1,Math.ceil(hi))]?.[0]??0;document.getElementById('range').textContent=`${t0.toFixed(2)}–${t1.toFixed(2)} ms`;}
function reset(){lo=0;hi=Math.max(1,D.length-1);draw()} document.getElementById('reset').onclick=reset;
for(const [id] of specs){const c=document.getElementById(id);c.addEventListener('wheel',e=>{e.preventDefault();const r=c.getBoundingClientRect(),q=(e.clientX-r.left)/r.width,span=hi-lo,f=e.deltaY<0?.78:1.28,nspan=Math.max(8,Math.min(D.length-1,span*f)),center=lo+q*span;lo=Math.max(0,center-q*nspan);hi=Math.min(D.length-1,lo+nspan);lo=Math.max(0,hi-nspan);draw()},{passive:false});c.addEventListener('pointerdown',e=>{c.setPointerCapture(e.pointerId);drag=[e.clientX,lo,hi]});c.addEventListener('pointermove',e=>{if(!drag)return;const span=drag[2]-drag[1],dx=(e.clientX-drag[0])/c.getBoundingClientRect().width*span;lo=Math.max(0,drag[1]-dx);hi=Math.min(D.length-1,drag[2]-dx);if(hi-lo<span){if(lo===0)hi=span;else lo=hi-span}draw()});c.addEventListener('pointerup',()=>drag=null);c.addEventListener('pointercancel',()=>drag=null);c.addEventListener('dblclick',reset)}
addEventListener('resize',draw);draw();
</script></body></html>
""";
        html = html.Replace("__DATA__", data).Replace("__META__", meta);
        File.WriteAllText(path,html,new UTF8Encoding(false));
    }
}
