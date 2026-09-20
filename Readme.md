# Mathematically Optimal Noise Rejection (Anti-chatter) Filter & Resampler
This is an Open Tablet Driver plugin based on spectral filtering. Noise has roughly uniform spectral power density but intended motion does not. 

- Up to ~6dB noise reduction at 0 latency **without** any smoothing, averaging, or deadzones. 
- ~10dB noise reduction with 2 sample buffer (compared to 4.77dB from a typical smoothing filter).
- Fully configurable filter strength and latency
- Reconstructive resampler (not interpolation), creates accurate trajectories without requiring a buffer or prediction. 

## Usage
Use with Open Tablet Driver 0.6.7 or newer, place the [plugin](https://github.com/Vestaia/ThePerfectFilter/releases/download/v1.0.0/AA.ThePerfectFilter.dll) at `C:\Users\<your username>\AppData\Local\OpenTabletDriver\Plugins\` and enable in your filters menu.

## Parameters
The defaults are optimized for 700hz custom Wacom firmware. I would recommend using the custom firmware if available for your tablet.
| Parameter | Symbol | Default | Effect on output |
|---|---:|---:|---|
| Window duration | $T_w$ | 100 ms | Sets the amount of position history used to construct the reconstruction. Longer windows provide more temporal context and finer modal resolution, while shorter windows make the model more local in time. |
| Maximum derivative order | $M$ | 8 | Sets maximum order for the finite differences operators used for constructing DCT eigenbasis. Maximum order influences eigenvalues. No reason to change this. |
| Lambda cutoff | $\lambda_{\mathrm{fixed}}$ | 1.50 | Sets the maximum eigenvalue retained when adaptive lambda is disabled. Lower values retain fewer modes and reject more noise but increase reconstruction error and implicit delay. Higher values retain more modes, improving endpoint tracking at the cost of admitting more noise. |
| Latency | $T_L$ | 5 ms | Sets how long output is delayed so that later samples can contribute to reconstruction of the reported position. Increasing latency generally permits substantially stronger noise rejection for the same trajectory accuracy. Zero latency forces reconstruction at the newest available sample. |
| Staleness timeout | $T_S$ | 25 ms | Stops emitting output when no physical tablet report has arrived for this long. This prevents continued extrapolation when a tablet leaves its last in-range report cached after the pen is lifted. Output resumes with the next physical report. |
| Outlier threshold | $Z_{\mathrm{outlier}}$ | 6 | Sets the threshold for rejecting isolated position deviations classified as outliers. Lower values reject smaller deviations more aggressively; higher values restrict replacement to more extreme deviations. |
| Use adaptive lambda | $A$ | True | Selects whether the lambda cutoff changes with reconstruction latency. When enabled, low-latency estimates use higher cutoffs to reduce endpoint error while estimates with more future information use lower cutoffs for stronger noise rejection. When disabled, all estimates use $\lambda_{\mathrm{fixed}}$. |
| Adaptive lambda, 0 ms | $\lambda_0^{A}$ | 1.50 | Sets the eigenvalue cutoff for reconstruction at the newest sample, where no future samples are available. This is the least aggressively filtered adaptive estimate because additional modes are required to reduce endpoint error. |
| Adaptive lambda, 2 ms | $\lambda_2^{A}$ | 1.20 | Sets the adaptive cutoff at 2 ms of reconstruction latency. Cutoffs at intermediate latencies are interpolated between adjacent adaptive-lambda parameters. |
| Adaptive lambda, 5 ms | $\lambda_5^{A}$ | 1.05 | Sets the adaptive cutoff at 5 ms of reconstruction latency. The additional future information allows a substantially lower-rank reconstruction than at the newest sample. |
| Adaptive lambda, 10 ms | $\lambda_{10}^{A}$ | 1.04 | Sets the adaptive cutoff at 10 ms of reconstruction latency. The small reduction relative to the 5 ms cutoff reflects the diminishing benefit of additional look-ahead at this timescale. |
| Adaptive lambda, 20 ms | $\lambda_{20}^{A}$ | 1.03 | Sets the adaptive cutoff at 20 ms of reconstruction latency. This is the lowest default adaptive cutoff and therefore the most aggressively filtered reconstruction in the default schedule. |
| CSV logging | — | False | Writes filter diagnostics and position data to CSV when enabled. It does not intentionally alter the reconstructed output and is disabled by default to avoid unnecessary I/O overhead. |
## Filter Design

```math
\mathbf{x}_t =
\begin{bmatrix}
x_t\\
y_t
\end{bmatrix}
```

```math
N =
\mathrm{round}\!\left(
\frac{f_s T_w}{1000}
\right)
```

```math
\mathbf{x}^{(N)}_t =
\begin{bmatrix}
\mathbf{x}_{t-N+1} &
\mathbf{x}_{t-N+2} &
\cdots &
\mathbf{x}_{t}
\end{bmatrix}
```

```math
D_m =
\begin{bmatrix}
\mathbf{d}_{m,0}^{\mathsf T}\\
\mathbf{d}_{m,1}^{\mathsf T}\\
\vdots\\
\mathbf{d}_{m,N-1}^{\mathsf T}
\end{bmatrix}
```

```math
\sum_{j\in S_i}
d_{m,i,j}(j-i)^p
=
\begin{cases}
m!, & p=m\\
0, & p\ne m
\end{cases}
\qquad
p=0,1,\ldots,m+1
```

```math
|S_i|=m+2
```

```math
\mathbf{d}_{m,i}
\leftarrow
\frac{\mathbf{d}_{m,i}}
{\left\|\mathbf{d}_{m,i}\right\|_2}
```

```math
L_m=D_m^{\mathsf T}D_m
```

```math
L_m
\leftarrow
\frac{L_m}
{\lambda_{\max}(L_m)}
```

```math
L =
\sum_{m=0}^{M}
\alpha_m L_m
```

```math
L\mathbf{q}_k
=
\lambda_k\mathbf{q}_k
```

```math
\mathbf{q}_i^{\mathsf T}\mathbf{q}_j
=
\delta_{ij}
```

```math
\lambda_0
\le
\lambda_1
\le
\cdots
\le
\lambda_{N-1}
```

```math
Q =
\begin{bmatrix}
\mathbf{q}_0 &
\mathbf{q}_1 &
\cdots &
\mathbf{q}_{N-1}
\end{bmatrix}
```

```math
Q^{\mathsf T}Q=I
```

```math
g_k(\lambda_c)
=
\begin{cases}
1, & \lambda_k\le\lambda_c\\
0, & \lambda_k>\lambda_c
\end{cases}
```

```math
G(\lambda_c)
=
\mathrm{diag}
\left(
g_0,g_1,\ldots,g_{N-1}
\right)
```

```math
P_{\lambda_c}
=
QG(\lambda_c)Q^{\mathsf T}
```

```math
\hat{\mathbf{x}}^{(N)}_t
=
P_{\lambda_c}\mathbf{x}^{(N)}_t
```

```math
\hat{x}_i
=
\sum_{k=0}^{N-1}
g_k q_k[i]
\sum_{j=0}^{N-1}
q_k[j]x_j
```

```math
\hat{y}_i
=
\sum_{k=0}^{N-1}
g_k q_k[i]
\sum_{j=0}^{N-1}
q_k[j]y_j
```

```math
h_{\lambda_c,i}[j]
=
\sum_{k=0}^{N-1}
g_k q_k[i]q_k[j]
```

```math
\hat{x}_i
=
\sum_{j=0}^{N-1}
h_{\lambda_c,i}[j]x_j
```

```math
\hat{y}_i
=
\sum_{j=0}^{N-1}
h_{\lambda_c,i}[j]y_j
```

```math
V_{\lambda_c,i}
=
\left\|
\mathbf{h}_{\lambda_c,i}
\right\|_2^2
```

```math
\frac{\sigma_{\mathrm{out}}^2}
{\sigma_{\mathrm{in}}^2}
=
V_{\lambda_c,i}
=
\sum_{k=0}^{N-1}
g_k^2 q_k[i]^2
```

```math
R_{\mathrm{dB}}
=
-10\log_{10}
\left(
V_{\lambda_c,i}
\right)
```

```math
\tau_r =
\frac{1000r}{f_s}
```

```math
\lambda_c(\tau)
=
\lambda_a
+
\frac{\tau-\tau_a}{\tau_b-\tau_a}
\left(
\lambda_b-\lambda_a
\right),
\qquad
\tau_a\le\tau\le\tau_b
```

```math
(\tau_a,\lambda_a),(\tau_b,\lambda_b)
\in
\left\{
(0,\lambda_0^{A}),
(2,\lambda_2^{A}),
(5,\lambda_5^{A}),
(10,\lambda_{10}^{A}),
(20,\lambda_{20}^{A})
\right\}
```

```math
\lambda_c(\tau)
=
\lambda_{20}^{A},
\qquad
\tau\ge20
```

```math
\lambda_c(\tau)
=
\begin{cases}
\lambda_{\mathrm{adaptive}}(\tau),
& A=1\\
\lambda_{\mathrm{fixed}},
& A=0
\end{cases}
```

```math
i_r=N-1-r
```

```math
\mathbf{h}_r
=
\mathbf{h}_{\lambda_c(\tau_r),i_r}
```

```math
V_r
=
\left\|\mathbf{h}_r\right\|_2^2
```

```math
\tilde{w}_r
=
\frac{1}{V_r}
```

```math
w_r
=
\frac{V_r^{-1}}
{\displaystyle\sum_{s=0}^{R}V_s^{-1}}
```

```math
\hat{\mathbf{x}}_t^{(r)}
=
\sum_{j=0}^{N-1}
h_r[j]\,
\mathbf{x}_{t-r-N+1+j}
```

```math
\hat{\mathbf{x}}_t
=
\sum_{r=0}^{R}
w_r
\hat{\mathbf{x}}_t^{(r)}
```

```math
R =
\mathrm{round}\!\left(
\frac{f_s T_L}{1000}
\right)
```

```math
T_L\ge0
```

```math
z_t
=
\frac{
\left\|
\mathbf{x}_t-\hat{\mathbf{x}}_t
\right\|_2
}
{\sigma_t}
```

```math
\mathbf{x}^{*}_t
=
\begin{cases}
\hat{\mathbf{x}}_t,
& z_t>Z_{\mathrm{outlier}}\\
\mathbf{x}_t,
& z_t\le Z_{\mathrm{outlier}}
\end{cases}
```


