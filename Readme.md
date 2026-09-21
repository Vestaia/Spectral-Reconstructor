# The Perfect Filter

The Perfect Filter is an OpenTabletDriver trajectory filter and resampler. It reconstructs tablet motion in a curvature-ordered eigenbasis, retaining low-complexity motion while attenuating high-complexity input noise.

- Preserves constant-position and constant-velocity trajectories in its lowest eigenspace.
- Supports configurable filtering strength and output latency.
- Runs on OpenTabletDriver's 1000 Hz asynchronous output pipeline.
- Uses linear extrapolation only when the output playback head is ahead of the latest tablet report.

## Usage
Use with OpenTabletDriver 0.6.7 or newer. Download the plugin from the [latest release](https://github.com/Vestaia/ThePerfectFilter/releases/latest), place it in `C:\Users\<your username>\AppData\Local\OpenTabletDriver\Plugins\`, and enable it in the Filters settings.

## License

This project is licensed under the [Sustainable Use License 1.0](LICENSE.md) (`SUL-1.0`). It may be used and modified for personal, noncommercial, or internal business purposes. Redistribution must be free of charge and for noncommercial purposes.

## Parameters
The defaults are optimized for 700hz custom Wacom firmware. I would recommend using the custom firmware if available for your tablet.
| Parameter | Symbol | Default | Effect on output |
|---|---:|---:|---|
| Window duration | $T_w$ | 100 ms | Sets the amount of position history used to construct the reconstruction. Longer windows provide more temporal context and finer modal resolution, while shorter windows make the model more local in time. |
| Filter strength | $s$ | 1 | Inversely scales retained modes above the affine minimum. Strength 1 uses 12 modes at zero look-ahead, 0 retains all modes, and larger values retain fewer modes. |
| Latency | $T_L$ | 5 ms | Sets how long output is delayed so that later samples can contribute to reconstruction of the reported position. Increasing latency generally permits substantially stronger noise rejection for the same trajectory accuracy. Zero latency forces reconstruction at the newest available sample. |
| Staleness timeout | $T_S$ | 25 ms | Stops emitting output when no physical tablet report has arrived for this long. This prevents continued extrapolation when a tablet leaves its last in-range report cached after the pen is lifted. Output resumes with the next physical report. |
| Outlier threshold | $Z_{\mathrm{outlier}}$ | 6 | Sets the threshold for rejecting isolated position deviations classified as outliers. Lower values reject smaller deviations more aggressively; higher values restrict replacement to more extreme deviations. |
| Use adaptive modes | $A$ | True | Selects whether retained mode count follows the sample-look-ahead schedule 12, 9, 8, 7, 6, then 5. The schedule is scaled by Filter strength. |
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
L = L_0 + \sum_{m=2}^{M} L_m
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
g_k(K)
=
\begin{cases}
1, & 0\le k<K\\
0, & K\le k<N
\end{cases}
```

```math
G(K)
=
\mathrm{diag}
\left(
g_0,g_1,\ldots,g_{N-1}
\right)
```

```math
P_K
=
QG(K)Q^{\mathsf T}
```

```math
\hat{\mathbf{x}}^{(N)}_t
=
P_K\mathbf{x}^{(N)}_t
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
h_{K,i}[j]
=
\sum_{k=0}^{N-1}
g_k q_k[i]q_k[j]
```

```math
\hat{x}_i
=
\sum_{j=0}^{N-1}
h_{K,i}[j]x_j
```

```math
\hat{y}_i
=
\sum_{j=0}^{N-1}
h_{K,i}[j]y_j
```

```math
V_{K,i}
=
\left\|
\mathbf{h}_{K,i}
\right\|_2^2
```

```math
\frac{\sigma_{\mathrm{out}}^2}
{\sigma_{\mathrm{in}}^2}
=
V_{K,i}
=
\sum_{k=0}^{N-1}
g_k^2 q_k[i]^2
```

```math
R_{\mathrm{dB}}
=
-10\log_{10}
\left(
V_{K,i}
\right)
```

```math
\tau_r =
\frac{1000r}{f_s}
```

```math
K_{\mathrm{base}}(r)
=
\begin{cases}
12,&r=0\\
9,&r=1\\
8,&r=2\\
7,&r=3\\
6,&r=4\\
5,&r\ge5
\end{cases}
```

```math
K_s(r)
=
2+\frac{K_{\mathrm{base}}(r)-2}{s^2}
```

```math
K(r)
=
\begin{cases}
N,&s=0\\
\operatorname{clamp}_{[2,N]}\!\left(\operatorname{round}(K_s(r))\right),&s>0
\end{cases}
```

```math
i_r=N-1-r
```

```math
\mathbf{h}_r
=
\mathbf{h}_{K(r),i_r}
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


