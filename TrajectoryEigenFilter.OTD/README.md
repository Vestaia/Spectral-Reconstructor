# Spectral Reconstructor – OpenTabletDriver plugin

Spectral Reconstructor rejects tablet noise without conventional smoothing and reconstructs the 1000 Hz output trajectory rather than merely interpolating between raw reports. Its curvature-ordered spectral model preserves simple intentional motion while suppressing high-complexity input noise, providing strong noise reduction with responsive tracking.

The plugin runs on OTD's fixed 1000 Hz output scheduler while estimating the tablet input rate independently. Filtering is performed in uniformly spaced sample-index space; report arrival timestamps are used for rate estimation and stream-reset detection rather than as acquisition timestamps.

## Defaults

- Window duration: 100 ms (recommended 70–150 ms)
- Filter strength: 1 (12 modes at zero look-ahead; 0 retains all modes; 1.5 is recommended for 133 Hz tablets)
- Adaptive modes: enabled; sample-look-ahead schedule 12, 9, 8, 7, 6, then 5
- Outlier threshold: 6 (recommended 5–8)
- Latency: 4 input samples (recommended 1–10 samples)
- Staleness timeout: 25 ms
- CSV logging: disabled

Adaptive mode counts are scaled as `2 + (base modes - 2) / strength^2`, rounded deterministically to the nearest integer, and clamped from 2 through the window length. Strength 1 reproduces the base schedule; strength 0 retains all modes. Eigenvalues remain available for diagnostics and ordering but are not reconstruction thresholds.

The maximum difference order is fixed internally at 8 because it has no useful user-facing control effect.

Latency is clamped to 0–20 input samples. Zero latency remains fully functional but is not recommended because it reduces noise rejection and tolerance of input-timing jitter. The fixed 1000 Hz scheduler interpolates within the reconstructed trajectory and linearly extrapolates when its playback head is ahead of the newest report. The scheduler stops emitting after the configured staleness timeout and resumes on the next physical report, without relying on a tablet-specific out-of-range report.

Model construction occurs off the realtime path and the previous model remains active until its replacement is ready. CSV logging is intended for diagnostics only.
