# AA The Perfect Filter — OpenTabletDriver plugin

The plugin runs on OTD's fixed 1000 Hz output scheduler while estimating the tablet input rate independently. Filtering is performed in uniformly spaced sample-index space; report arrival timestamps are used for rate estimation and stream-reset detection rather than as acquisition timestamps.

## Defaults

- Window duration: 100 ms (recommended 70–150 ms)
- Lambda cutoff: 1.5 when adaptive lambda is disabled
- Depth / maximum finite-difference order: 8 (recommended 8)
- Outlier threshold: 6 (recommended 5–8)
- Latency: 5 ms (recommended 0–20 ms)
- Adaptive lambda: enabled
- Adaptive anchors: 0 ms 1.50; 2 ms 1.20; 5 ms 1.05; 10 ms 1.04; 20 ms 1.03
- CSV logging: disabled

Adaptive lambda linearly interpolates between the latency anchors. At 20 ms and above the 20 ms value is used. Turning adaptive lambda off makes the normal Lambda cutoff setting authoritative.

Latency is clamped to 0–20 ms. The fixed 1000 Hz scheduler interpolates within the reconstructed trajectory when buffered samples are available and linearly extrapolates the latest reconstructed segment when the scheduler runs ahead of the newest tablet report. Extrapolation is automatic and is not exposed as a separate option.

Model construction occurs off the realtime path and the previous model remains active until its replacement is ready. CSV logging is intended for diagnostics only.
