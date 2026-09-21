# AA The Perfect Filter — OpenTabletDriver plugin

The plugin runs on OTD's fixed 1000 Hz output scheduler while estimating the tablet input rate independently. Filtering is performed in uniformly spaced sample-index space; report arrival timestamps are used for rate estimation and stream-reset detection rather than as acquisition timestamps.

## Defaults

- Window duration: 100 ms (recommended 70–150 ms)
- Filter strength: 1.0 (0 retains all modes; 1 reproduces the original cutoffs)
- Outlier threshold: 6 (recommended 5–8)
- Latency: 5 ms (recommended 0–20 ms)
- Staleness timeout: 25 ms
- CSV logging: disabled

Adaptive lambda uses the original latency anchors: 0 ms 1.50; 2 ms 1.20; 5 ms 1.05; 10 ms 1.04; 20 ms 1.03. Filter strength scales the distance of each cutoff above 1 inversely: strength 1 leaves the ramp unchanged, strength 0 retains all modes, and higher values exclude more modes. Adaptive cutoffs have a hard minimum of 1.03. Disabling adaptive lambda applies the same strength mapping to a fixed lambda of 1.50.

The maximum difference order remains fixed internally at 8 because it had no useful user-facing effect. The experimental future-boundary fit is not used by the production filter.

Latency is clamped to 0–20 ms. Zero latency remains fully functional but is not recommended because it reduces noise rejection and tolerance of input-timing jitter. The fixed 1000 Hz scheduler interpolates within the reconstructed trajectory. When it runs ahead of the newest report, an even-boundary continuation is used through the DCT-II reflection boundary half a sample later; only requests beyond that point use the linear delayed-input fallback. The scheduler stops emitting after the configured staleness timeout and resumes on the next physical report, without relying on a tablet-specific out-of-range report.

Model construction occurs off the realtime path and the previous model remains active until its replacement is ready. CSV logging is intended for diagnostics only.
