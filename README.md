# sr_tracking

MATLAB toolbox for single-molecule localization and tracking in fluorescence microscopy data.

## Features

- **Two detection pipelines**:
  - **Wavelet** (default) — a-trous wavelet filter + Crocker-Grier peak finding + Poisson MLE Gaussian fitting
  - **LLR** — Log-likelihood ratio test + radial symmetry / integrated Gaussian fitting
- **Two noise models** for the LLR pipeline:
  - `'gaussian'` — least-squares integrated Gaussian (fast, suitable for high-SNR data)
  - `'poisson'` — radial symmetry initial guess + Poisson MLE Gaussian refinement
- Particle tracking via `simpletracker`
- ROI-based trajectory culling with watershed segmentation
- Drift correction via maximum spanning forest
- Visualization of tracks, bright-field images, and intensity projections

## Quick Start

```matlab
% Create an SMD object
smd = SMD('/path/to/data/', 'movie.tif', 0.02, 50);

% --- Option A: Wavelet pipeline (default) ---
smd.localize();
smd.track();
smd.all_plot();

% --- Option B: LLR pipeline ---
smd.detection_method = 'llr';
smd.noise_model = 'gaussian';   % or 'poisson'
smd.psf_sigma = 1.3;            % PSF width in pixels
smd.localize();
smd.track();
smd.all_plot();
```

## Class Properties

### User-Configurable

| Property | Default | Description |
|----------|---------|-------------|
| `dirname` | — | Parent directory for raw data |
| `fname` | — | Image stack filename (.tif or .nd2) |
| `camera` | `'EMCCD'` | Camera type |
| `pixelsize` | `0.16` | Pixel size in microns |
| `exposure_time` | — | Exposure time (seconds) |
| `frame_rate` | — | Acquisition frame rate (Hz) |
| `localization_box` | `3` | Half-width of fitting sub-image (pixels) |
| `localization_threshold` | `1.75` | Wavelet detection threshold |
| `detection_method` | `'wavelet'` | `'wavelet'` or `'llr'` |
| `noise_model` | `'gaussian'` | `'gaussian'` or `'poisson'` (LLR only) |
| `psf_sigma` | `1.3` | PSF standard deviation (pixels) |
| `llr_params` | `struct('k',1,'w',9,'t',20)` | LLR kernel width, window, threshold |
| `tracking_params` | `struct(...)` | Tracking radius, gap frames, min length, motion type |

### Read-Only (use getter methods)

| Property | Getter | Description |
|----------|--------|-------------|
| `spots` | `get_spots()` | Cell array of localized spots per frame |
| `tracks` | `get_tracks()` | Cell array of trajectories |
| `relative_tracks` | `get_relative_tracks()` | Drift-corrected trajectories |
| `maxint` | `get_maxint()` | Maximum intensity projection |

## Pipeline Methods

| Method | Description |
|--------|-------------|
| `localize()` | Detect and fit spots in image stack |
| `track()` | Link spots into trajectories |
| `cull_tracks()` | Filter trajectories by ROI membership |
| `get_roi()` | Detect or draw ROIs (nuclear boundaries) |
| `image_regions()` | Load bright-field reference images |
| `remove_relative_motion()` | Drift correction via spanning forest |
| `all_plot()` | 2x2 visualization panel |

## Dependencies

- MATLAB R2019b+ (for `arguments` blocks)
- Image Processing Toolbox (`imtophat`, `wiener2`, `imgaussfilt`, `imregionalmax`, `strel`)
- `simpletracker` (included or on path)
- `tiffread2` (included)
- `wavelet_filter`, `pkfnd`, `cntrd` (a-trous / Crocker-Grier utilities)

## References

- Radial symmetry center: Parthasarathy, *Nature Methods* 9, 724–726 (2012)
- Gaussian MLE fitting: Abraham et al., *Opt. Express* 17, 23352 (2009)
- Log-likelihood ratio: based on Serge et al., *Nature Methods* 5, 687–694 (2008)

## License

MIT
