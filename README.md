# SMD

MATLAB toolbox for single-molecule localization and tracking in fluorescence microscopy data.

## Features

- **Two detection pipelines**:
  - **LLR** (default) — Log-likelihood ratio test + radial symmetry / integrated Gaussian fitting
  - **Wavelet** — a-trous wavelet filter + Crocker-Grier peak finding + Poisson MLE Gaussian fitting
- **Two noise models** for the LLR pipeline:
  - `'gaussian'` — least-squares integrated Gaussian (fast, suitable for high-SNR data)
  - `'poisson'` — radial symmetry initial guess + Poisson MLE Gaussian refinement
- Particle tracking via `quot` tracker (Hungarian linker)
- ROI-based trajectory culling with watershed segmentation
- Drift correction via maximum spanning forest
- ND2 and TIFF image stack support
- Visualization of tracks, bright-field images, and intensity projections

## Setup

Add the toolbox and all helper directories to the MATLAB path:

```matlab
smd_root = '/path/to/SMD';
addpath(smd_root);
addpath(fullfile(smd_root, 'helpers', 'crocker_grier'));
addpath(fullfile(smd_root, 'helpers', 'simpletracker'));
addpath(fullfile(smd_root, 'helpers', 'drift'));
```

Or add all at once:

```matlab
addpath(genpath('/path/to/SMD'));
```

## Quick Start

```matlab
% Create an SMD object
smd = SMD('/path/to/data/', 'movie.tif', 0.02, 50);

% --- Default pipeline: LLR detection + quot tracker ---
smd.localize();
smd.track();
smd.all_plot();

% --- Option: Wavelet pipeline ---
smd.detection_method = 'wavelet';
smd.localize();
smd.track();
smd.all_plot();

% --- LLR tuning (optional) ---
smd.noise_model = 'poisson';    % radialcenter + Poisson MLE (more accurate, slower)
smd.psf_sigma = 1.3;            % PSF width in pixels
smd.llr_params.t = 20.0;        % detection score threshold
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
| `detection_method` | `'llr'` | `'llr'` or `'wavelet'` |
| `noise_model` | `'gaussian'` | `'gaussian'` or `'poisson'` (LLR only) |
| `psf_sigma` | `1.3` | PSF standard deviation (pixels) |
| `llr_params` | `struct('k',1,'w',9,'t',20)` | LLR kernel width, window, threshold |
| `tracking_params` | `struct(...)` | Tracking radius, gap frames, min length, motion type, tracker (`'quot'` or `'simpletracker'`) |

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

## Repository Structure

```
SMD/
├── @SMD/                       Core class definition and methods
│   ├── SMD.m                   Class definition, constructor, getters
│   ├── localize.m              Spot detection + fitting (wavelet & LLR)
│   ├── track.m                 Trajectory linking
│   ├── cull_tracks.m           ROI-based track filtering
│   ├── get_roi.m               ROI detection / manual drawing
│   ├── image_regions.m         Bright-field image loading
│   ├── remove_relative_motion.m  Drift correction
│   └── all_plot.m              Visualization
├── helpers/
│   ├── crocker_grier/          Peak finding & wavelet detection
│   │   ├── pkfnd.m             Local maxima finder (Crocker & Grier)
│   │   ├── cntrd.m             Sub-pixel centroid refinement
│   │   └── wavelet_filter.m    A-trous wavelet filter (Izeddin et al. 2012)
│   ├── simpletracker/          Frame-to-frame particle linking
│   │   ├── simpletracker.m     Main tracker (J.-Y. Tinevez, 2011)
│   │   ├── hungarianlinker.m   Hungarian algorithm linker
│   │   └── nearestneighborlinker.m  Nearest-neighbor gap closing
│   └── drift/                  Drift correction utilities
│       ├── maxSpanningForest.m
│       └── relativeTracksFromForest.m
├── calcBG.m                    Local background estimation (LMS filter)
├── gaussfit2DMLE.m             Poisson MLE 2D Gaussian fitting (Parthasarathy)
├── tiffread2.m                 Multi-frame TIFF reader
├── llr.m                       Log-likelihood ratio spot detection
├── mle_amp_setup.m             LLR kernel setup
├── ls_int_gaussian.m           Integrated Gaussian fitting (Levenberg-Marquardt)
├── radialcenter.m              Radial symmetry center finding (Parthasarathy, GPL)
└── README.md
```

## Dependencies

- MATLAB R2019b+ (for `arguments` blocks)
- Image Processing Toolbox (`imtophat`, `wiener2`, `imgaussfilt`, `imregionalmax`, `strel`)
- [Bio-Formats for MATLAB](https://www.openmicroscopy.org/bio-formats/) — required for reading ND2 files (provides `BioformatsImage` and `getPlane`). Not needed for TIFF-only workflows.

All other dependencies are included in `helpers/`.

## Third-Party Acknowledgments

| Component | Author | License |
|-----------|--------|---------|
| `radialcenter.m` | Raghuveer Parthasarathy (U. Oregon, 2011-2012) | GPL v3 |
| `gaussfit2DMLE.m` | Raghuveer Parthasarathy (2012) | — |
| `simpletracker` | Jean-Yves Tinevez (2011-2012) | — |
| `pkfnd.m`, `cntrd.m` | Eric Dufresne, Daniel Blair, John Crocker | — |
| `wavelet_filter.m` | Based on Izeddin et al., *Opt. Express* 20, 2081-2095 (2012) | — |

## References

- Radial symmetry center: Parthasarathy, *Nature Methods* 9, 724–726 (2012)
- Gaussian MLE fitting: Abraham et al., *Opt. Express* 17, 23352 (2009)
- Log-likelihood ratio: Serge et al., *Nature Methods* 5, 687–694 (2008)
- Wavelet detection: Izeddin et al., *Opt. Express* 20, 2081–2095 (2012)
- Crocker-Grier algorithm: Crocker & Grier, *J. Colloid Interface Sci.* 179, 298–310 (1996)

## License

MIT
