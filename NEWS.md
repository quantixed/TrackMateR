# TrackMateR

# TrackMateR 0.3.11

- Better handling of more than one channel in TrackMate XML files

# TrackMateR 0.3.10

- Fixed an issue where sparse TrackMate XML files caused `compareDatasets()` to fail
- These files are skipped and a warning is issued. They can be loaded and analyzed manually if desired (limit is more than 3 tracks with at least one with 10 frames)

# TrackMateR 0.3.9

- Warning for TrackMate XML files that contain tracks with more than one point per frame
- Added xy scaling for shape descriptors

# TrackMateR 0.3.8

- Added report of intensity and duration to workflow

# TrackMateR 0.3.7

- Added estimator of diffusion co-efficient and histogram of D (per track)

# TrackMateR 0.3.6

- Ability to load a ground truth csv dataset
- Compare multiple ground truth datasets

# TrackMateR 0.3.5

- More flexibility to automated workflows (axes on MSD plots can be configured)
- Minor fixes to plots in `compareDatasets()`
- Better handling in jump distance calculation
- Test datasets available for multiple comparisons
- Improved documentation

# TrackMateR 0.3.4

- Fixes for parallelisation. Windows users still do not benefit from parallel processing that Mac/Linux users enjoy.

# TrackMateR 0.3.3

- More flexibility to automated workflows by allowing parameters to be passed via ellipsis to `reportDataset()` and `compareDatasets()`
- Better resilience to failures in jump distance fitting, and reading XML

# TrackMateR 0.3.2

- Fixes for edge cases and issues

# TrackMateR 0.3.1

- Fixes for parallelisation issue and better handling of fitting errors

# TrackMateR 0.3.0

- First release on GitHub
