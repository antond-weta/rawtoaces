This branch contains a proposed version 2 of rawtoaces.

See [here](v1/README.md) for the original README.

As this is pretty much a complete rewrite, it is not backwards compatible.
- Both libraries' public interfaces have been changed.
- The command line tool's command line parameters have been changed.
Read below for more info.

# Main changes in this branch:
- Removes dependency on boost:
  - boost::filesystem is replaced by std::filesystem
  - boost::json is replaced by nlohmann_json
  - boost::unittest is replaced by OpenImageIO
- Replaces libraw with OpenImageIO for reading raw files.
- Replaces AcesContainer with OpenImageIO for writing .exr files.
- Replaces the custom command line parser with OpenImageIO.
- Implements a custom IDT matrix solver, making the dependency on Eigen and Ceres_solver optional.

These changes are (probably) temporary, to make the transition easier:
- The targets are renamed, so both v1 and v2 can be built simultaneously:
  - rawtoaces_idt -> raw2aces_core
  - rawtoaces_util -> raw2aces_util
  - rawtoaces -> raw2aces
- The new version uses "_rta2" instead of "_aces" as the output file suffix.

Since this version uses a different command line parser, there are some changes to the command line parameters, so the new version uses slightly different syntax. See below for some examples of how to use different white balance and matrix methods:

## White balance methods:
| old | new |
|-----|-----|
| rawtoaces --wb-method 0            | raw2aces --wb-method metadata                      |
| rawtoaces --wb-method 1 d65        | raw2aces --wb-method illuminant --illuminant d65   |
| rawtoaces --wb-method 2            | raw2aces --wb-method box                           |
| rawtoaces --wb-method 3 \<x y w h> | raw2aces --wb-method box --wb-box \<x y w h>       |
| rawtoaces --wb-method 4 \<r g b g> | raw2aces --wb-method custom --custom-wb \<r g b g> |

## Matrix methods:
| old | new |
|-----|-----|
| rawtoaces --mat-method 0                                            | raw2aces --mat-method spectral                                                   |
| rawtoaces --mat-method 1                                            | raw2aces --mat-method metadata                                                   |
| rawtoaces --mat-method 2                                            | raw2aces --mat-method Adobe                                                      |
| rawtoaces --mat-method 3 \<m1r m1g m1b m2r m2g m2b m3r m3g m3b>     | raw2aces --mat-method custom --custom-mat \<m1r m1g m1b m2r m2g m2b m3r m3g m3b> |

## Other changes
- raw2aces will not overwrite output files unless executed with "--overwrite"
- custom output directory can be specified with "--output-dir"
- missing output directories provided in "--output-dir" will be created automatically if executed with "--create-dirs"
- uses RAWTOACES_DATA_PATH environment variable instead of AMPAS_DATA_PATH

# Build

## Minimal build
To build only the new version wih minimal dependencies, install OpenImageIO and nlohmann-json, and configure like this:

    cmake \
    -D RTA_ENABLE_EIGEN=0 \
    -D RTA_ENABLE_CERES=0 \
    -D RTA_BUILD_V1=0

## Full build
To build both the new and original versions, install the new dependencies (OpenImageIO, nlohmann-json), and follow the 
normal build instructions [here](v1/README.md) 


