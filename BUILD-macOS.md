# Building Structure-from-Sherds++ on macOS (Apple Silicon)

Upstream targets Linux via Docker (or Windows via vcpkg) and pins Ceres 1.14 and
PCL 1.9. This document covers building against current Homebrew packages on
macOS. Verified on macOS 15.6.1 / arm64 with AppleClang.

Note the licence: SfS++ is **CC-BY-NC-SA-4.0**, i.e. non-commercial use only and
share-alike for derivatives. The ICCV 2021 predecessor is MIT.

## 1. Dependencies

```sh
brew install cmake pcl ceres-solver
```

Verified with `pcl 1.15.1`, `ceres-solver 2.2.0`, `vtk 9.7.0`, `libomp 20.1.2`.
The bundled `ceres-solver-1.14.0.tar.gz` is only needed by the Dockerfile and is
not used here.

## 2. Build

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(sysctl -n hw.ncpu)"
```

Produces `build/SfSpp`.

### Port changes

Identical in kind to the ICCV-2021 repository - the same four MSVC/Linux-isms:

| File | Problem | Fix |
|---|---|---|
| `CMakeLists.txt` | Target/glob wiring assumed in-tree Ceres and missed `class/*.cpp`; AppleClang ships no OpenMP runtime. | Single `SfSpp` target, corrected globs, Homebrew `libomp` hints. |
| `main.cpp` | `<experimental/filesystem>` is gone from libc++. | `<filesystem>`. |
| `class/axis_estimation.cpp` | Ceres 2.2 removed `LocalParameterization` / `HomogeneousVectorParameterization` / `SetParameterization`. | Version-guarded `SFS_SET_UNIT_SPHERE` using `SphereManifold<3>`. |
| `class/ranking_system.h` | `StateManager` takes `string&` but is handed a temporary. | `const string&`. |

## 3. SFS_AUTOSAVE

PCL dispatches keyboard shortcuts by VTK key *symbol*, which follows the active
keyboard layout. Under a non-Latin layout (Cyrillic, Greek, ...) no letter key
matches, so `s` never saves and `o`/`f`/`g` never fire - only `spacebar` works,
because it is not a letter.

```sh
cd build && SFS_AUTOSAVE=exit ./SfSpp    # save the rank-1 result, then quit
cd build && SFS_AUTOSAVE=1 ./SfSpp       # save it, then stay in the viewer
```

## 4. Data

SfS++ expects the same four preprocessed inputs per fragment as its predecessor
(breakline, axis, inner/outer surfaces, mesh). Unlike the 2021 release, the code
that produces them **is** published, in a separate repository:
<https://github.com/DominicoRyu/SfSpp_preprocessing>.

`download.sh` fetches the authors' 142-fragment dataset via `gdown`; paths are
selected by the `#define POT_*` block at the top of `class/data_path.h` and
rooted at `path` in the same file.

Our own scans live in `data/test_fragments/` and are gitignored.

## 5. Running our own scans end to end

The preprocessing repository (ported alongside, see its own `BUILD-macOS.md`)
produces surfaces, breaklines and axes; `tools/extract_axis` replaces its
MATLAB-only axis stage. `SFS_DATA_ROOT` points this binary at the result, and
`#define POT_TEST` in `class/data_path.h` selects our four fragments.

```sh
# 1. surfaces  (parameters rescaled for ~550 mm sherds)
cd ../SfSpp_preprocessing/build
SFSPP_SAMPLING_RADIUS=2.5 SFSPP_NORMAL_NEIGHBORS=16 SFSPP_SMOOTHNESS_DEG=7 ./MeshPreprocessing

# 2. axes  (must precede EdgeLineExtraction: rim classification needs them)
cd ../../structure-from-sherds-pp/build
for i in 01 02 03 04; do
  ./extract_axis ../../SfSpp_preprocessing/Temp/Data/Pot_A/Pot_A_Piece_${i}_Surface_0.xyz \
                 ../../SfSpp_preprocessing/Temp/Data/Pot_A/Pot_A_Piece_${i}_Surface_1.xyz \
                 ../../SfSpp_preprocessing/Dataset/Axes/Pot_A_Piece_${i}_Axis.xyz
done

# 3. breaklines
cd ../../SfSpp_preprocessing/build && ./EdgeLineExtraction

# 4. reassembly
cd ../../structure-from-sherds-pp/build
SFS_DATA_ROOT=../../SfSpp_preprocessing/Dataset/ SFS_AUTOSAVE=exit ./SfSpp
```

Note the ordering: upstream's README runs axis extraction last, but
`EdgeLineExtraction` reads the axis while classifying rim segments, so it has to
come second.
