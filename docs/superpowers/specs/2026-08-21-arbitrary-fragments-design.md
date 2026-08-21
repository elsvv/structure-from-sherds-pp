# Reassembling arbitrary fragment collections — design

**Date:** 2026-08-21
**Status:** approved, ready for planning

## Goal

Take a folder of raw 3D-scanned pot fragments and reassemble them, without
hand-tuning parameters per vessel. The fragment count is not known in advance —
it must work for four fragments and for a hundred, with no rework in between.

## Definition of done

Measured, not eyeballed. One command reproduces this table:

Measured against the **SfS++ dataset** (`Dataset/SfS_pp`, ten pots A–J), which is
a strict superset of the ICCV 2021 collection — Pot A is byte-identical in both.
Targets are the accuracies published on the SfS++ project page.

| Pot | Sherds | Target | Now |
|---|---|---|---|
| A | 8 | 8/8 | **0/8** |
| B | 9 | 9/9 | not run |
| C | 4 | 4/4 | not run |
| E | 31 | 31/31 | not run |
| H | 11 | 10/11 | not run |
| D | 29 | 24/29 | not run |
| F, G, I, J | 7, 7, 30, 19 | reported, not gated | not run |

F, G, I and J are not gated because the project page and `class/data_path.h`
disagree on their sherd counts (7 vs 6, 7 vs 9, 30 vs 27, 19 vs 11), so there is
no unambiguous target to hold them to.

Two conditions on top of the numbers:

1. **No per-vessel tuning.** The same invocation produces the whole table. If
   Pot D needs its own parameters, the pipeline will not work unattended on new
   scans either.
2. **The chain starts at raw meshes**, not at the authors' preprocessed files.
   Feeding their files is what already works and is not the problem.

## What is already established

Verified during diagnosis, so the design does not have to re-litigate it:

| Component | State | Evidence |
|---|---|---|
| SfS++ reassembly | **sound** | 8/8, mean 4.4 mm, on the authors' preprocessed Pot A |
| macOS port | **sound** | same run; both repositories build and run |
| `extract_axis` | **sound** | 0.35° median vs the authors' axes, on their surfaces |
| Base detection | **sound** | Pot A piece 1 → `info 2`, matching the reference header |
| Breakline extraction | **close** | segment counts match the reference almost exactly; point counts +5–10% |
| Segment-header writing | **sound** | byte-compatible with the reference format |
| Surface segmentation | **suspect** | axes derived from our surfaces: 3.5° vs 0.35° |
| Rim detection | **suspect** | 0 rims found where the reference marks 5 of 6 |

The failure is therefore inside preprocessing, and within preprocessing it is
one or both of the two suspects. **Which one is fatal has not been established** —
that is the first thing the work must determine, before any fixing.

## Components to build

### 1. Harness — first, because everything else is measured through it

One command over a folder of raw meshes: convert, preprocess, extract axes,
extract breaklines, reassemble, report. For the reference pots it additionally
scores against ground truth and prints the table above.

It is built first not for tidiness but because **the diagnosis below is its
first use**, and there are dozens of such runs ahead. Much of it already exists
as loose scripts written during diagnosis — mesh/PLY→dataset converters, the
Kabsch scoring against ground truth, the renderers — and needs consolidating.

Required outputs per run:

- placed / total, and per-fragment error when ground truth exists
- the intermediate counts that diagnosis needs: clusters per fragment, breakline
  points and segments, axes per fragment, matches before and after pruning
- a rendering of the result

Crucially it must distinguish **three** outcomes, not two:

- assembled correctly
- assembled incorrectly (fragments placed, joins false)
- **no joins found** — the honest answer when the fragments were never
  neighbours, which no amount of code can change

### 2. Diagnosis — which defect is fatal

With the harness in place, run Pot A substituting one component at a time from
the authors' files:

| Run | Surfaces | Axes | Breaklines | Answers |
|---|---|---|---|---|
| A | ours | **theirs** | ours | is the axis error fatal |
| B | ours | ours | **theirs** | is the missing rim fatal |
| C | **theirs** | ours | ours | is segmentation fatal directly |

Three runs. Fix what they point at — not both defects blind.

Note that surfaces feed the pipeline twice: once as the input to axis
estimation, and again during registration (`surface_in`/`surface_out` are used
to refine the axis while sherds are added). Run C therefore tests something
runs A and B do not.

### 3. Scale inference

Roughly fifteen thresholds are absolute lengths in millimetres, calibrated for
~200 mm pots. `SFS_LENGTH_SCALE` and the `SFSPP_*` variables added during
diagnosis are a stopgap: nobody will hand-tune them per vessel.

They must be derived from the data. **Two reference quantities, not one**, because
the thresholds measure different things:

- **Wall thickness** — the median nearest-surface distance between the inner and
  outer surface of a fragment. Computed without any threshold of its own, and it
  is what the fracture band actually spans. Governs the thresholds that live
  across the fracture: the correspondence window (5 mm / 20 mm in
  `MakeCorWOTree`), the inlier threshold, the overlap area.
- **Fragment extent** — governs the thresholds that live along the profile:
  profile-curve bin size, edge-line resampling distance.

Measured: the authors' Pot A has a median wall thickness of 3.69 mm (piece 1,
the base fragment, is 7.26 mm — take the median across the collection, not per
fragment). Our fragments are 8.25 mm, a ratio of 2.2, while their extents differ
by a factor of 4–5. The two ratios genuinely disagree, which is why one global
knob is wrong.

**Correctness check:** on the reference pots the inferred scales must come out
at 1.0 within a few percent, or the change breaks what already works.

### 4. Fix preprocessing

Driven by the diagnosis. Two candidates, in the order they are likely to matter:

**Rim detection.** The shipped C++ heuristic finds nothing. `check_base_and_rim.m`
carries a different, stricter criterion — a segment is a rim when radius and
height stay near-constant along it (`std < 1 mm` on both, plus bounded first
differences and ≥20 points). The base half of that same file was ported during
diagnosis and produced the correct answer, so the approach is known to work;
this is the unported half. Small.

Rim pruning is what removes roughly a third of the candidate matches, all of
them false — which is why our run starts from 224 matches where the reference
starts from 76.

**Surface segmentation.** Larger and less predictable. The two largest clusters
cover only 46–81% of points even on the authors' own data.

### 5. Fallback: segment the mesh, not the point cloud

If diagnosis blames segmentation and the authors' code cannot be repaired
cheaply, there is a better approach available to us than the original.

The authors discard the mesh and run region growing over a **point cloud** using
PCA normals from 5 neighbours spanning about a millimetre — which is where the
noise comes from (measured: 46–68% of neighbouring normal pairs exceed the 4.5°
growing threshold). We have the triangulated mesh. Segmenting it by dihedral
angle across shared edges is far more robust: face normals are exact and
connectivity is known.

A sherd is a thin shell — two smooth faces separated by a band of fracture
surface meeting them at sharp angles. A dihedral-angle cut separates them almost
trivially. Roughly two hundred lines.

Not the starting point. But if it is needed, the route is short.

## Risk that no code removes

The four fragments in `data/potB` may simply never have been neighbours. Four
sherds from a vessel over half a metre across are unlikely to adjoin. This is a
property of the material, not of the software, which is why the harness reports
"no joins found" as a distinct outcome. If the reference pots reach their targets
and these four still do not join, that result means what it says.

## Sequencing

1. Harness (consolidate existing scripts, add the reference scoring table)
2. Diagnosis — the three substitution runs
3. Fix whatever they identify
4. Scale inference, verified to be a no-op on the reference pots
5. Full reference sweep A–E
6. Run `potB` with everything in place
