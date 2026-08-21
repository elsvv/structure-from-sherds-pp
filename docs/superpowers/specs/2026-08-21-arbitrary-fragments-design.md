# Reassembling arbitrary fragment collections — design

**Date:** 2026-08-21
**Status:** revised after independent review — see
`docs/superpowers/notes/2026-08-21-plan-review.md`
**Revision:** every number below was measured against the code and the data on
disk. Where a figure in the previous revision was wrong, the correction is
marked `(was: …)` so the earlier reasoning stays traceable.

## Goal

Take a folder of raw 3D-scanned pot fragments and reassemble them, without
hand-tuning parameters per vessel. The fragment count is not known in advance —
it must work for four fragments and for a hundred, with no rework in between.

## Metric

Success is the **Sherd Accuracy (SA)** defined in the SfS++ paper (arXiv
2502.13986, §V-B): the number of sherds with **at least one correct edge**,
where edge `(i, j)` is correct when the *relative* transform `T_i⁻¹ T_j` matches
the ground-truth relative transform within `τ_R = 20°` (excluding axis deviation)
and `τ_t = 50 mm`, evaluated over the ground-truth adjacency matrix
`Ground Truth/Pot_<X>_simple_graph.txt`. It is deliberately **not** anchored on
the base fragment — the paper's pipeline does not use the base for
reconstruction.

**The binary already computes this.** `class/data_structure.cpp:920 CountResult`
implements exactly that test (`rad_threshold = 0.35` rad = 20.05°,
`t_threshold = 50.0`) and `SaveAcc` writes `Result/1. Acc.txt`. Nothing needs to
be reimplemented; it only needs to be reachable without a keypress (today it sits
behind the interactive `f` key — see the plan's Task 3).

Two consequences worth stating plainly:

- The Python Kabsch scorer built during diagnosis is a **diagnostic, not the
  gate**. It measures absolute vertex displacement after one global rigid
  alignment of the whole assembly, which is a different quantity: it applies no
  rotation test, its single global fit is dragged off by whichever fragments are
  misplaced, and SA's "at least one correct edge" is far more permissive than
  "every vertex within 20 mm". On Pot A both read 8/8, which is why the
  difference stayed invisible; it surfaces the moment a partly-correct pot
  (D, H, I, J) is scored. Keep the scorer — "how far off is this fragment, in
  millimetres" is genuinely useful and SA does not tell you — but never present
  its output as the published number.
- *(was: "mean placement error < 20 mm after one global rigid alignment" — that
  threshold and that metric appear in neither paper.)*

## Definition of done

Measured, not eyeballed. One command reproduces this table.

Measured against the **SfS++ dataset** (`Dataset/SfS_pp`, ten pots A–J), a
superset of the ICCV 2021 collection — Pot A is geometrically identical in both
(the two files differ in one line, `mtllib`; same 100 064 `v`, 200 124 `f`,
600 372 `vt`). Targets are the accuracies published on the SfS++ project page,
back-solved over the number of fragments the authors actually evaluate.

| Pot | evaluated | target SA | equals | now |
|---|---|---|---|---|
| A | 8 | 100 % | 8/8 | **0/8** |
| B | 9 | 100 % | 9/9 | not run |
| C | 4 | 100 % | 4/4 | not run |
| D | 28 | 82.1 % | 23/28 | not run |
| E | 31 | 100 % | 31/31 | not run |
| F | 6 | 100 % | 6/6 | not run |
| G | 7 | 100 % | 7/7 | not run |
| H | 11 | 90.9 % | 10/11 | not run |
| I | 27 | 92.6 % | 25/27 | not run |
| J | 11 | 72.7 % | 8/11 | not run |

**All ten are gateable.** *(was: D at 24/29, and F, G, I, J "reported, not
gated".)* The two reasons the previous revision gave were both wrong:

- `24/29` mixed sources: `24/28` is ICCV 2021 Table 2, while 82.1 % is SfS++.
  82.1 % of 29 has no integer solution; over 28 it is exactly 23/28.
- The apparent page-vs-code disagreements dissolve on inspection. **Pot G**: the
  project page's Pot G card is a verbatim copy of the Pot B card (a copy-paste
  bug); the page's own table image says 7 sherds, agreeing with `SHARD_NUMBER`.
  **Pot J**: `POT_J` declares `SHARD_NUMBER 19` but supplies only 12 array
  initialisers — a C++ bug, not a data disagreement (`POT_C_J` has the same
  defect, 23 vs 16). **Pots D, F, I**: the paper excludes very small fragments
  from evaluation, and the dataset says exactly which — they are the **all-zero
  rows** of `Pot_<X>_simple_graph.txt`, and precisely the fragments with no
  `_T.txt`. Verified on disk:

```
pot  graph NxN  edges  GT _T.txt  excluded piece(s)
A    8x8         15    8          -
B    9x9         15    9          -
C    4x4          5    4          -
D   29x29        69   28          22
E   31x31        62   31          -
F    7x7          9    6          7
G    7x7         10    7          -
H   11x11        17   11          -
I   30x30        65   27          28, 29, 30
J   12x12        21   11          9
```

The mechanism for removing them already exists — `shard_on_off` — and the
authors applied it only to Pot J. Setting it `false` for D 22, F 7 and I 28–30 is
a prerequisite for any run to match the published percentage, because
`CountResult`'s `total_sherd` counts enabled shards.

Gate on **SA, not Edge Accuracy**: Pot D's published EA of 83.3 % does not
back-solve over the 69 edges in the graph (it does over 66). Two further page
errors, recorded so nobody re-derives from them: the per-pot sherd counts sum to
143, not the 142 advertised; and Pot F's "# Edges: 7" disagrees with the 9 edges
in `Pot_F_simple_graph.txt`.

Two conditions on top of the numbers:

1. **No per-vessel tuning.** The same invocation produces the whole table. If
   Pot D needs its own parameters, the pipeline will not work unattended on new
   scans either.
2. **The chain starts at raw meshes**, not at the authors' preprocessed files.
   Feeding their files is what already works and is not the problem.

## What is already established

Verified during diagnosis and re-verified during review, so the design does not
have to re-litigate it:

| Component | State | Evidence |
|---|---|---|
| SfS++ reassembly | **sound** | 8/8, pooled mean 3.74 mm, on the authors' preprocessed Pot A |
| macOS port | **sound** | same run; both repositories build and run |
| `extract_axis` | **sound** | 0.35° median vs the authors' axes, on their surfaces |
| Breakline extraction | **sound, not merely close** | segment counts match the reference *exactly* (5,3,4,4,4,2,3,3); point counts within +3.2 % |
| Segment-header writing | **sound** | byte-compatible with the reference format |
| Base detection | **works, but inert** | ports correctly, and `NO_BASE_INFO` discards the result before reassembly |
| `surface_fr` | **never produced** | by us or by the authors' Pot A — an unlisted difference from the published configuration |
| Surface segmentation | **suspect** | axes derived from our surfaces: 3.5° vs 0.35° |
| Rim detection | **suspect** | 1 rim found where the reference marks 5 |

Three corrections to the previous revision:

- **Mean error is 3.74 mm pooled over all vertices** *(was: 4.4 mm)*. The 4.4
  figure is the unweighted mean of the eight per-fragment means (4.37). Both are
  defensible; the harness reports the pooled figure, and any gate must say which.
- **Rim detection finds one rim, piece 04** *(was: "0 rims found")*. The
  reference marks five. So the criterion is not dead, it is far too strict.
- **Breakline point counts differ by at most +3.2 %** *(was: +5–10 %)*. Measured
  deltas: +0.5, +2.4, +3.2, 0, 0, 0, 0, +1.7 %. This *strengthens* the case that
  breakline geometry is fine and only the rim flag is wrong — it is the best
  evidence available for prioritising the rim fix.

**On the base flag.** `main.cpp:32` has `#define NO_BASE_INFO` active — it is
present in the first upstream commit — and it zeroes `is_seg_base_` for every
shard immediately after loading (`main.cpp:105-109`). Base detection therefore
has no effect on reassembly whatsoever, which matches the paper's claim of being
base-agnostic. Independently: the SfS++ reference breaklines for Pot A carry
**no base flag at all** (`info` is never 2), while the older ICCV files mark
piece 1 as base. Nobody should spend time on base detection.

## Components to build

### 1. Arbitrary fragment count — the stated goal, and currently unimplemented

Nothing in the pipeline makes the fragment count arbitrary today.
`class/data_path.h` is 6 790 lines of 30 mutually exclusive `#ifdef` blocks; each
enumerates `SHARD_NUMBER` plus seven fully spelled string arrays by hand. There
is no globbing anywhere in the reassembly. A new collection of *n* fragments
means writing roughly 9*n* lines and recompiling. On the preprocessing side the
code *does* glob, but bakes the pot letter into a compile-time `#define POT_A`,
so one binary serves one pot.

Two changes remove both limits:

- **Generate `class/data_path.h` from the staged directory listing** — a small
  Python generator writing one preset from whatever fragments are present.
  `data_path.h` is included by `main.cpp` alone and `SHARD_NUMBER` appears
  nowhere else, so regenerating it recompiles one translation unit: measured at
  about 5 seconds on this machine, which makes "generate, rebuild, run" a
  perfectly ordinary step in the harness rather than a per-collection porting
  job.
- **`SFSPP_POT_ID` in the preprocessing `data_path.h`**, so `potID` is read at
  runtime with the compiled value as fallback and one preprocessing binary
  serves every collection.

This is listed first because it is the goal, and because it is also the cleanest
fix for the reference sweep, which cannot run without it.

### 2. Harness — everything else is measured through it

One command over a folder of raw meshes: convert, generate the preset, build,
preprocess, extract axes, extract breaklines, reassemble, report. For the
reference pots it additionally reads back `Result/1. Acc.txt` and prints the
table above.

Required outputs per run:

- Sherd Accuracy and Edge Accuracy from the binary, when ground truth exists
- per-fragment millimetre displacement from the Kabsch diagnostic, same condition
- the intermediate counts diagnosis needs: clusters per fragment, breakline
  points and segments, candidate axes per fragment, matches before and after
  pruning
- *(dropped: "a rendering of the result". The reassembly's viewer is
  interactive-only and keyboard-layout dependent — which is why `SFS_AUTOSAVE`
  exists — and no task was ever going to build a renderer. The saved result OBJs
  can be opened in any viewer; join geometry is judged numerically instead, by
  nearest-surface distance between placed fragments.)*

Crucially it must distinguish **three** outcomes, not two:

- assembled correctly
- assembled incorrectly (fragments placed, joins false)
- **no joins found** — the honest answer when the fragments were never
  neighbours, which no amount of code can change

### 3. Diagnosis — which defect is fatal

With the harness in place, run Pot A substituting one component at a time from
the authors' files. Four runs, not three:

| Run | Surfaces | Axes | Breaklines | Answers |
|---|---|---|---|---|
| A | ours | **theirs** | ours | is the axis error fatal |
| B1 | ours | ours | ours + **their rim flags** | is the missing rim flag fatal |
| B2 | ours | ours | **theirs, whole file** | is anything else about our breaklines fatal |
| C | ours → **theirs** | theirs (re-derived) | ours, or re-derived | is segmentation fatal |

Three corrections to the previous revision's table:

- **Run C's axes are not "ours".** Substituting surfaces requires re-running
  `extract_axis` on them, so Run C's axes are derived from theirs and Run C is
  confounded with Run A. Stated so nobody reads the result as a clean
  single-factor test.
- **Run B did not test the rim.** Symlinking the whole `_Breakline_0.pcd`
  replaces geometry, segment ranges *and* the rim flag together. Since segment
  counts already match exactly and point counts are within 3.2 %, the sharp
  experiment is B1: keep our breakline files and overwrite only the third header
  field with the reference values. That isolates the flag completely and costs
  one `sed`. B2 remains as the catch-all.
- **Surfaces feed the breaklines too.** `EdgeLineExtraction` consumes surfaces
  and the axis, so substituting surfaces without re-running it leaves breaklines
  derived from *our* surfaces — a fourth, unnamed mixed state. Either re-run
  `EdgeLineExtraction` in Run C and say so, or state explicitly that Run C tests
  only the registration-time use of `surface_in` / `surface_out`.

Fix what they point at — not both defects blind.

### 4. Length scale — an envelope, not a proportion

The previous revision proposed inferring a scale factor from wall thickness and
fragment extent, checked by "on the reference pots the inferred scales must come
out at 1.0 within a few percent". Measured across all ten pots, that design
fails its own check on nine of them.

Wall thickness (median nearest-surface distance, `Surface_0` → `Surface_1`) and
extent, per pot:

| pot | n | median thickness | min–max | median max-extent |
|---|---|---|---|---|
| A | 8 | **3.67** | 3.57–7.75 | 93.4 |
| B | 9 | 3.44 | 2.98–5.45 | 78.7 |
| C | 7 | 6.07 | 3.51–8.69 | 86.8 |
| D | 32 | **8.20** | 4.10–13.44 | 79.8 |
| E | 34 | 6.79 | 3.87–11.31 | 73.5 |
| F | 7 | 4.72 | 4.30–6.74 | **35.6** |
| G | 7 | **2.80** | 2.42–5.09 | 54.4 |
| H | 11 | 5.87 | 4.53–11.27 | 72.1 |
| I | 30 | 7.51 | 3.74–11.25 | 86.3 |
| J | 19 | 5.53 | 4.74–10.28 | 40.9 |

Per-pot median thickness spans **2.80 to 8.20 mm — a factor of 2.93**; median
extent spans **35.6 to 93.4 mm — a factor of 2.63**. Units are consistent across
all ten. And **all ten reach their published accuracy at
`SFS_LENGTH_SCALE = 1.0`**. So the thresholds are demonstrably *not*
proportional to wall thickness: they are fixed, and they work across a 2.9×
range. A thickness-proportional scale would return 1.0 on Pot A by construction
and 0.76–2.23 elsewhere, perturbing nine pots that currently work.

**Reframed goal:** infer a scale only when a collection falls *outside* the
envelope the fixed thresholds are known to cover — roughly **2.4–13.4 mm**
per-fragment thickness and **22–211 mm** per-fragment extent. Inside the
envelope, scale is 1.0 by definition. The correctness check becomes achievable
and meaningful — "scale = 1.0 for every one of the ten reference pots" — and it
is testable on all ten instead of tautologically on one.

Corrections to the constants and figures behind the old design:

- `REFERENCE_THICKNESS_MM = 3.69` is sound: measured 3.6731 (3.6736 at full
  resolution), seven of eight Pot A fragments in a tight 3.57–3.94 band.
- `REFERENCE_EXTENT_MM = 119.0`, commented "median bbox diagonal-max across the
  8 fragments", **is neither** *(was: presented as a median)*. Pot A's median
  max-extent is **93.42** and its median bbox diagonal is **126.50**; 119.0 is
  `Pot_A_Piece_01`'s max extent — the *largest* fragment, at the 89th percentile
  of all 164 reference fragments. Calibrating on it over-estimates every other
  pot by up to 2.63×.
- Drop the claim that "piece 1, the base fragment, is 7.26 mm". Not reproducible
  under any definition tried (S0→S1 median 7.81, S1→S0 7.06, symmetric 7.38) and
  piece 1's distance distribution is bimodal, so a single number for it is not
  well defined.
- Pot A is **atypically thin**: 8 of 10 pots are thicker and the dataset-wide
  median is 6.68 mm, 1.8× Pot A. It is the worst available choice of calibration
  pot — another reason to prefer an envelope over a ratio.
- **Our own fragments measure ~37.6 mm, a ratio of 10.2×** *(was: 8.25 mm, ratio
  2.2)*. Independently re-measured during this revision on
  `DatasetB/Surfaces/Pot_A/*_Surface_{0,1}.xyz`: 37.56, 37.21, 37.95, 38.96 —
  median 37.8. The segmentation is not obviously to blame: `Surface_0` and
  `Surface_1` have opposite mean normals (dot ≈ −1.0), exactly as the reference
  does, so they really are the two faces of the shell. The extent ratio is 5.4×
  (93.4 → 508.5), close to the "4–5" previously stated. So the two ratios
  disagree by 1.9×, not by the factor previously claimed — and see "The scale of
  our own scans" below, which may explain both.
- **`SFS_LENGTH_SCALE` is a partial knob and must be described as one.**
  `sfsLen()` / `sfsArea()` wrap exactly **12 literals** across 11 lines *(was:
  "roughly fifteen thresholds")*, while roughly **45** comparable millimetre
  constants in the same three files are untouched — `RejectOutlier` distance
  gates, `CauchyLoss` scales on distance/axis/rim residuals, `isConverge`
  translation terms, `t_threshold` in the ranking system, the LCS feature
  quantisation steps in `FeatureCompGraphBuilding`. On a collection whose scale
  genuinely differs, moving 12 of 57 leaves the pipeline internally inconsistent
  — arguably worse than moving none. Either scale all of them or scope the knob
  honestly.
- **Do not infer thickness from `Surface_0`/`Surface_1`.** That is the output of
  surface segmentation, which this same design lists as a prime suspect for the
  total failure. Calibrating the thresholds that segmentation's output is judged
  by, on a measurement taken from that output, is circular. Derive thickness
  from the mesh (shell-to-shell distance across the triangulation, needing no
  segmentation), or defer the whole component until the diagnosis has cleared
  segmentation. The plan does the latter.

### 5. Fix preprocessing

Driven by the diagnosis. Two candidates, in the order they are likely to matter:

**Rim detection.** The shipped C++ heuristic finds one rim on Pot A where the
authors' own files mark five. `check_base_and_rim.m` carries a different,
stricter criterion — a segment is a rim when radius and height stay
near-constant along it (`std < 1 mm` on both, plus bounded first differences and
≥20 points). The base half of that same file was ported during diagnosis and
matches to the last threshold, so the approach is known to work; this is the
unported half. Small.

Rim pruning is what removes roughly a third of the candidate matches, all of
them false — which is why our run starts from 224 matches where the reference
starts from 76.

**It does nothing for our own fragments.** They have no rim to detect. The rim
fix is for the reference pots and for the diagnosis; it is not on the path to
`potB`. See "Fragments without rim or base" below.

**Surface segmentation.** Larger and less predictable. The two largest clusters
cover only 46–81 % of points even on the authors' own data.

### 6. Fallback: segment the mesh, not the point cloud

If diagnosis blames segmentation and the authors' code cannot be repaired
cheaply, there is a better approach available to us than the original.

The authors discard the mesh and run region growing over a **point cloud** using
PCA normals from 5 neighbours spanning about a millimetre — which is where the
noise comes from (measured: 46–68 % of neighbouring normal pairs exceed the 4.5°
growing threshold). We have the triangulated mesh. Segmenting it by dihedral
angle across shared edges is far more robust: face normals are exact and
connectivity is known.

A sherd is a thin shell — two smooth faces separated by a band of fracture
surface meeting them at sharp angles. A dihedral-angle cut separates them almost
trivially. Roughly two hundred lines.

Not the starting point. But if it is needed, the route is short.

## Scope and assumptions

### Single-vessel input is an assumption, and the mixed case is the general one

A folder of raw scanned sherds is not guaranteed to come from one vessel. The
dataset ships `Transformation/ABC`, `ABCDE`, `ABCDE(10,5)`, `DE`, and
`data_path.h` ships `POT_A_B_C`, `POT_D_E`, `POT_E_I`, `POT_C_J`,
`POT_A_B_F_G_H`, `POT_A_B_C_D_E` and `POT_All`. The project page publishes five
mixed experiments (Dishes 97.6 %, C+J 80.0 %, E+I 96.6 %, ABCDE 92.5 %, All
87.3 %) — the harder half of the published result, and the case that matches
"a box of sherds" most closely. `NUM_MIXED_SHERD` already exists to express it.

This design **assumes single-vessel input** and does not gate on mixed
collections, for one reason: cost. Paper Table IX puts the All-mixed experiment
at 42.7 h for (k,b) = (10,5) and 137 h for (20,10), against single-pot runs at
(5,3). Gating on that needs a budget decision first. The obligation this
imposes on the harness is narrower and must be met: it must not report "no
joins" when the real answer is "these are two pots".

### Fragments without rim or base — our actual case

`main.cpp:31-32` carries both switches:

```cpp
//#define NO_RIM_INFO
#define NO_BASE_INFO
```

`NO_BASE_INFO` is active, so base flags never reach reassembly (above).
`NO_RIM_INFO` is commented out, and it is the exact switch for "fragments with
no rim" — which is our four scans. The honest test for our data is therefore:
**does the pipeline reach its targets on a reference pot with `NO_RIM_INFO`
defined?** Paper Table X's "Non-enhanced" ablation suggests it degrades but does
not collapse. That run costs one rebuild and three small pots (A, B, C) and it
tells us, before spending anything on `potB`, whether a rimless four-fragment
collection is in scope at all.

### The scale of our own scans is unverified, and everything depends on it

Measured bounding boxes:

```
set              file                   verts     max extent   diagonal
test_fragments   FY234007_reduced.ply    615160     546.31       687.86
test_fragments   FY234021_reduced.ply    526278     470.28       635.07
test_fragments   FY234094_reduced.ply    670739     679.17       790.49
test_fragments   FY234104_reduced.ply    354074     370.91       486.43
potB             FY234007.ply           7689470     546.59       688.27
potB             FY234021.ply           6578460     470.51       624.95
potB             FY234094.ply           8384209     679.60       785.44
potB             FY234104.ply           4425897     371.12       481.16
```

Reference fragments: median max extent 74.7, dataset maximum 210.6. Ours are
371–680 in the same nominal units, with a 37.6-unit wall. The PLY headers say
`Created in RealityCapture` (three files) and `by Geomagic Studio` (one) —
photogrammetry output, which is **unscaled unless a scale constraint was
supplied at reconstruction**. The coordinates may not be millimetres at all.

**One physical measurement of one sherd with a caliper resolves this**, and until
it exists a "no joins found" verdict on `potB` is uninterpretable: it could
equally mean "the correspondence window is 5 units on a 500-unit object". This
is a precondition on the `potB` run, not a nicety.

Separately: `test_fragments` and `potB` are the same four scans decimated ~12×,
but their z-extents disagree (145.15 vs 89.76 for FY234021; 131.39 vs 91.64 for
FY234094; 103.41 vs 71.47 for FY234104). The decimated copies carry stray points
30–55 units off-surface in z. They must be cleaned before `test_fragments` is
used for any geometric statistic.

## Risk that no code removes

The four fragments in `data/potB` may simply never have been neighbours. Four
sherds from a vessel over half a metre across are unlikely to adjoin. This is a
property of the material, not of the software, which is why the harness reports
"no joins found" as a distinct outcome. If the reference pots reach their targets
and these four still do not join, that result means what it says.

## Sequencing

1. Repository hygiene, then the generated `data_path.h` + `SFSPP_POT_ID` — the
   arbitrary-count component, and the prerequisite for everything measured
2. Wire `CountResult` into the non-interactive path; fix `shard_on_off`,
   `POT_J`/`POT_C_J` shard counts and `POT_CTRL`'s ground-truth paths
3. Harness: staging, the Kabsch diagnostic, the runner and the report
4. Diagnosis — the four substitution runs
5. Fix whatever they identify (rim criterion, or mesh segmentation)
6. The rimless run on A, B, C — decides whether `potB` is in scope
7. Full reference sweep A–J against the corrected targets
8. Length-scale envelope, verified to be a no-op on all ten pots
9. Measure the physical scale of the `potB` scans, then run them
