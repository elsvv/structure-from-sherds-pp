# Review of the arbitrary-fragments design and plan

**Reviewer:** independent
**Date:** 2026-08-21
**Documents reviewed:**
- `docs/superpowers/specs/2026-08-21-arbitrary-fragments-design.md`
- `docs/superpowers/plans/2026-08-21-arbitrary-fragments.md`

Everything below was checked against the code and the data on disk. Where I could
not check something I say so.

---

## 1. Verdict

**Not fit to execute as written.** The diagnostic instinct is right and the two
established results survive scrutiny — I reproduced 8/8 on the authors' preprocessed
Pot A and 0/8 on ours, from the artefacts already on disk. But the plan measures the
wrong quantity, against a target table that is partly copied from the wrong paper,
using a harness whose reference sweep cannot run at all. Specifically: the published
accuracy is a *relative-pose, per-edge* metric that the reassembly binary **already
implements** (`class/data_structure.cpp:920 CountResult`, 20° / 50 mm) and that the
plan reimplements in Python as something else entirely (absolute vertex error after a
global Kabsch fit, 20 mm, no rotation test); the Pot D target `24/29` is the ICCV 2021
numerator over the SfS++ denominator and is unreachable by construction; F, G, I and J
are left ungated on the strength of a copy-paste bug on the project page, when the
dataset itself determines their targets unambiguously; Task 8 fails for four
independent reasons before it produces a number; and the face-expansion rule that the
plan states as universal is in fact conditional and is *false* for every mesh the
harness itself stages. Separately, the scale-inference design contradicts its own
correctness criterion — the ten reference pots differ in wall thickness by 2.9×, and
all ten already meet their published accuracy at `SFS_LENGTH_SCALE=1.0`, so inferring
a scale from thickness would perturb nine of them. Tasks 1, 2 (minus its gate), 4 and
5 are close to sound and can proceed after the corrections below; Tasks 3, 7 and 8
need redesign, not repair.

---

## 2. Defects

Ordered by severity. **Wrong** = I verified it is incorrect. **Unverified** = the plan
asserts it and nobody has checked; I say which way I lean.

---

### D1 — WRONG (blocking). The harness measures a different quantity from the published targets, and the correct metric is already in the binary.

The plan's Global Constraints say:

> Success threshold per fragment, from the paper: mean placement error < 20 mm after
> one global rigid alignment of the whole assembly.

That threshold is in neither paper, and the metric is in neither paper.

**What the published numbers actually are.** The SfS++ paper (arXiv 2502.13986, §V-B)
defines Sherd Accuracy as *"the number of correctly reconstructed sherds with at least
one correct edge. We set threshold values of τ_R = 20° for 3D rotation (excluding axis
deviation) and τ_t = 50 mm for translation."* An edge (i,j) is correct when the
**relative** transform `T_i⁻¹ T_j` matches the ground-truth relative transform within
those thresholds, over the GT adjacency matrix `M`. It is explicitly not base-anchored:
*"this pipeline does not rely on the base fragment for reconstruction, which makes it
challenging to assess success based on specific sherds (previously, we used base
sherds as a reference)."*

**The binary already computes it.** `class/data_structure.cpp:920`:

```cpp
tuple<int, int, int, int> CountResult(MatrixXd& GT_graph, vector<Trans>& GT_trans,
                                      MatrixXd& graph, vector<Trans>& T_result,
                                      vector<bool>& right_sherd)
{
    double rad_threshold = 0.35, t_threshold = 50.0;   // 0.35 rad = 20.05 deg
    ...
    GT_T = GT_T_a * GT_T_b;      // relative GT pose
    T    = T_a * T_b;            // relative result pose
    if (isSimilarTrans(GT_T, T, rad_threshold, t_threshold)) { k++; ... }
    ...
    cout << "********* Sherd accuracy : " << k_sherd << " / " << total_sherd << ...
    cout << "********* Edge accuracy  : " << k / 2   << " / " << total / 2   << ...
}
```

and writes `Result/1. Acc.txt` via `SaveAcc`. It consumes
`Ground Truth/Pot_<X>_simple_graph.txt` (loaded at `main.cpp:111-143`) and the
`_T.txt` transforms.

**Why the plan's metric is not a substitute.** `score_result` pools every fragment's
vertices, fits one Kabsch transform over the whole assembly, and thresholds each
fragment's mean vertex error at 20 mm. Three consequences:

1. No rotation test. A fragment rotated about its own centroid but centred correctly
   is scored by displacement only.
2. The single global fit is dragged by whichever fragments are misplaced. On a
   part-correct assembly (D, H, I, J — the only pots where the target is not 100%)
   the alignment is biased by the failures, so correct fragments can be pushed over
   20 mm. The published SA cannot go down this way; it is per-edge and local.
3. SA counts a sherd correct with **at least one** correct edge. That is far more
   permissive than "every vertex within 20 mm of ground truth". The two numbers are
   not comparable in either direction.

For Pot A everything is correct so both metrics say 8/8 — which is why this has not
surfaced yet. It will surface the moment the sweep reaches D, H, I or J.

**One real obstacle:** `CountResult` is reached only from the `vis.fine_` branch
(`main.cpp:441`), i.e. the interactive `f` key. `SFS_AUTOSAVE=exit` goes
`vis.first_` → `vis.save_` → `SaveResult` → `break` (`main.cpp:255-258, 355-364`) and
never scores. The `fine_` branch at `main.cpp:388-445` does not re-optimise; it
reconstructs `T_result` from the stored per-state graph transforms and then scores.

**Correction.** Extend the autosave path to run the `fine_` reconstruction and
`CountResult` + `SaveAcc` before exiting (e.g. `SFS_AUTOSAVE=score` or unconditionally
under autosave), and have the harness parse `Result/1. Acc.txt` or the
`********* Sherd accuracy : k / n` line. Keep the Python Kabsch scorer, but demote it
to a *diagnostic* — it is genuinely useful for "how far off is this fragment in mm",
which SA does not tell you — and stop presenting its output as the published metric.

---

### D2 — WRONG (blocking). The target table is wrong for D and needlessly ungated for F, G, I and J. The dataset determines all ten unambiguously.

The plan gates on `{"A": 8, "B": 9, "C": 4, "D": 24, "E": 31, "H": 10}` and drops
F, G, I, J because *"the project page and `class/data_path.h` disagree on their sherd
counts (7 vs 6, 7 vs 9, 30 vs 27, 19 vs 11)"*.

**Pot D `24/29` is from the wrong paper.** ICCV 2021 Table 2 reports Pot D with **28**
sherds and `24/28` at (b=10, k=20). SfS++ reports Pot D as 82.1 % SA. The plan has
taken the ICCV numerator and the SfS++ denominator. 82.1 % of 29 = 23.8 — no integer
solution. Over 28 it is exactly 23/28 = 82.14 %.

**Pot G does not disagree.** The `7 vs 9` comes from the HTML card on
`https://sj-yoo.info/sfs/`, which for Pot G is a verbatim copy of the Pot B card
(`Type: Dish, # Sherds: 9, # Edges: 15, Acc. : 100 %`). The page's own table image
(paper Table III) says Pot G = 7 sherds, 10 edges — agreeing with `SHARD_NUMBER 7`.

**F, I and J do not disagree either — the code just forgot to disable the excluded
fragments.** The paper excludes very small fragments from evaluation. The dataset
encodes exactly which, in two redundant ways, and I checked both:

```
pot  graph NxN  edges  GT _T.txt files  missing piece(s)  isolated rows in graph
A    8x8         15    8                []                []
B    9x9         15    9                []                []
C    4x4          5    4                []                []
D   29x29        69   28                [22]              [22]
E   31x31        62   31                []                []
F    7x7          9    6                [7]               [7]
G    7x7         10    7                []                []
H   11x11        17   11                []                []
I   30x30        65   27                [28, 29, 30]      [28, 29, 30]
J   12x12        21   11                [9]               [9]
```

The excluded fragments are precisely the all-zero rows of
`Ground Truth/Pot_<X>_simple_graph.txt`, and precisely the ones with no `_T.txt`.
The mechanism the authors use to remove them is `shard_on_off` — and **they applied it
only to Pot J**: `POT_J` has `shard_on_off[8] = false`
(`Pot_J_Piece_09_Breakline_0.pcd`), `POT_C_J` has index 12 false, `POT_All` has index
144 false — all three the same J piece 9. `POT_D`, `POT_F` and `POT_I` are all-true.

Back-solving the page's percentages over the connected sub-count gives clean integers
for every pot:

```
A: 100%  x 8  = 8      F: 100%  x 6  = 6
B: 100%  x 9  = 9      G: 100%  x 7  = 7
C: 100%  x 4  = 4      H: 90.9% x 11 = 9.999  -> 10
D: 82.1% x 28 = 22.988 -> 23     I: 92.6% x 27 = 25.002 -> 25
E: 100%  x 31 = 31     J: 72.7% x 11 = 7.997  -> 8
```

**Corrected target table (Sherd Accuracy, SfS++, k=5 b=3):**

| Pot | evaluated | target SA | equals |
|---|---|---|---|
| A | 8 | 100 % | 8/8 |
| B | 9 | 100 % | 9/9 |
| C | 4 | 100 % | 4/4 |
| D | 28 | 82.1 % | 23/28 |
| E | 31 | 100 % | 31/31 |
| F | 6 | 100 % | 6/6 |
| G | 7 | 100 % | 7/7 |
| H | 11 | 90.9 % | 10/11 |
| I | 27 | 92.6 % | 25/27 |
| J | 11 | 72.7 % | 8/11 |

All ten are gateable. Nothing needs to be "reported but not gated".

**Prerequisite:** set `shard_on_off = false` for D piece 22, F piece 7, I pieces
28–30 in `class/data_path.h`, otherwise `CountResult`'s `total_sherd` (which counts
`shard_on_off` = true, `data_structure.cpp:933-936`) uses 29 / 7 / 30 and no run can
ever match the published percentage. `POT_All` needs the same five.

Two page errors worth recording so nobody re-derives from them: the page's per-pot
sherd counts sum to 143, not the 142 it advertises (Pot D printed as 29 where the
evaluation uses 28), and Pot F's "# Edges: 7" disagrees with the 9 edges in
`Pot_F_simple_graph.txt`. Pot D's EA of 83.3 % also does not back-solve over 69 edges
(it does over 66) — so gate on SA, not EA.

---

### D3 — WRONG (blocking). Task 8 cannot run. Four independent reasons.

**(a) The source glob matches nothing for pots C–J.** `_pot_source` globs
`Pot_{pot}_Piece_*_Mesh.obj`. Only Pots A and B ship full-resolution meshes; the other
147 of 164 files are `_Mesh_DS.obj`:

```
$ ls Dataset/SfS_pp/Mesh | grep -c _DS      -> 147
$ ls Dataset/SfS_pp/Mesh | grep -c '\.obj$' -> 164
Pot_C_Piece_01_Mesh_DS.obj      # does not match Pot_C_Piece_*_Mesh.obj
```

`class/data_path.h` agrees: `POT_A`/`POT_B` use `_Mesh.obj`, `POT_C`…`POT_J` use
`_Mesh_DS.obj`. `score_result` has the same bug — it builds
`Pot_{pot_id}_Piece_{i:02d}_Mesh.obj`, absent for C–J.

**(b) The preprocessing binaries are compiled for one pot.** `SfSpp_preprocessing/data_path.h:5`
is `#define POT_A`, and `potID` is a compile-time constant feeding
`getPointDatasetPath(potID)` etc. `run_pipeline` sets `SFSPP_DATASET_ROOT` and calls
`./MeshPreprocessing`, which reads `<root>/Point/Pot_A/`. `stage_dataset` writes
`<root>/Point/Pot_B/` for pot B. The sweep would preprocess nothing for nine of ten
pots. `SFSPP_*` gives you a configurable *root*, not a configurable pot.

**(c) Staged layout and per-pot presets disagree on directory nesting.**
Preprocessing writes **nested**: `Breaklines/Pot_A/`, `Surfaces/Pot_A/`,
`Mesh/Pot_A/` (`data_path.h:57-70`, confirmed on disk in `DatasetA/`). The
`POT_A`…`POT_J` presets read **flat**: `path + "Breaklines/Pot_A_Piece_01_Breakline_0.pcd"`,
`path + "Mesh/Pot_A_Piece_01_Mesh.obj"`. Only `POT_CTRL` and `POT_TEST` use the nested
form. So `SfSpp_A` … `SfSpp_J`, built exactly as Task 8 specifies, will not find a
single file the harness produced. (`Axes/` is flat everywhere — that one is consistent.)

**(d) `--pot-id` is largely cosmetic.** `POT_TEST` (the default target `SfSpp`) and
`POT_CTRL` are both hard-wired to `Pot_A_Piece_NN` names. Any `--pot-id` other than `A`
stages files the binary never reads. Task 9 quietly relies on this by passing
`--pot-id A` for `potB`; Task 8 does not.

**Correction.** Either (i) resolve the mesh suffix per file and normalise every
collection to a single staged layout, add `SFSPP_POT_ID` as an environment override in
the preprocessing `data_path.h` so one binary serves all pots, and add per-pot presets
that read the nested layout; or (ii) drop per-pot presets entirely and generate
`class/data_path.h` from a manifest (see Gap G1, which you need anyway).

---

### D4 — WRONG (blocking). Face expansion is conditional, and false for every mesh the harness stages.

The plan states this as an unconditional property:

> Result OBJs written by the reassembly are face-expanded: vertex `3k+j` is corner `j`
> of face `k` of the source mesh. Any comparison against a source mesh must expand the
> source the same way (`V[F.reshape(-1)]`), or correspondence is lost.

It is true only when the *input* OBJ carries per-corner attributes. `SaveResult` →
`Visualize::SaveMesh` → `pcl::io::saveOBJFile` on a mesh loaded by
`pcl::io::loadPolygonFileOBJ`; VTK duplicates points on load when the OBJ has per-corner
`vt` indices. The authors' meshes have them (`Pot_A_Piece_01_Mesh.obj`: 100 064 `v`,
200 124 `f`, **600 372 `vt`**). Meshes written by the plan's `write_obj` have none.

Measured, on the run already on disk:

```
staged source  DatasetB/Mesh/Pot_A/Pot_A_Piece_01_Mesh.obj : v=615160  f=1230314  vt=0
result         DatasetB/Result/..._Top_1_OBJ_1.obj          : v=615160  f=1230314
                                                              (not 3 x 1230314)

authors' mesh  Dataset/SfS_pp/Mesh/Pot_A_Piece_01_Mesh.obj  : v=100064  f=200124  vt=600372
result         DatasetRef/Result/..._Top_1_OBJ_1.obj         : v=600372  (= 3 x 200124)
```

So `score_result`'s unconditional `expand_faces` is correct for reference runs and
wrong for every run staged by `stage_dataset` — which is Task 9 (the user's own
fragments) and every pot in Task 8. The consequence is the exact silent-correspondence-
loss failure the plan warns about, with the sign flipped.

**Correction.** Pick the layout from the data, not from a rule:

```python
Vr, _ = read_mesh(path)
Vo, Fo = read_mesh(src)
if len(Vr) == 3 * len(Fo):
    src_pts = expand_faces(Vo, Fo)
elif len(Vr) == len(Vo):
    src_pts = Vo
else:
    raise ValueError(f"{path}: {len(Vr)} verts matches neither {len(Vo)} nor {3*len(Fo)}")
```

and delete the `n = min(len(src), len(res))` truncation, which currently hides exactly
this mismatch instead of raising on it.

---

### D5 — WRONG (blocking). The Task 2 gate expects a number the specified code does not produce, and instructs the worker to stop.

Task 2 Step 5:

> Expected output: `8 / 8 mean 4.4` (± 0.1)
> If this does not reproduce, stop — the scorer is wrong and every later measurement
> would be worthless.

I implemented `score_result` verbatim from the plan and ran it against the artefacts
already on disk:

```
DatasetRef  8 / 8  placed 8  mean 3.739
   per-fragment: {1: 1.52, 2: 4.85, 3: 5.08, 4: 4.56, 5: 4.09, 6: 5.22, 7: 3.84, 8: 5.79}
DatasetA    0 / 8  placed 8  mean 89.463
   per-fragment: {1: 74.51, 2: 68.86, 3: 142.40, 4: 96.11, 5: 87.29, 6: 85.21, 7: 75.26, 8: 146.66}
```

The 8/8 and 0/8 both reproduce (see Confirmations C1). The mean does not: **3.739, not
4.4**. The discrepancy is a definition, not an error — `mean_error` pools all vertices,
so the large, well-placed piece 1 dominates; the unweighted mean of the eight
per-fragment means is `34.95 / 8 = 4.369 ≈ 4.4`, which is where the plan's figure comes
from.

As written, a correct implementation trips the stop condition on its first use.

**Correction.** Expect `8 / 8 mean 3.74 (± 0.05)`, or change `mean_error` to
`float(np.mean(list(out["per_fragment"].values())))` and keep 4.37. State which
definition is meant. The design's "mean error ~4.4 mm" needs the same fix.

---

### D6 — WRONG (design-level). Scale inference contradicts its own correctness criterion, and both reference constants are mismeasured.

The design's check is:

> **Correctness check:** on the reference pots the inferred scales must come out at 1.0
> within a few percent, or the change breaks what already works.

I measured wall thickness (median nearest-neighbour distance, `Surface_0` → `Surface_1`,
3000 sampled points per fragment, voxel-grid NN validated against exact brute force on
Pot A to 3.6e-15) and bounding-box extent for all 163 surface pairs and 164 meshes:

| pot | n | median thickness | min–max | median max-extent | vs Pot A thickness |
|---|---|---|---|---|---|
| A | 8 | **3.6731** | 3.566–7.749 | 93.42 | 1.00 |
| B | 9 | 3.4386 | 2.980–5.446 | 78.73 | 0.94 |
| C | 7 | 6.0673 | 3.506–8.688 | 86.79 | 1.65 |
| D | 32 | **8.1993** | 4.095–13.436 | 79.79 | **2.23** |
| E | 34 | 6.7863 | 3.874–11.314 | 73.51 | 1.85 |
| F | 7 | 4.7183 | 4.297–6.742 | **35.58** | 1.29 |
| G | 7 | **2.7997** | 2.416–5.091 | 54.38 | **0.76** |
| H | 11 | 5.8715 | 4.525–11.270 | 72.08 | 1.60 |
| I | 30 | 7.5089 | 3.738–11.245 | 86.27 | 2.04 |
| J | 19 | 5.5321 | 4.735–10.276 | 40.94 | 1.51 |

Per-pot median thickness spans **2.80 (G) to 8.20 (D) — a factor of 2.93**. Per-pot
median max-extent spans **35.58 (F) to 93.42 (A) — a factor of 2.63**. Units are
consistent across all ten pots (common frame, z ∈ 287–587 everywhere); no pot is in
different units.

All ten pots already reach their published accuracy with `SFS_LENGTH_SCALE = 1.0`.
Therefore an inferred thickness scale would come out at 1.0 on Pot A **by construction**
and at 0.76–2.23 on the other nine — violating the stated criterion on nine of ten
reference pots, and changing the thresholds on nine pots that currently work. The
design's own test, applied honestly, refutes the design's approach.

Task 7 Step 5 only ever runs the check on Pot A, and passes `[119.0]` as the extent
list, so `extent_scale` is 1.0 tautologically. The two `assert ... - 1.0 < 1e-9` unit
tests are likewise tautological — they feed the reference constants back in.

**Both constants are also mismeasured:**

- `REFERENCE_THICKNESS_MM = 3.69` — measured 3.6731 (3.6736 at full resolution).
  Close enough; the number is fine.
- `REFERENCE_EXTENT_MM = 119.0`, commented *"median bbox diagonal-max across the 8
  fragments"* — Pot A's median max-extent is **93.42** and its median bbox diagonal is
  **126.50**. 119.0 is neither: it is `Pot_A_Piece_01`'s max extent (119.164), i.e. the
  **largest** fragment, and it sits at the 89th percentile of all 164 reference
  fragments. Calibrating on it over-estimates every other pot, by 2.63× for Pot F.
- The design's *"piece 1, the base fragment, is 7.26 mm"* is not reproducible under any
  definition I tried: S0→S1 median 7.807, S1→S0 median 7.055, symmetric pooled 7.383.
  Piece 1's distance distribution is broad and bimodal (p25 6.12, p75 10.18), so a
  single number for it is not well defined. Don't cite 7.26 as measured.
- Pot A is *atypically thin*: 8 of 10 pots are thicker, and the dataset-wide median is
  6.68 mm, **1.8× Pot A**. It is the worst available choice of calibration pot.

**The user's own fragments do not measure 8.25 mm either.** The design states *"Our
fragments are 8.25 mm, a ratio of 2.2, while their extents differ by a factor of 4–5."*
Applying the design's own definition to the only surfaces that exist for those
fragments — `SfSpp_preprocessing/DatasetB/Surfaces/Pot_A/*_Surface_{0,1}.xyz`, identical
to `TempB/Data/Pot_A/` — with the same validated NN routine:

```
piece 1: 37.3142     piece 3: 37.8390
piece 2: 36.7888     piece 4: 39.1160     median 37.5766   ratio to Pot A: 10.23x
```

Not 8.25, and not a ratio of 2.2 — **37.6, a ratio of 10.2**. The segmentation is not
obviously to blame for this: `Surface_0` and `Surface_1` have opposite mean normals
(dot −0.998 and −1.000 for pieces 1 and 2), exactly as the reference does (−1.000,
−0.999), so they really are the two faces of the shell. Measured extent ratio is
93.42 → 508.5, i.e. **5.4×**, which is close to the design's "4–5".

So under the design's own two-scale scheme the real figures are thickness ×10.2 and
extent ×5.4 — a 1.9× disagreement between the two, on numbers roughly 4.6× and 1.2×
larger than the design cites. I could not reproduce 8.25 under any definition I tried.
Either the figure came from a different run or a different measurement; it needs
re-deriving before anything is built on it.

**And there is a circularity.** The design proposes to infer the scale from
`Surface_0`/`Surface_1` — the output of surface segmentation, which the same design
lists as one of the two suspects for the total failure. A thickness measured from a
possibly-broken segmentation cannot calibrate the thresholds that segmentation's output
is then judged by. Derive thickness from the mesh instead (shell-to-shell distance
across the triangulation, which needs no segmentation), or defer Task 7 until after the
diagnosis clears segmentation.

**Correction.** The evidence says the thresholds are not thickness-proportional across
the range the authors tested — they are fixed and they work from 2.8 mm to 8.2 mm walls.
Reframe the goal: infer a scale only when the collection falls *outside* the range the
fixed thresholds are known to cover (roughly 2.4–13.4 mm thickness, 22–211 mm extent).
Then the correctness check becomes achievable and meaningful — "scale = 1.0 for any
collection inside the reference envelope" — and it is testable on all ten pots instead
of one. If a proportional scale is kept anyway, the gate must be *"the reference sweep
still reaches the same targets with inference enabled"*, run on all ten pots, not an
identity check on Pot A.

`_nearest_distances` also has a correctness bug worth fixing regardless: cell size is
`span/40` (≈3 mm for a 120 mm fragment) and the search covers only the 27 adjacent
cells, so a true nearest neighbour at 3.7 mm can be missed, and points whose 27-cell
neighbourhood is empty are silently dropped from the median. Expand the ring until
`best_so_far <= ring * cell`.

---

### D7 — WRONG (design-level). The ablation cannot isolate what the table claims it isolates.

The design's table:

| Run | Surfaces | Axes | Breaklines | Claimed answer |
|---|---|---|---|---|
| A | ours | theirs | ours | is the axis error fatal |
| B | ours | ours | theirs | is the missing rim fatal |
| C | theirs | ours | ours | is segmentation fatal directly |

Three problems.

**Run C's axes are not "ours".** Task 4 Step 4 re-runs `extract_axis` on the substituted
surfaces, which makes Run C's axes derived-from-theirs. Run C is therefore
`surfaces=theirs + axes=theirs-derived`, confounded with Run A. Not fatal to the
experiment (if C repairs and A does not, the registration use of the surfaces is
implicated) but the table is wrong and will mislead whoever reads the result.

**Run B does not test the rim.** It symlinks the authors' whole `_Breakline_0.pcd` —
geometry, segment ranges *and* the rim/base flags in the header. If B repairs the
assembly you learn "our breaklines are bad in some respect", not "the missing rim is
fatal". The design states it as the latter.

This one is cheap to fix properly, because the geometry is already essentially
identical. Comparing headers on disk:

```
piece   segs  points  info      segs  points  info
        (ours, DatasetA)        (theirs, SfS_pp)
01       5     209     2          5     208     2
02       3     173     0          3     169     1
03       4     161     0          4     156     1
04       4     151     1          4     151     1
05       4     147     0          4     147     0
06       2     150     0          2     150     1
07       3      65     0          3      65     0
08       3      61     0          3      60     1
```

Segment counts are **identical** for all eight; point counts differ by 0–3.2 %. The only
substantive difference is the `info` byte. So the sharpest Run B is: keep our breakline
files and overwrite only the third header field with the reference values. That
separates "rim flag" from "breakline geometry" completely, and takes one `sed`.

**Surfaces feed the breaklines too.** `EdgeLineExtraction` consumes surfaces and the
axis (`edgeline_extraction.cpp:2894-2897, 2991`). Substituting surfaces without
re-running `EdgeLineExtraction` leaves breaklines derived from *our* surfaces — a
fourth, unnamed mixed state. Either re-run EdgeLine in Run C and say so, or state
explicitly that Run C tests only the registration-time use of `surface_in`/`surface_out`.

**Two factual corrections to the design's evidence table while you are there:**

- *"Rim detection: 0 rims found"* — **false**. Our `DatasetA` finds one rim, piece 04
  (`info 1`). The reference marks five (pieces 2, 3, 4, 6, 8). So it is 1 of 5, not 0 of 5.
- *"Breakline extraction: point counts +5–10 %"* — **overstated**. Measured deltas are
  +0.5 %, +2.4 %, +3.2 %, 0 %, 0 %, 0 %, 0 %, +1.7 %. Maximum +3.2 %. This *strengthens*
  the case that breakline geometry is fine and only the flag is wrong — worth stating
  correctly, because it is the best evidence you have for prioritising Task 5.

---

### D8 — WRONG. Steps that fail as written.

**(a) `$I` is undefined.** Task 4 Step 2:

```bash
R="/Users/.../Dataset/SfS_pp"
...
substitute('/tmp/ablate-axes/Dataset', 'A', 'axes', '$I')
```

`$I` is never set in that block (only `$R` is). Inside the double-quoted `python3 -c`
string the shell expands it to the empty string, so `reference_root=''` and the
`FileNotFoundError(src)` guard fires on a relative path. Should be `'$R'`. This is the
same class of bug the commit `94ea633` claims to have fixed; one instance remains.

**(b) `pytest` is not installed and the given install command is blocked.**

```
$ /opt/homebrew/bin/python3 -m pytest --version
/opt/homebrew/opt/python@3.14/bin/python3.14: No module named pytest

$ /opt/homebrew/bin/python3 -m pip install --user --dry-run pytest
error: externally-managed-environment  (PEP 668)
```

The plan's fallback — `python3 -m pip install --user pytest` — fails on this
Homebrew Python. Every "Run test to verify it fails/passes" step in Tasks 1, 2, 3 and 7
is blocked at Task 1 Step 2. Also note the interpreter is **3.14.7**, not the 3.13 the
Tech Stack line claims.

Options: `pip install --user --break-system-packages pytest`; a venv at
`tools/.venv`; or rewrite the tests for `unittest` (available, stdlib) and run
`python3 -m unittest discover`. Pick one and put it in Global Constraints, not in a
parenthesis inside Task 1.

**(c) `substitute` hard-codes eight fragments.** `names = [... for i in range(1, 9)]`.
Fine for Task 4 (Pot A only), but the docstring and signature promise generality and
the function takes `pot_id`. Derive the list from the destination directory instead, or
rename it to make the Pot-A-only scope explicit.

**(d) `Result/` is never cleared between runs.** `run_pipeline` creates it with
`exist_ok=True` and then counts `stats["placed"] = len([f for f in os.listdir(result_dir)
if "_Top_1_OBJ_" in f])`. Filenames are `<Y>_<M>_<D>_<H>_<Min>_Top_1_OBJ_<i>.obj`
(`class/visualize.cpp:520-526`), so a second run in a different minute *adds* files
rather than replacing them. `placed` inflates, and `score_result`'s glob returns two
files with the same `OBJ_<i>` suffix and its sort key cannot tell them apart — it would
score a mixture of two runs. Task 4's ablation steps do `rm -rf .../Result` manually,
which shows the author knows; `run_pipeline` must do it for every run.

**(e) `--src "$R/Mesh"` picks up all 164 meshes.** The plan already flags this in a note
and tells the worker to copy eight files to a scratch folder. Make it the instruction,
not the fallback — as written, Step 5 of Task 3 says to run the broken form first.

---

### D9 — WRONG. `run_pipeline`'s interface does not match its implementation.

The declared return is
`{"clusters": {name: int}, "breakline_points": {name: int}, "breakline_segments": {name: int}, ...}`
— dicts keyed by fragment name. `parse_stats` returns **lists**
(`[int(x) for x in re.findall(...)]`), unkeyed and in log order. `cli.py` then prints
them as lists. Either fix the interface text or zip them against `names`; as written a
reader cannot tell which fragment a count belongs to, which is the whole point of
collecting them.

Also `re.search(r"Score : (\d+)", log)` matches the **first** occurrence in the
concatenated log, and there are two producers: `class/visualize.cpp:358`
(`"Result : r/n, Score : ..."`) and `class/visualize.cpp:414` (`"Score : ..."`). Anchor
the regex to line start.

---

### D10 — WRONG. Task 5's verification vector is incorrect.

Task 5 Step 3: *"The reference is `2, 1, 1, 1, 1, 1, 0, 0` for pieces 1–8."*

Measured from `Dataset/SfS_pp/Breaklines/Pot_A_Piece_<NN>_Breakline_0.pcd`, second
header line, third field:

```
piece   1  2  3  4  5  6  7  8
info    2  1  1  1  0  1  0  1
```

The correct vector is `2, 1, 1, 1, 0, 1, 0, 1`. As written the check would reject a
correct port (pieces 5, 6, 7, 8 all disagree with the plan's expectation) or accept a
wrong one.

Two supporting facts, both consistent with the corrected vector: `info 2` is the base
(piece 1), `info 1` is a rim — five rims, matching ICCV 2021 Table 2's *"# rim sherds
(detected / actual): 5/6"* for Pot A and paper Table III's 1 base / 6 rim.

---

### D11 — WRONG (magnitude). `SFS_LENGTH_SCALE` covers 12 thresholds, not "roughly fifteen", and ~45 comparable ones are unscaled.

The design says *"Roughly fifteen thresholds are absolute lengths in millimetres"* and
implies `SFS_LENGTH_SCALE` is the stopgap covering them.

Counted: `sfsLen()` / `sfsArea()` (defined `class/data_structure.h:33-45`) wrap exactly
**12 literals** across 11 lines — 7 lengths, 5 areas — in `feature_matching.cpp` and
`ranking_system.cpp`.

Still unscaled in the same three files, ~45 sites, including:

- `reconstruction.cpp` `RejectOutlier(cor, dist_TH, ...)` first argument, an mm gate:
  lines 460, 524, 587, 1187, 1298 (20 mm); 660, 685 (10 mm)
- `reconstruction.cpp` `CauchyLoss(5.0)` on `loss_dist`: 1215, 1417, 1610, 1783, 1970;
  `CauchyLoss(2.0)` on `loss_axis` and `loss_rim` (both length-valued residuals):
  1217/1218, 1419/1420, 1612/1613, 1785/1786, 1972/1973
- `reconstruction.cpp` `isConverge` translation term: 1279, 1475, 1684 (2.0 mm);
  1861 (1.0); 2042 (0.5); `FillMatchedPoints` `dist < 3.0` at 1105
- `ranking_system.cpp:20` `t_threshold = 20.0`; `:256`, `:315` `isSimilarTrans(..., 20.0)`;
  `:622` `isConverge(..., 10)`; `:1410`, `:1439` `CountPCInlier(..., 7.0, 1.5)`;
  `:1732` `RefineAxis(..., 0.5, ...)`
- `feature_matching.cpp:579, 581, 588` (`cor.len > 4`, `= 4`, `> 2`); `:1211, 1344`
  (`dist > 10`); `:1235-1236, 1365-1366` (rim-height slack 10); `:1596`
  `isSimilar(..., 20)`
- `feature_matching.cpp:1518-1521` — LCS feature quantisation steps, **unconditional**
  in `FeatureCompGraphBuilding`, the function the incremental builder actually calls
  (`ranking_system.cpp:792, 830`)
- outside the three files: `class/filter.cpp:248` (`> 1.0`),
  `class/data_structure.h:271, 275` (`t_threshold = 10` default)

This matters for Task 7: setting `SFS_LENGTH_SCALE` moves 12 of ~57 length-valued
constants. On a collection whose scale genuinely differs, the other 45 stay put and the
pipeline ends up internally inconsistent — arguably worse than leaving all of them
alone. Say so, and either scale all of them or scope `SFS_LENGTH_SCALE` honestly as a
partial knob.

---

### D12 — WRONG (minor but will corrupt history). Git hygiene in the preprocessing repo.

```
$ git -C /Users/vaceslaveliseev/@dev/SfSpp_preprocessing status --porcelain
 M data_path.h
 M edgeline_extraction.cpp
?? DatasetA/   ?? DatasetB/   ?? DatasetRef/   ?? TempA/   ?? TempB/
```

Two problems. First, `data_path.h` and `edgeline_extraction.cpp` carry **uncommitted**
work — the `SFSPP_DATASET_ROOT` / `SFSPP_TEMP_ROOT` plumbing that the whole plan depends
on. Task 5 Step 5 does `git add edgeline_extraction.cpp && git commit -m "Port the
reference rim criterion..."`, which would silently fold that plumbing into the rim
commit, and `data_path.h` would still never be committed by any task.

Second, that repo's `.gitignore` covers `Dataset/` and `Temp/` but **not** `DatasetA/`,
`DatasetB/`, `DatasetRef/`, `TempA/`, `TempB/` — 1.5 GB of untracked, un-ignored
artefacts. The plan's `git add -A` warning is scoped to the pp repo only; it applies
here more urgently.

**Correction.** Add a Task 0: commit the existing preprocessing changes on their own,
and extend `SfSpp_preprocessing/.gitignore` with `Dataset*/` and `Temp*/`.

---

### D13 — WRONG (minor, but the comment is misleading). `POT_CTRL`'s ground-truth paths do not resolve.

`class/data_path.h:133-146`:

```cpp
// Real ground truth: this run is scored against the authors' own transformations.
string gt_T_path[SHARD_NUMBER] = { path + "GroundTruth/Transformation/Pot_A_Piece_1_T.txt", ... };
string gt_graph_path[1]        = { path + "GroundTruth/Pot_A_simple_graph.txt" };
```

Against `Dataset/SfS_pp` the folder is `Ground Truth/` (with a space) and there is no
`Transformation/` subdirectory, so none of the nine paths resolve. Against the ICCV
layout that `DatasetRef/GroundTruth` symlinks to, the `_T.txt` paths resolve but the
graph does not — the file is at `GroundTruth/Transformation/Pot_A_simple_graph.txt`,
not `GroundTruth/Pot_A_simple_graph.txt`:

```
$ ls "ICCV Data/GroundTruth/Pot_A_simple_graph.txt"
No such file or directory
```

`main.cpp:122-142` reads the graph with an `ifstream` and no existence check, so
`num_raw = 0`, `GT_graph` stays all-zero, and the accuracy report is silently disabled.
The comment claiming this run is scored against real ground truth is false today.
`POT_A`…`POT_J` use the correct `Ground Truth/` form; only `POT_CTRL` and `POT_TEST` are
wrong.

---

### D14 — WRONG (upstream bug the plan inherits). `POT_J` and `POT_C_J` declare more shards than they define.

```
block         SHARD_NUMBER  file_path entries  shard_on_off entries
POT_J                   19                 12                   12
POT_C_J                 23                 16                   16
```

Both arrays are declared `[SHARD_NUMBER]` and given fewer initialisers, so the tail is
value-initialised to empty strings and `false`. `GT_graph` is sized 19×19 while
`Pot_J_simple_graph.txt` is 12×12. Any run of `POT_J` loads seven fragments with empty
paths. This is exactly the *"19 vs 11"* the plan reads as a page/code disagreement — it
is a C++ bug, and the effective count is 12 defined minus 1 disabled = **11**, matching
the page. Correcting `SHARD_NUMBER` to 12 (and `POT_C_J` to 16) removes the last of the
four claimed disagreements.

---

### D15 — UNVERIFIED, but I expect it to bite. Staged PCD size.

`stage_dataset` defaults to `target_points=2_000_000` and writes ASCII PCD via
`np.savetxt`. On disk, the existing staging of the user's four fragments produced
1.75 M points → **97 MB per fragment**, 390 MB for four. At the plan's 2 M default and a
100-fragment collection that is roughly 11 GB of scratch per run, written line-by-line
from Python. The reference pots are far smaller (Pot A piece 1: 100 064 points, 5.4 MB)
so this will not show up until Task 9.

Worth measuring before the sweep, and worth considering binary PCD (`DATA binary`) —
PCL reads it and it is ~7× smaller and far faster to write.

---

## 3. Gaps

Requirements for the stated goal that no task covers.

### G1 — Nothing makes the fragment count arbitrary. This is the actual goal, and it is unaddressed.

The design opens with *"the fragment count is not known in advance — it must work for
four fragments and for a hundred, with no rework in between."* Every task in the plan
operates on a fixed, hand-written preset.

`class/data_path.h` is **6 790 lines** of 30 mutually exclusive `#ifdef` blocks. Each
enumerates, by hand, `SHARD_NUMBER` plus seven fully spelled string arrays
(`file_path`, `obj_path`, `axis_path`, `surface_in`, `surface_out`, `surface_fr`,
`gt_T_path`, `gt_graph_path`, `shard_on_off`). There is no globbing anywhere in the
reassembly. A new collection of *n* fragments means writing ~9n lines and recompiling —
and Task 8 makes this worse by adding ten more binaries.

Meanwhile the preprocessing side *does* glob (`mesh_processing.cpp:1706-1712`,
`edgeline_extraction.cpp:3172-3182`) but bakes the pot letter into a compile-time
`#define` (D3b).

Nothing in the plan generates `data_path.h`, replaces it with a runtime manifest, or
adds `SFSPP_POT_ID`. Without one of those, the deliverable is "Pot A reassembles" plus
a benchmark, and the headline requirement is untouched. This should be a task, and
arguably the *first* one after the harness — it also happens to be the cleanest fix for
D3.

Rough shape: emit `class/data_path.h` from the staged directory listing (a ~60-line
Python generator), or change the reassembly to read a manifest file at startup so the
binary stops being per-collection. The former is smaller and needs no C++ change beyond
deleting the presets you stop using.

### G2 — Mixed-pot reassembly is ignored, and it is the harder half of the published result.

The dataset ships `Transformation/ABC`, `ABCDE`, `ABCDE(10,5)`, `DE`, and `data_path.h`
ships `POT_A_B_C` (21), `POT_D_E` (60), `POT_E_I` (61), `POT_C_J` (23),
`POT_A_B_F_G_H` (42), `POT_A_B_C_D_E` (81), `POT_All` (148). The page publishes five
mixed experiments (Dishes 97.6 %, C+J 80.0 %, E+I 96.6 %, ABCDE 92.5 %, All 87.3 %).

This matters for the goal, not just for completeness: **a folder of raw scanned sherds
is not guaranteed to come from one vessel.** The "no joins found" outcome the design
introduces is the right instinct, but the mixed case is the general case, and
`NUM_MIXED_SHERD` already exists to express it. At minimum the design should say
explicitly that single-vessel input is an assumption, and the harness should not report
"no joins" when the real answer is "these are two pots".

Note the runtime cost before committing: paper Table IX gives the All-mixed experiment
at (k,b) = (10,5) → 42.7 h, (20,10) → 137 h. The single-pot numbers are at (5,3). Do not
gate on mixed results without deciding the budget.

### G3 — The user's fragments have neither rim nor base, and nothing tests that path.

`main.cpp:31-32`:

```cpp
//#define NO_RIM_INFO
#define NO_BASE_INFO
```

`NO_BASE_INFO` is **active upstream** (present in the first commit) and zeroes
`is_seg_base_` for every shard at `main.cpp:105-109`. So base detection has no effect on
reassembly at all — which makes the design's *"Base detection: sound"* row true but
irrelevant, and matches the paper's claim of being base-agnostic. Worth saying, so
nobody spends time there.

`NO_RIM_INFO` exists and is commented out. It is the exact switch for "fragments with no
rim", and no task uses it. Two consequences:

1. Task 5 (port the rim criterion) improves Pot A but does **nothing** for the user's
   own fragments, which have no rim to detect. The plan sequences the rim fix ahead of
   Task 9 as though it helps; it does not. Say so.
2. The honest test for the user's data is: does the pipeline reach its targets on a
   reference pot with `NO_RIM_INFO` defined? Paper Table X's "Non-enhanced" ablation
   suggests it degrades but does not collapse. That run is cheap on Pots A/B/C and would
   tell you, before Task 9, whether a rimless four-fragment collection is even in scope.

Add it as a task. It is one `#define`, one rebuild, three small pots.

### G4 — The user's fragments are 4–9× larger than any reference fragment, and their units are unverified.

Measured bounding boxes:

```
set              file                 verts     max extent   diagonal
test_fragments   FY234007_reduced.ply   615160     546.31      687.86
test_fragments   FY234021_reduced.ply   526278     470.28      635.07
test_fragments   FY234094_reduced.ply   670739     679.17      790.49
test_fragments   FY234104_reduced.ply   354074     370.91      486.43
potB             FY234007.ply          7689470     546.59      688.27
potB             FY234021.ply          6578460     470.51      624.95
potB             FY234094.ply          8384209     679.60      785.44
potB             FY234104.ply          4425897     371.12      481.16
```

Reference fragments: median max extent 74.7, dataset maximum 210.6. The user's are
371–680 in the same nominal units. The PLY headers say `Created in RealityCapture`
(three files) and `by Geomagic Studio` (one) — photogrammetry output, which is
**unscaled unless a scale constraint was supplied at reconstruction**. So the coordinates
may not be millimetres at all.

Nothing in the plan checks this, and every millimetre threshold in the pipeline depends
on it. Before Task 9 runs, one physical measurement of one sherd resolves it. Without
that, a "no joins found" verdict on `potB` is uninterpretable — it could equally mean
"the correspondence window is 5 units on a 500-unit object".

Also: `test_fragments` and `potB` are the same four scans, decimated ~12×, but their
z-extents disagree (145.15 vs 89.76 for FY234021; 131.39 vs 91.64 for FY234094; 103.41
vs 71.47 for FY234104). The decimated copies carry stray points 30–55 units off-surface
in z. Clean them before using `test_fragments` for any geometric statistic — including
the extent that Task 7 would feed into `infer_scales`.

### G5 — `surface_fr` (fractured-surface points) is never produced, and the design does not mention it.

`data_path.h` declares `surface_fr[]` (`_Surface_0_FracturedSurfacePts.pcd`) for every
preset, and `POT_TEST` carries a comment saying `ReadPCD` tolerates its absence and
"leaves the fractured-surface term out of the optimisation". Our preprocessing produces
none, for any collection. So every one of our runs — including the ablations — is
missing an optimisation term that the authors' Pot A run also lacks (the ICCV data has
none either). That makes it a non-difference for the Pot A diagnosis, but it is an
unlisted difference from the *published* configuration and belongs in the design's
"what is already established" table.

### G6 — The harness's required "rendering of the result" has no task.

The design lists it as a required output per run and Task 9 Step 3 asks the worker to
judge *"whether the edge-on renders show a continuous wall or two shells crossing at an
angle"*. No task builds a renderer, and the reassembly's own viewer is interactive-only
(and the port notes say keyboard shortcuts are layout-dependent, which is why
`SFS_AUTOSAVE` exists). Either add the task or drop the requirement.

### G7 — Runtime budget is understated.

Task 8 Step 4: *"Expect the full sweep to run for a couple of hours."* ICCV 2021 Table 2,
on the authors' hardware, gives preprocessing alone at 47.8 min (D) + 37.2 min (E), and
reassembly at (b=3,k=5) 22.5 + 17.0 min. Pot I is a third pot of that size. Plus the
other seven pots and their preprocessing. Several hours is optimistic; a full day is
more realistic on a laptop. Sequence the sweep so partial results are usable (the plan
already orders small pots first — good), and checkpoint per pot so an interruption does
not lose the run.

### G8 — `Ground Truth Axes/` (141 files) is unused.

The dataset ships the authors' ground-truth axes separately from `Axes/`. `extract_axis`
was validated at 0.35° against `Axes/`, not against `Ground Truth Axes/`. If those are
the manually-corrected axes, they are the better validation target and a free way to
check the axis stage across all ten pots rather than Pot A alone. Worth one command
before deciding.

---

## 4. Confirmations

Claims I checked that held up. Do not re-examine these.

**C1 — The two headline results are real. Reproduced.** I implemented the plan's
`score_result` verbatim and ran it against the artefacts already on disk:

```
DatasetRef (authors' preprocessed Pot A)  ->  8 / 8   mean 3.739 mm
DatasetA   (our preprocessing, same Pot A) ->  0 / 8   mean 89.463 mm
```

Per-fragment errors on `DatasetA` are 68–147 mm — not marginal misplacement, total
failure. The diagnosis that the defect is inside preprocessing is well supported. (The
mean value itself is wrong in the plan — see D5.)

**C2 — Both runs used the same meshes, so the comparison is sound.**
`DatasetA/Mesh/Pot_A/Pot_A_Piece_01_Mesh.obj` and
`DatasetRef/Mesh/Pot_A/...` (a symlink) are both md5 `10075a5f…`, the ICCV file.

**C3 — Pot A geometry is identical between the ICCV and SfS++ collections**, though not
byte-identical as the plan states. The two files differ in exactly one line:

```
800503c800503
< mtllib 18-0702-03-01.mtl
---
> mtllib Pot_A_Piece_01_Mesh.mtl
```

Same 100 064 `v`, 200 124 `f`, 600 372 `vt`, identical first vertex. `Pot_A_Piece_1_T.txt`
and `Pot_A_simple_graph.txt` are identical between the two collections. The plan's
"byte-identical" is a small overstatement with no practical consequence; change the word
to "geometrically identical" and move on.

**C4 — `Mesh.zip` and `Point.zip` exist, the sizes are right, and `Point.zip` does
contain what the plan assumes. Verified without downloading it.**

Both Drive IDs in `SfSpp_preprocessing/download.sh` resolve. Exact sizes from a
one-byte ranged request:

```
Mesh.zip   185,832,298 bytes  (177.2 MiB)
Point.zip  2,663,152,436 bytes (2.48 GiB)
```

I then fetched only the last 2 MB of `Point.zip` and parsed its ZIP central directory
(154 entries, 11 directories + 143 `.pcd`):

```
Point/Pot_A  8    Point/Pot_F  6
Point/Pot_B  9    Point/Pot_G  7
Point/Pot_C  4    Point/Pot_H 11
Point/Pot_D 29    Point/Pot_I 27
Point/Pot_E 31    Point/Pot_J 11     total 143 pcd, 3.73 GiB uncompressed
```

Layout is `Point/Pot_<X>/Pot_<X>_Piece_<NN>_Point.pcd` — exactly what
`getPointDatasetPath(potID)` expects. **The plan's assertion is safe to rely on** for
A, B, C, D, E, G, H. It is *not* safe for F, I and J: `Point.zip` contains only the
evaluated subset (F 6, I 27, J 11), while `POT_F`/`POT_I` list 7 / 30 breakline files.
Those extras have no published point cloud, so the full chain cannot be run for them
from published raw inputs without first disabling those fragments — which you should be
doing anyway per D2.

Cheap confirmation command, for the record:

```bash
UA="Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 Chrome/126.0 Safari/537.36"
ID=1bBuqBIFnOQug9O0aYTkYbDwzIAaGMTsx; SIZE=2663152436
curl -s -A "$UA" -H "Range: bytes=$((SIZE-2000000))-$((SIZE-1))" \
  "https://drive.usercontent.google.com/download?id=$ID&export=download&confirm=t" \
  -o /tmp/point_tail.bin
# then parse the ZIP64 EOCD + central directory from the tail
```

**C5 — The pipeline-order claim is correct and the upstream README is wrong.**
`edgeline_extraction.cpp:2894-2897` builds the axis path and checks existence; the axis
drives rim classification (2928-2961) and `detectBaseType` (2991). The preprocessing
`README.md:74-99` lists axis extraction third, after EdgeLineExtraction. One nuance to
add: the axis is **optional**, not required — with no axis file the code prints "Axis
information is not available." (2964) and falls back to the geometric rim verdict with
`base_type = 0`. So a wrong order degrades silently rather than failing.

**C6 — All fifteen `SFSPP_*` environment variables exist and are read.** No extras, none
missing. `SFSPP_DATASET_ROOT`/`SFSPP_TEMP_ROOT` at `data_path.h:54-55`;
`SFSPP_NORMAL_NEIGHBORS` (5), `SFSPP_SAMPLING_RADIUS` (0.6), `SFSPP_MIN_CLUSTER` (50),
`SFSPP_RG_NEIGHBORS` (10), `SFSPP_SMOOTHNESS_DEG` (4.5), `SFSPP_CURVATURE` (1.5) in
`mesh_processing.cpp`; the seven `SFSPP_BASE_*` at `edgeline_extraction.cpp:2707-2713`.
Parse failures fall back to defaults silently. The reassembly side reads exactly four:
`SFS_AUTOSAVE`, `SFS_LENGTH_SCALE`, `SFS_DATA_ROOT`, `SFS_AXIS_COST_TOLERANCE`.

**C7 — The 14-character truncation is real.** `mesh_processing.cpp:338` (not 318):
`fileNameOnly = filePath.stem().string().substr(0, 14);`. It does **not** throw on
shorter stems (`pos = 0` is always valid, `count` is clamped) — the failure mode is
silent truncation and collision for stems *longer* than 14. So the constraint is "≤ 14
chars, or unique in the first 14", which `Pot_A_Piece_01` satisfies. Note two other
naming schemes coexist: `substr(0, find("_Point"))` at `mesh_processing.cpp:1720`, and
`substr(0, find("_Mesh"))` at `edgeline_extraction.cpp:2464, 2894` — the latter is what
builds the axis filename, so a rename breaks the schemes inconsistently.

**C8 — `MakeCorWOTree`'s 5 mm / 20 mm window is as described.**
`feature_matching.cpp:528, 531`. Under 5 mm a breakline point joins the overlap region
unconditionally; between 5 and 20 mm only if the surface normals agree and the
connecting vector is within 30° (`dir >= 0.86`). `OverlapCheck_3d` then returns false if
either region has fewer than 10 points (`:642-644`) — which is the mechanism the design's
rationale describes.

**C9 — The `check_base_and_rim.m` rim criterion is quoted correctly.** All five numbers
confirmed at `AxisExtraction/check_base_and_rim.m:13-16, 116, 121-122`. Three semantic
notes for the port: it is `std`, not variance, despite the variable names
`mean_r_err`/`mean_h_err`; the threshold is on `r`, not `r²`
(`compute_l2p_distance.m:12` returns a true perpendicular distance); and `h` is the raw
projection `dir · p` with no origin subtraction (harmless for `std` and `diff`). The
argument order is `[direction; position]`, i.e. `vt(4:6)` is the direction. There is also
a rim tie-break at lines 123-130 that the plan does not mention.

The **base** half was ported and matches: `edgeline_extraction.cpp:2696-2798`
`detectBaseType`, all seven thresholds, the 12×30° histogram, the `atan2d(x,y)` argument
order, and the N−1 std. The rim half genuinely is unported — the current C++ rim decision
is `isBreaklineSegARim` (`:2250`) plus a `100*stdev/mean < 50` rule at `:2947`. Task 5's
premise is sound.

One caveat: `check_base_and_rim.m` calls `align_point_cloud` (line 57) and
`read_axis.m:14` calls `load_root_dir` with a hard-coded `'D:/SFS_BB_temp/Plt_A'` —
neither function exists in the repository. The MATLAB rim path is not runnable as
shipped, so there is nothing to diff the port against except the criterion itself.

**C10 — Task 5's edit-site line numbers are right where it matters.**
`isBreaklineSegARim` is called at `edgeline_extraction.cpp:2927` and `:3051` — exactly as
the plan says, and those are the only two call sites. Its definition is at `:2250`;
signature `bool isBreaklineSegARim(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_BreakLineSeg,
string sampledDataWithNormals, int segNo)`. One correction: `processFragmentData` is at
`:2831`, not ~2650. The helper block does run 2652–2829 and a comment block ends at 2650,
so 2650 remains a valid insertion point — it is just ~180 lines above
`processFragmentData`, not adjacent to it.

**C11 — `extract_axis`'s interface is as documented.** `tools/extract_axis.cpp:236-245`:
`extract_axis <Surface_0.xyz> <Surface_1.xyz> <out_Axis.xyz> [ransac_iters]`, inner
surface first, iteration count optional (default 1000). It errors and writes nothing when
an input is missing (exit 1). One correction for the harness: the output has **1 to 10
lines**, not a fixed count — one per surviving candidate after 10° de-duplication and the
`SFS_AXIS_COST_TOLERANCE` filter. Real files in `DatasetA/Axes/` have 1, 1, 1, 1, 1, 1, 2, 3
lines. `run_pipeline`'s `axes[name] = sum(1 for line in ...)` is therefore measuring
something meaningful, but the plan should say what.

**C12 — Every log string `parse_stats` depends on is really emitted.**
`"Total number : "` `main.cpp:211`; `"Total number pruned : "` `main.cpp:217`;
`"Number of clusters is equal to "` `mesh_processing.cpp:1433`; `"... segment(s), ...
points, rim ..."` `edgeline_extraction.cpp:3005-3007`. Only the `Score :` regex needs
anchoring (D9).

**C13 — `Pot_A` breakline segment counts already match the reference exactly.** Ours vs
theirs: 5/5, 3/3, 4/4, 4/4, 4/4, 2/2, 3/3, 3/3. Point counts within 0–3.2 %. The design's
"breakline extraction: close" is right, and understated — see D7 for how to exploit it.

**C14 — Pot A's median wall thickness is 3.6731 mm**, so `REFERENCE_THICKNESS_MM = 3.69`
is accurate to 0.5 %. Seven of eight fragments sit in a tight 3.566–3.939 band. (The
*extent* constant is a different story — D6.)

---

## 5. Recommended changes

Concrete edit list. Design document first, then plan.

### `2026-08-21-arbitrary-fragments-design.md`

1. **Replace the "Definition of done" table** with the corrected ten-pot table from D2,
   expressed as Sherd Accuracy percentages with the integer fraction alongside. Delete
   the paragraph about F, G, I, J being ungated; replace it with a note that the excluded
   fragments are the all-zero rows of `Pot_<X>_simple_graph.txt` and are disabled via
   `shard_on_off`.
2. **Add a "Metric" subsection** stating that success is `CountResult`'s Sherd Accuracy
   (relative pose, 20° / 50 mm, ≥1 correct edge over the GT adjacency graph), that the
   binary already computes it, and that the Python Kabsch scorer is a *diagnostic* for
   millimetre displacement, not the gate. Delete "mean placement error < 20 mm after one
   global rigid alignment".
3. **Fix three numbers in "What is already established"**: mean error 3.74 mm (pooled) or
   4.37 mm (mean of per-fragment means) — say which; rim detection finds **1** rim, not 0;
   breakline point counts are within **+3.2 %**, not +5–10 %.
4. **Add two rows to that table**: `NO_BASE_INFO` is active upstream so base flags never
   reach reassembly (making the base-detection row informational); `surface_fr` is
   produced by nobody, for any collection.
5. **Rewrite the ablation table** per D7: correct Run C's axes column to "theirs
   (re-derived)"; split Run B into B1 (rim flags only, by rewriting the header byte) and
   B2 (whole breakline files); note that breaklines are downstream of surfaces so Run C
   leaves a mixed state unless `EdgeLineExtraction` is re-run.
6. **Rewrite "Scale inference"** per D6. Replace the "must come out at 1.0" criterion with
   an envelope criterion: measured across A–J, thickness spans 2.80–8.20 mm (2.93×) and
   extent 35.6–93.4 mm (2.63×) *with fixed thresholds and published accuracy throughout*,
   so scale 1.0 must hold for anything inside that envelope and inference applies only
   outside it. Correct `REFERENCE_EXTENT_MM`: Pot A's median max-extent is 93.42 and its
   median diagonal 126.50; 119.0 is the largest single fragment. Drop the 7.26 mm figure.
   Correct the user-fragment figures: measured thickness is **37.58**, ratio **10.2×**,
   not 8.25 / 2.2×; extent ratio 5.4×. Add that `SFS_LENGTH_SCALE` moves 12 of ~57
   length constants (D11), and that inferring thickness from `Surface_0`/`Surface_1`
   calibrates on the output of a suspect component.
7. **Add a section on the unknown scale of the user's own fragments** (G4), with the
   measured extents, the RealityCapture/Geomagic provenance, and the requirement for one
   physical measurement before Task 9 is interpretable. Note the stray z-points in
   `test_fragments`.
8. **Add "arbitrary fragment count" as a first-class component** (G1). It is the stated
   goal and currently nothing addresses it.
9. **Add a paragraph on mixed collections** (G2) — at minimum, state single-vessel input
   as an explicit assumption and note that `NUM_MIXED_SHERD` exists.
10. **Add `NO_RIM_INFO` to the rimless discussion** (G3), and state plainly that Task 5
    does not help the user's own fragments.

### `2026-08-21-arbitrary-fragments.md`

11. **Global Constraints:** change "Python 3.13" to 3.14.7; state the test-runner
    decision (venv, `--break-system-packages`, or `unittest`) — `pytest` is not installed
    (D8b); change "byte-identical" to "geometrically identical, differing only in the
    `mtllib` line" (C3); replace the face-expansion rule with the conditional form (D4);
    replace the success threshold with the SA definition (D1); add `Mesh.zip`
    185,832,298 B and `Point.zip` 2,663,152,436 B with the verified per-pot content and
    the F/I/J caveat (C4); add the `git add -A` warning for the preprocessing repo (D12).
12. **New Task 0:** commit the pending `data_path.h` / `edgeline_extraction.cpp` changes
    in the preprocessing repo on their own; extend its `.gitignore` with `Dataset*/`
    and `Temp*/` (D12).
13. **New Task 0b (or fold into 1):** generate `class/data_path.h` from a staged
    directory listing, and add `SFSPP_POT_ID` to the preprocessing `data_path.h` so one
    preprocessing binary serves every collection (G1, D3b). This replaces Task 8 Step 1
    entirely and is the prerequisite for the sweep.
14. **Task 1:** note that `write_obj` emits no `vt`, and that this is what makes result
    OBJs non-expanded (D4) — a comment in the file, so nobody "fixes" it later without
    understanding the consequence. Consider binary PCD and measure the ASCII size before
    Task 9 (D15).
15. **Task 2:** change the expected gate to `8 / 8 mean 3.74 (± 0.05)`, or redefine
    `mean_error` as the mean of per-fragment means and keep 4.37 (D5). Make
    `expand_faces` conditional and raise on a mismatch instead of truncating (D4).
16. **New task between 2 and 3 — wire the published metric into the binary.** Extend the
    `SFS_AUTOSAVE` path to run the `vis.fine_` reconstruction and `CountResult` +
    `SaveAcc` before exiting; parse `Result/1. Acc.txt` in the harness. Set
    `shard_on_off = false` for D piece 22, F piece 7, I pieces 28–30 (and the same five in
    `POT_All`); fix `POT_J` `SHARD_NUMBER` 19 → 12 and `POT_C_J` 23 → 16 (D14); fix
    `POT_CTRL`'s `gt_T_path` / `gt_graph_path` to the `Ground Truth/` layout (D13). This
    is the task the whole benchmark depends on.
17. **Task 3:** clear `Result/` at the start of every `run_pipeline` (D8d); make the
    "copy eight meshes to a scratch folder" step the instruction rather than a fallback
    (D8e); reconcile the declared dict-returning interface with the list-returning
    implementation, or key the lists by fragment name (D9); anchor the `Score :` regex
    (D9); document that `axes[name]` counts 1–10 candidate axes, not a fixed number (C11).
18. **Task 4:** fix `'$I'` → `'$R'` (D8a); derive the fragment list from the directory
    rather than `range(1, 9)` (D8c); add Run B1 (header-byte-only substitution) per D7;
    correct Run C's description; state whether `EdgeLineExtraction` is re-run in Run C.
19. **Task 5:** correct the reference vector to `2, 1, 1, 1, 0, 1, 0, 1` (D10); correct
    `processFragmentData` to line 2831 (C10); correct the current-state claim from "0 rims"
    to "1 rim, piece 04" (D7); commit only the rim hunk, not the whole file (D12).
20. **Task 7:** rewrite against the envelope criterion (D6); fix
    `REFERENCE_EXTENT_MM`; replace the two tautological unit tests with a test that
    `infer_scales` returns 1.0 for *each* of the ten reference pots under the envelope
    rule; fix `_nearest_distances` to expand the search ring until
    `best <= ring * cell` (D6); state that `SFS_LENGTH_SCALE` covers 12 of ~57 constants
    and decide what to do about the other 45 (D11). Move Task 7 **after** the diagnosis
    and derive thickness from the mesh rather than from `Surface_0`/`Surface_1`, to
    break the circularity of calibrating on the output of a suspect component (D6).
21. **Task 8:** rewrite on top of the generated `data_path.h` (item 13). Until then it
    cannot run — see D3 (a)–(d). Resolve the `_Mesh` / `_Mesh_DS` suffix per file in both
    `_pot_source` and `score_result`. Gate on all ten pots using the corrected targets
    (D2). Revise the runtime estimate (G7) and checkpoint per pot.
22. **New task before 9 — the rimless run.** Build with `NO_RIM_INFO` and sweep Pots A,
    B, C. This is the only evidence that says whether the user's rimless fragments are in
    scope at all, and it costs one rebuild and three small pots (G3).
23. **Task 9:** add a step that establishes the physical scale of the `potB` scans before
    running, and state that a "no joins" verdict is uninterpretable without it (G4). Note
    that `test_fragments` carries stray z-outliers and should not be used for the extent
    statistic. Add the `--pot-id A` requirement explicitly and explain why (`POT_TEST` is
    hard-wired to `Pot_A_Piece_NN`), rather than leaving it as an unexplained argument.
    Replace the `O(n²)` pairwise-distance snippet in Step 2: on the existing `potB`
    result (615 160 vertices per fragment) `M[i][:,None,:] - M[j][None,::7,:]` at
    stride 17/7 allocates a `36186 x 5170 x 3` float64 array — **4.5 GB per pair**, six
    pairs. It will thrash or MemoryError as written. Use the same voxel-grid nearest-
    neighbour helper Task 7 already needs.
24. **Either build the renderer or drop it** from the design's required outputs and from
    Task 9 Step 3 (G6).
25. **Add one cheap task:** validate `extract_axis` against `Ground Truth Axes/` across
    all ten pots, not `Axes/` on Pot A alone (G8).

---

## Appendix — how to reproduce the measurements in this review

```bash
# D2: reconcile the published percentages with the ground-truth graphs
GT="Dataset/SfS_pp/Ground Truth"
python3 -c "
import numpy as np, glob, re
for p in 'ABCDEFGHIJ':
    G=np.loadtxt(f'$GT/Pot_{p}_simple_graph.txt')
    keep=[i for i in range(G.shape[0]) if G[i].sum()>0]
    print(p, G.shape[0], 'evaluated', len(keep), 'edges', int(G[np.ix_(keep,keep)].sum())//2)
"

# D4: face expansion is conditional
grep -c '^vt ' Dataset/SfS_pp/Mesh/Pot_A_Piece_01_Mesh.obj          # 600372
grep -c '^vt ' ../SfSpp_preprocessing/DatasetB/Mesh/Pot_A/Pot_A_Piece_01_Mesh.obj  # 0
grep -c '^v '  ../SfSpp_preprocessing/DatasetB/Result/*_Top_1_OBJ_1.obj            # 615160 = source verts

# D5 / C1: the plan's scorer against the artefacts on disk
#   (see the transcript; score_result copied verbatim from the plan)

# D10: reference rim flags
for i in 01 02 03 04 05 06 07 08; do
  sed -n 2p "Dataset/SfS_pp/Breaklines/Pot_A_Piece_${i}_Breakline_0.pcd"
done

# D8b: pytest
/opt/homebrew/bin/python3 -m pytest --version
/opt/homebrew/bin/python3 -m pip install --user --dry-run pytest
```

Thickness and extent measurements (D6, G4) used a voxel-grid nearest-neighbour search
with ring expansion until `best <= ring * cell`, validated against exact chunked brute
force on all eight Pot A fragments (max absolute difference 3.6e-15). Median of 3000
sampled `Surface_0` points per fragment, seed fixed. Scripts are in the session
scratchpad (`thickness.py`, `bbox.sh`, `aggregate.py`, `ply_bbox.py`); total runtime
about 46 s for all 163 reference surface pairs and 164 meshes. The user-fragment
thickness (37.58) used the same routine on
`SfSpp_preprocessing/DatasetB/Surfaces/Pot_A/Pot_A_Piece_0{1,2,3,4}_Surface_{0,1}.xyz`.
