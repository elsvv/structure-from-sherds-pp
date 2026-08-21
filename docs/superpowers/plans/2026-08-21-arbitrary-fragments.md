# Arbitrary Fragment Reassembly — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Revision 2 (2026-08-21).** Rewritten after independent review
(`docs/superpowers/notes/2026-08-21-plan-review.md`). Five blocking defects are
fixed here: the plan now gates on the metric the papers publish and the binary
already computes; the target table is corrected and covers all ten pots; the
sweep runs on a generated preset instead of ten binaries that could not find
their files; face expansion is decided from the data instead of asserted; and the
length-scale component is rebuilt around an envelope rather than a ratio that
contradicted its own test. Task numbering changed — the old Task *n* is named
where it moved.

**Goal:** Reassemble a folder of raw scanned pot fragments without per-vessel
parameter tuning, verified against the authors' published Sherd Accuracy for all
ten SfS++ pots.

**Architecture:** The reassembly stops being per-collection: `class/data_path.h`
is generated from the staged directory listing, and the preprocessing binaries
take their pot id from the environment. A Python harness (`sfsbench`) then drives
the C++ binaries end to end, and reads the accuracy the binary itself writes.
Diagnosis identifies which preprocessing defect is fatal; that defect is fixed;
the millimetre thresholds are addressed last, once segmentation is cleared.

**Tech Stack:** Python 3.14.7 (numpy 2.4.1 only), C++17, PCL 1.15, CGAL 6.2,
Ceres 2.2, CMake, Homebrew.

## Global Constraints

- Three repositories, all on branch `macos-port`:
  - `/Users/vaceslaveliseev/@dev/structure-from-sherds-pp` — reassembly (`build/SfSpp`, `build/SfSpp_ctrl`, `build/extract_axis`). Remote: `origin` = user's fork. **Harness lives here.**
  - `/Users/vaceslaveliseev/@dev/SfSpp_preprocessing` — preprocessing (`build/MeshPreprocessing`, `build/EdgeLineExtraction`). **No remote yet** — commit locally only.
  - `/Users/vaceslaveliseev/@dev/structure-from-sherds` — ICCV 2021 repo, holds the older `ICCV Data/`. **Read-only, and superseded as the reference by `Dataset/SfS_pp` below.**
- Python is `/opt/homebrew/bin/python3` (**3.14.7**). Only `numpy` may be
  imported. `scipy`, `trimesh`, `open3d` are NOT installed and must NOT be added.
- **Test runner: a venv at `tools/.venv`.** `pytest` is *not* installed for the
  Homebrew interpreter and `pip install --user pytest` fails with
  `externally-managed-environment` (PEP 668). Create it once, before Task 3:

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
/opt/homebrew/bin/python3 -m venv --system-site-packages .venv
.venv/bin/pip install -q pytest
echo "tools/.venv/" >> ../.gitignore
```

  `--system-site-packages` inherits numpy 2.4.1 rather than rebuilding it. Every
  `pytest` invocation below is `.venv/bin/pytest`. (Fallbacks if the venv is
  unavailable: `pip install --break-system-packages pytest`, or rewrite the
  tests for stdlib `unittest`.)
- Rebuild after any C++ change: `cmake --build build -j10` in that repository.
- **Reference dataset: the SfS++ collection**, `structure-from-sherds-pp/Dataset/SfS_pp`
  (already downloaded, 713 MB, gitignored). Ten pots A–J with `Mesh/`, `Axes/`,
  `Breaklines/`, `Surfaces/`, `Ground Truth/` and `Transformation/`.
  - It supersedes the ICCV 2021 data. `Pot_A_Piece_01_Mesh.obj` is
    **geometrically identical** in both — same 100 064 `v`, 200 124 `f`,
    600 372 `vt`, identical first vertex — but not byte-identical: the two files
    differ in exactly one line, `mtllib`.
  - **Mesh filenames carry two suffixes.** Pots A and B ship `_Mesh.obj`; pots
    C–J ship `_Mesh_DS.obj` — 147 of the 164 files. Any glob must accept both.
    Conversely `mesh_processing.cpp:1722` builds its mesh filename as
    `<stem>_Mesh.obj` unconditionally, so **staging must always write
    `_Mesh.obj`**, whatever the source was called.
  - The ICCV 2021 folder `structure-from-sherds/ICCV Data` is used only where a
    step explicitly says so.
- Ground truth: `Dataset/SfS_pp/Ground Truth/Pot_<X>_Piece_<n>_T.txt` — **note the
  space** in `Ground Truth`, and that the piece number is **not** zero-padded here
  while mesh filenames are (`Pot_A_Piece_01_Mesh.obj`).
- Raw pipeline inputs for the reference pots are published separately and are
  fetched by `SfSpp_preprocessing/download.sh`. Both Drive ids resolve; sizes
  verified by ranged request: `Mesh.zip` **185,832,298 B** (177.2 MiB),
  `Point.zip` **2,663,152,436 B** (2.48 GiB). `Point.zip`'s central directory was
  parsed without downloading the body: 143 `.pcd` under
  `Point/Pot_<X>/Pot_<X>_Piece_<NN>_Point.pcd`, exactly the layout
  `getPointDatasetPath(potID)` expects, distributed
  A 8, B 9, C 4, D 29, E 31, F **6**, G 7, H 11, I **27**, J **11**.
  **Caveat:** for F, I and J the archive holds only the *evaluated* subset, while
  `POT_F`/`POT_I` list 7 / 30 breakline files. The extra fragments have no
  published point cloud, so the full chain cannot be run for them from published
  raw inputs — which is fine, because they are the ones being disabled anyway
  (Task 2).
- **Success metric: Sherd Accuracy**, as defined in SfS++ §V-B and implemented at
  `class/data_structure.cpp:920 CountResult` — relative pose `T_i⁻¹ T_j` against
  the ground-truth relative pose, thresholds 20° (0.35 rad) and 50 mm, over the
  adjacency graph, a sherd counting as correct with **at least one** correct edge.
  The binary writes it to `Result/1. Acc.txt`. The Python Kabsch scorer is a
  **diagnostic** reporting millimetre displacement; it is never the gate.
  *(Removed: "mean placement error < 20 mm after one global rigid alignment" —
  that metric and that threshold are in neither paper.)*
- **Face expansion is conditional, not universal.** VTK duplicates points on load
  only when the input OBJ carries per-corner `vt` indices. The authors' meshes do
  (`Pot_A_Piece_01_Mesh.obj`: 600 372 `vt`), so their results come back with
  `3 × n_faces` vertices. Meshes written by our stager have no `vt`, so their
  results come back with `n_vertices`. Measured on the runs already on disk:

```
staged source DatasetB/Mesh/Pot_A/Pot_A_Piece_01_Mesh.obj  v=615160 f=1230314 vt=0
result        DatasetB/Result/..._Top_1_OBJ_1.obj          v=615160        (not 3f)
authors' mesh Dataset/SfS_pp/Mesh/Pot_A_Piece_01_Mesh.obj  v=100064 f=200124 vt=600372
result        DatasetRef/Result/..._Top_1_OBJ_1.obj        v=600372  (= 3 x 200124)
```

  Any comparison must **pick the layout from the vertex count** and raise on a
  mismatch. *(The previous revision stated expansion as unconditional, which is
  correct for reference runs and wrong for every run the harness itself stages.)*
- Pipeline order is fixed and non-obvious: **MeshPreprocessing → extract_axis → EdgeLineExtraction → SfSpp**. `EdgeLineExtraction` reads the axis while classifying rim segments, so the axis must exist before it runs. Upstream's README has the axis last; that is wrong. Note the axis is *optional*, not required — with no axis file the code prints "Axis information is not available." and falls back to the geometric rim verdict, so a wrong order degrades silently rather than failing.
- Commit after every task. **Never use `git add -A` in either repository.** In
  the pp repo `data/` and `ICCV Data/` are large and gitignored but scratch
  output is not; in the preprocessing repo `DatasetA/ DatasetB/ DatasetRef/
  TempA/ TempB/` are 1.5 GB of untracked, un-ignored artefacts until Task 0
  fixes the ignore file.

---

### Task 0: Repository hygiene in the preprocessing repo

*(New. Blocks every later commit that touches that repo.)*

`SfSpp_preprocessing` has uncommitted work — the `SFSPP_DATASET_ROOT` /
`SFSPP_TEMP_ROOT` plumbing in `data_path.h` and `edgeline_extraction.cpp` — that
the entire plan depends on, and 1.5 GB of un-ignored scratch directories. Left
alone, the first `git add edgeline_extraction.cpp` in a later task silently folds
unrelated plumbing into that commit, and `data_path.h` never gets committed at all.

**Files:**
- Modify: `/Users/vaceslaveliseev/@dev/SfSpp_preprocessing/.gitignore`

- [ ] **Step 1: Confirm what is pending**

```bash
cd /Users/vaceslaveliseev/@dev/SfSpp_preprocessing
git status --porcelain
git diff --stat
```

Expected: `M data_path.h`, `M edgeline_extraction.cpp`, and five untracked
`Dataset*/` `Temp*/` directories.

- [ ] **Step 2: Extend the ignore file**

Append to `.gitignore` (it currently ignores `Dataset/` and `Temp/` only, which
does not cover the suffixed working copies):

```gitignore
Dataset*/
Temp*/
```

- [ ] **Step 3: Commit the pending plumbing on its own**

```bash
cd /Users/vaceslaveliseev/@dev/SfSpp_preprocessing
git add .gitignore data_path.h edgeline_extraction.cpp
git commit -m "Make the dataset and scratch roots configurable

SFSPP_DATASET_ROOT and SFSPP_TEMP_ROOT replace the hardwired container mount
points so the same binaries can be pointed at another collection without a
rebuild; defaults are unchanged. Also ignore the suffixed working copies
(DatasetA/, TempB/, ...) that the diagnosis runs produce."
git status --porcelain
```

Expected: clean tree, no untracked `Dataset*/`.

---

### Task 1: Make the fragment count arbitrary

*(New — this is the stated goal of the whole design and nothing addressed it.
It also replaces the old Task 8 Step 1, which added ten per-pot binaries that
could not have found their files.)*

Today `class/data_path.h` is 6 790 lines of 30 mutually exclusive `#ifdef`
blocks, each spelling out `SHARD_NUMBER` plus seven string arrays by hand; there
is no globbing anywhere in the reassembly. On the preprocessing side the code
*does* glob but bakes the pot letter into `#define POT_A`, so one binary serves
one pot. Both limits are removed here.

`class/data_path.h` is included by `main.cpp` alone and `SHARD_NUMBER` appears in
no other file, so regenerating the header recompiles one translation unit —
measured at ~5 s on this machine. "Generate, rebuild, run" is therefore an
ordinary harness step, not a porting job.

**Files:**
- Create: `tools/sfsbench/preset.py`
- Create: `tools/sfsbench/tests/test_preset.py`
- Modify: `class/data_path.h` — add the `POT_GEN` hook
- Modify: `/Users/vaceslaveliseev/@dev/SfSpp_preprocessing/data_path.h` — `SFSPP_POT_ID`
- Modify: `CMakeLists.txt` — a `SfSpp_gen` target defining `POT_GEN`

**Interfaces:**
- Produces:
  - `discover(dataset_root: str, pot_ids: list[str]) -> list[dict]` — one record
    per fragment, in reassembly order, each `{"pot": str, "name": str,
    "breakline": str, "mesh": str, "axis": str, "surface_in": str,
    "surface_out": str, "surface_fr": str}`. Paths are absolute. Resolves the
    `_Mesh.obj` / `_Mesh_DS.obj` suffix per file.
  - `excluded_pieces(gt_graph_path: str) -> set[int]` — the 1-based indices of
    all-zero rows, i.e. the fragments the paper excludes from evaluation.
  - `generate(fragments: list[dict], out_header: str, gt_root: str | None) -> None`
    — writes `class/data_path_generated.h`.

- [ ] **Step 1: Write the failing test**

```python
# tools/sfsbench/tests/test_preset.py
import os, tempfile
import numpy as np
from sfsbench.preset import discover, excluded_pieces, generate


def _touch(path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    open(path, "w").close()


def _fake_collection(root, pot, n, mesh_suffix="_Mesh.obj"):
    for i in range(1, n + 1):
        name = f"Pot_{pot}_Piece_{i:02d}"
        _touch(os.path.join(root, "Breaklines", f"Pot_{pot}", f"{name}_Breakline_0.pcd"))
        _touch(os.path.join(root, "Mesh", f"Pot_{pot}", f"{name}{mesh_suffix}"))
        _touch(os.path.join(root, "Surfaces", f"Pot_{pot}", f"{name}_Surface_0.xyz"))
        _touch(os.path.join(root, "Surfaces", f"Pot_{pot}", f"{name}_Surface_1.xyz"))
        _touch(os.path.join(root, "Axes", f"{name}_Axis.xyz"))


def test_discover_finds_every_fragment_in_order():
    with tempfile.TemporaryDirectory() as d:
        _fake_collection(d, "A", 12)
        frags = discover(d, ["A"])
        assert len(frags) == 12
        assert [f["name"] for f in frags][:2] == ["Pot_A_Piece_01", "Pot_A_Piece_02"]
        assert frags[-1]["name"] == "Pot_A_Piece_12"


def test_discover_accepts_the_downsampled_mesh_suffix():
    with tempfile.TemporaryDirectory() as d:
        _fake_collection(d, "C", 4, mesh_suffix="_Mesh_DS.obj")
        frags = discover(d, ["C"])
        assert len(frags) == 4
        assert frags[0]["mesh"].endswith("_Mesh_DS.obj")


def test_discover_spans_several_pots_for_a_mixed_collection():
    with tempfile.TemporaryDirectory() as d:
        _fake_collection(d, "A", 3)
        _fake_collection(d, "B", 2)
        frags = discover(d, ["A", "B"])
        assert [f["pot"] for f in frags] == ["A", "A", "A", "B", "B"]


def test_excluded_pieces_are_the_all_zero_rows():
    with tempfile.TemporaryDirectory() as d:
        g = np.array([[0, 1, 0, 0],
                      [1, 0, 1, 0],
                      [0, 1, 0, 0],
                      [0, 0, 0, 0]])
        p = os.path.join(d, "g.txt")
        np.savetxt(p, g, fmt="%d")
        assert excluded_pieces(p) == {4}


def test_generated_header_declares_the_right_count():
    with tempfile.TemporaryDirectory() as d:
        _fake_collection(d, "A", 5)
        out = os.path.join(d, "data_path_generated.h")
        generate(discover(d, ["A"]), out, gt_root=None)
        text = open(out).read()
        assert "#define SHARD_NUMBER 5" in text
        assert "#define NUM_MIXED_SHERD 1" in text
        assert text.count("_Breakline_0.pcd") == 5
        assert "shard_on_off" in text


def test_generated_header_disables_excluded_fragments():
    with tempfile.TemporaryDirectory() as d:
        _fake_collection(d, "F", 7)
        gt = os.path.join(d, "gt")
        os.makedirs(gt)
        g = np.ones((7, 7), dtype=int)
        np.fill_diagonal(g, 0)
        g[6, :] = 0
        g[:, 6] = 0                       # piece 7 excluded, as in the real Pot F
        np.savetxt(os.path.join(gt, "Pot_F_simple_graph.txt"), g, fmt="%d")
        out = os.path.join(d, "data_path_generated.h")
        generate(discover(d, ["F"]), out, gt_root=gt)
        body = open(out).read()
        on_off = body.split("shard_on_off")[1]
        assert on_off.count("true") == 6
        assert on_off.count("false") == 1
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
.venv/bin/pytest sfsbench/tests/test_preset.py -v
```

Expected: FAIL — `ModuleNotFoundError: No module named 'sfsbench.preset'`

- [ ] **Step 3: Write the generator**

```python
# tools/sfsbench/preset.py
"""Generate class/data_path_generated.h from a staged collection.

The reassembly enumerates every fragment path in a hand-written #ifdef block, so
a new collection of n fragments used to mean writing ~9n lines and recompiling.
data_path.h is included only by main.cpp and SHARD_NUMBER appears nowhere else,
so regenerating this header recompiles exactly one translation unit (~5 s).
"""
import glob, os
import numpy as np

MESH_SUFFIXES = ("_Mesh.obj", "_Mesh_DS.obj")


def _mesh_for(mesh_dir, name):
    for suffix in MESH_SUFFIXES:
        p = os.path.join(mesh_dir, name + suffix)
        if os.path.exists(p):
            return p
    raise FileNotFoundError(f"{name}: no mesh in {mesh_dir} "
                            f"(tried {', '.join(MESH_SUFFIXES)})")


def discover(dataset_root, pot_ids):
    """Every fragment of every named pot, in reassembly order.

    Layout is the one MeshPreprocessing writes: Breaklines/, Mesh/ and Surfaces/
    nested per pot, Axes/ flat. That is also what POT_CTRL and POT_TEST read;
    the shipped POT_A..POT_J presets read a flat layout instead and are left
    alone.
    """
    out = []
    for pot in pot_ids:
        bl_dir = os.path.join(dataset_root, "Breaklines", f"Pot_{pot}")
        mesh_dir = os.path.join(dataset_root, "Mesh", f"Pot_{pot}")
        srf_dir = os.path.join(dataset_root, "Surfaces", f"Pot_{pot}")
        axis_dir = os.path.join(dataset_root, "Axes")
        names = sorted(os.path.basename(p)[:14] for p in
                       glob.glob(os.path.join(bl_dir, f"Pot_{pot}_Piece_*_Breakline_0.pcd")))
        if not names:
            raise FileNotFoundError(f"no breaklines for pot {pot} under {bl_dir}")
        for name in names:
            out.append({
                "pot": pot,
                "name": name,
                "breakline":   os.path.join(bl_dir, f"{name}_Breakline_0.pcd"),
                "mesh":        _mesh_for(mesh_dir, name),
                "axis":        os.path.join(axis_dir, f"{name}_Axis.xyz"),
                "surface_in":  os.path.join(srf_dir, f"{name}_Surface_0.xyz"),
                "surface_out": os.path.join(srf_dir, f"{name}_Surface_1.xyz"),
                "surface_fr":  os.path.join(srf_dir, f"{name}_Surface_0_FracturedSurfacePts.pcd"),
            })
    return out


def excluded_pieces(gt_graph_path):
    """1-based indices of the fragments the paper excludes from evaluation.

    They are exactly the all-zero rows of Pot_<X>_simple_graph.txt, and exactly
    the ones with no _T.txt. Verified for all ten pots: D 22, F 7, I 28-30, J 9.
    """
    G = np.loadtxt(gt_graph_path)
    G = np.atleast_2d(G)
    return {i + 1 for i in range(G.shape[0]) if G[i].sum() == 0}


def _array(name, values, kind="string"):
    body = ",\n".join(f'\t{v}' for v in values)
    return f"{kind} {name}[SHARD_NUMBER] = {{\n{body}\n}};\n\n"


def generate(fragments, out_header, gt_root=None):
    """Write the preset. gt_root is Dataset/SfS_pp/Ground Truth, or None.

    With gt_root=None the ground-truth graph path is deliberately pointed at a
    file that does not exist: main.cpp reads it with a bare ifstream and no
    existence check, so num_raw comes out 0, GT_graph stays all-zero and the
    accuracy report is skipped. That is the correct behaviour for a collection
    with no ground truth, but it is silent, so the harness must not mistake a
    missing "1. Acc.txt" for a failed run.
    """
    pots = sorted({f["pot"] for f in fragments}, key=lambda p: [f["pot"] for f in fragments].index(p))
    n = len(fragments)

    def q(path):
        return '"' + path.replace("\\", "/") + '"'

    excluded = set()
    gt_T, gt_graph = [], []
    if gt_root:
        offset = 0
        for pot in pots:
            members = [f for f in fragments if f["pot"] == pot]
            graph = os.path.join(gt_root, f"Pot_{pot}_simple_graph.txt")
            gt_graph.append(q(graph))
            if os.path.exists(graph):
                excluded |= {offset + i for i in excluded_pieces(graph)}
            for i in range(1, len(members) + 1):
                # note: not zero-padded here, unlike the mesh filenames
                gt_T.append(q(os.path.join(gt_root, f"Pot_{pot}_Piece_{i}_T.txt")))
            offset += len(members)
    else:
        gt_T = [q(os.path.join("/nonexistent", f"{f['name']}_T.txt")) for f in fragments]
        gt_graph = [q("/nonexistent/simple_graph.txt")]

    lines = [
        "// GENERATED by tools/sfsbench/preset.py -- do not edit by hand.\n",
        f"// {n} fragment(s), pot(s) {', '.join(pots)}.\n\n",
        f"#define SHARD_NUMBER {n}\n",
        f"#define NUM_MIXED_SHERD {len(gt_graph)}\n\n",
    ]
    lines.append(_array("file_path",   [q(f["breakline"])   for f in fragments]))
    lines.append(_array("obj_path",    [q(f["mesh"])        for f in fragments]))
    lines.append(_array("axis_path",   [q(f["axis"])        for f in fragments]))
    lines.append(_array("surface_in",  [q(f["surface_in"])  for f in fragments]))
    lines.append(_array("surface_out", [q(f["surface_out"]) for f in fragments]))
    lines.append(_array("surface_fr",  [q(f["surface_fr"])  for f in fragments]))
    lines.append(_array("gt_T_path",   gt_T))
    body = ",\n".join(f"\t{g}" for g in gt_graph)
    lines.append(f"string gt_graph_path[{len(gt_graph)}] = {{\n{body}\n}};\n\n")

    flags = [f"\t{'false' if (i + 1) in excluded else 'true'},\t// {f['name']}"
             for i, f in enumerate(fragments)]
    lines.append("bool shard_on_off[SHARD_NUMBER] = {\n" + "\n".join(flags) + "\n};\n")

    os.makedirs(os.path.dirname(os.path.abspath(out_header)), exist_ok=True)
    with open(out_header, "w") as f:
        f.writelines(lines)
```

- [ ] **Step 4: Add the `POT_GEN` hook to `class/data_path.h`**

The generated header defines the same symbols the presets do, so it must be
mutually exclusive with them. Immediately **after** `string path = sfsDataRoot();`
and before the `POT_CTRL` block, insert:

```cpp
// POT_GEN: the preset is generated from the staged collection by
// tools/sfsbench/preset.py, so the fragment count stops being a compile-time
// constant written by hand. Only main.cpp includes this header and only main.cpp
// uses SHARD_NUMBER, so regenerating it recompiles one translation unit.
#ifdef POT_GEN
#include "data_path_generated.h"
#endif
```

and change the default-preset guard near the top so `POT_GEN` also suppresses
`POT_TEST`:

```cpp
#if !defined(POT_CTRL) && !defined(POT_GEN)
#define POT_TEST			// our own scans, see BUILD-macOS.md
#endif
```

Add `class/data_path_generated.h` to the pp repo's `.gitignore` — it is build
output, and committing it would put absolute scratch paths in history.

- [ ] **Step 5: Add the `SfSpp_gen` target**

In `CMakeLists.txt`, after the `SfSpp_ctrl` block:

```cmake
# --- SfSpp_gen --------------------------------------------------------------
# Same sources with POT_GEN: the preset comes from class/data_path_generated.h,
# written by tools/sfsbench/preset.py for whatever collection is staged. This is
# the binary the harness rebuilds and runs; SfSpp and SfSpp_ctrl stay as they are.
add_executable(SfSpp_gen ${PROJ_SRC})
target_include_directories(SfSpp_gen PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/class ${PCL_INCLUDE_DIRS})
target_link_libraries(SfSpp_gen PRIVATE Ceres::ceres ${PCL_LIBRARIES})
target_compile_definitions(SfSpp_gen PRIVATE ${PCL_DEFINITIONS} POT_GEN)
if(OpenMP_CXX_FOUND)
  target_link_libraries(SfSpp_gen PRIVATE OpenMP::OpenMP_CXX)
endif()
if(NOT MSVC)
  target_compile_options(SfSpp_gen PRIVATE -Wno-deprecated-declarations)
endif()
```

- [ ] **Step 6: Add `SFSPP_POT_ID` to the preprocessing repo**

In `/Users/vaceslaveliseev/@dev/SfSpp_preprocessing/data_path.h`, the `#ifdef`
chain currently ends in `static const std::string potID = "A";` and friends.
Replace each arm's assignment with a `#define SFSPP_DEFAULT_POT_ID "X"`, then
after the chain:

```cpp
// The pot id was a compile-time constant, so one binary served one pot and a
// sweep would have needed ten builds. SFSPP_POT_ID overrides it at run time;
// the compiled value stays as the default.
inline std::string sfsppPotID() {
    const char* v = std::getenv("SFSPP_POT_ID");
    return (v && *v) ? std::string(v) : std::string(SFSPP_DEFAULT_POT_ID);
}
static const std::string potID = sfsppPotID();
```

- [ ] **Step 7: Run tests and build**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
.venv/bin/pytest sfsbench/tests/test_preset.py -v
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp
PYTHONPATH=tools /opt/homebrew/bin/python3 -c "
from sfsbench.preset import discover, generate
R='/Users/vaceslaveliseev/@dev/SfSpp_preprocessing/DatasetRef'
generate(discover(R, ['A']), 'class/data_path_generated.h',
         gt_root='/Users/vaceslaveliseev/@dev/structure-from-sherds-pp/Dataset/SfS_pp/Ground Truth')
"
head -20 class/data_path_generated.h
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release > /dev/null && cmake --build build -j10 --target SfSpp_gen 2>&1 | grep -E "error:|Built target"

cd /Users/vaceslaveliseev/@dev/SfSpp_preprocessing
cmake --build build -j10 2>&1 | grep -E "error:|Built target"
SFSPP_POT_ID=B SFSPP_DATASET_ROOT=/tmp/nonexistent ./build/MeshPreprocessing 2>&1 | head -3
```

Expected: 7 passed; `#define SHARD_NUMBER 8` in the generated header;
`Built target SfSpp_gen`; both preprocessing targets built; and the last command
mentioning `Pot_B` in its path error, proving the override reaches the paths.

- [ ] **Step 8: Commit**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp
git add tools/sfsbench/preset.py tools/sfsbench/tests/test_preset.py \
        class/data_path.h CMakeLists.txt .gitignore
git commit -m "Generate the reassembly preset from the staged collection

class/data_path.h enumerated every fragment path by hand in one of 30 mutually
exclusive ifdef blocks, so an n-fragment collection meant writing ~9n lines and
recompiling. Only main.cpp includes it and only main.cpp uses SHARD_NUMBER, so a
generated header costs one translation unit (~5 s) per collection."

cd /Users/vaceslaveliseev/@dev/SfSpp_preprocessing
git add data_path.h
git commit -m "Take the pot id from SFSPP_POT_ID at run time

potID was a compile-time constant, so one preprocessing binary served one pot
and a ten-pot sweep would have needed ten builds."
```

---

### Task 2: Wire the published metric into the binary

*(New. The whole benchmark depends on it — without this the harness has no access
to the number the papers publish.)*

`CountResult` already computes Sherd Accuracy exactly as SfS++ defines it, and
`SaveAcc` writes `Result/1. Acc.txt`. But it is reached only from the `vis.fine_`
branch (`main.cpp:388`), i.e. the interactive `f` key, while `SFS_AUTOSAVE=exit`
goes `vis.first_` → `vis.save_` → `SaveResult` → `break` and never scores.

Three facts that make the change safe, all verified in the source:

1. The `fine_` branch does **not** re-optimise despite its log line. It restores
   `T_result` from the stored per-state graph transforms — `T_result[j] = T_axis[j]`
   composed with `graph_[i].T_[j]` — and the `R_fine`/`t_fine` composition is
   identity throughout. The `t_dummy[1] += 200*i` offsets are applied to the
   *viewer* and to `manager.shard_`, not to `T_result`.
2. `CountResult`'s `graph` parameter is **unused** — the only line that read it is
   commented out. This matters because `main.cpp:418` declares
   `MatrixXd graph_dummy(SHARD_NUMBER, SHARD_NUMBER)` and never assigns it, then
   reads it at :433; the resulting `graph` is uninitialised memory. Harmless
   today, but do not "fix" it by feeding `graph` into the metric.
3. `CountResult` **mutates** `right_sherd` (it rewrites the vector with the
   per-sherd verdict), so it must be handed a fresh copy initialised from
   `shard_on_off`, not the live one the viewer uses.

**Files:**
- Modify: `main.cpp` — extract the reconstruction, score on the autosave path
- Modify: `class/data_path.h` — `shard_on_off` for D/F/I, `POT_J`/`POT_C_J` counts, `POT_CTRL` ground-truth paths

- [ ] **Step 1: Extract the transform reconstruction**

Add above `main`, and call it from the existing `vis.fine_` branch in place of
the inline restore so there is exactly one implementation:

```cpp
// The result transforms are reconstructable from the stored per-state graph
// without touching the viewer: T_result[j] = T_axis[j] composed with the graph
// transform for j. The fine_ branch also displaces meshes by 200*i on y for
// display; that offset is deliberately not applied here.
static vector<Trans> ReconstructResultTransforms(StateManager& manager,
                                                 int state_index,
                                                 const vector<Trans>& T_axis)
{
	vector<Trans> T_result = T_axis;
	const int num_graph = manager.out_state_[state_index].graph_.size();
	for (int i = 0; i < num_graph; i++) {
		for (int j = 0; j < SHARD_NUMBER; j++) {
			if (manager.out_state_[state_index].graph_[i].node_[j]) {
				Matrix3d R_prev;
				Vector3d t_prev;
				manager.out_state_[state_index].graph_[i].T_[j].Output(R_prev, t_prev);
				T_result[j].Input(R_prev, t_prev);
			}
		}
	}
	return T_result;
}
```

- [ ] **Step 2: Score on the autosave path**

In the `else if (vis.save_)` branch, after `SaveResult` and before the
`autosave_exit` break:

```cpp
if (autosave) {
	// CountResult is the metric the papers publish (relative pose, 20 deg /
	// 50 mm, a sherd correct with at least one correct edge). It was reachable
	// only through the interactive 'f' key, so an autosave run produced result
	// meshes and no score. It mutates right_sherd, hence the copy.
	vector<Trans> T_scored = ReconstructResultTransforms(manager, count_move_state, T_axis);
	vector<bool> right_scored(shard_on_off, shard_on_off + SHARD_NUMBER);
	MatrixXd graph_unused(SHARD_NUMBER, SHARD_NUMBER);   // CountResult ignores it
	graph_unused.setZero();
	int k_sherd, t_sherd, k_edge, t_edge;
	tie(k_sherd, t_sherd, k_edge, t_edge)
		= CountResult(GT_graph, GT_trans, graph_unused, T_scored, right_scored);
	pair<int, int> sherd_acc = make_pair(k_sherd, t_sherd);
	pair<int, int> edge_acc = make_pair(k_edge, t_edge);
	SaveAcc(path_result, sherd_acc, edge_acc, time_result);
}
```

Scoring goes **after** `SaveResult` so the saved meshes are unaffected, and it
touches neither `pc_origin` nor `manager.shard_`.

- [ ] **Step 3: Disable the fragments the paper excludes**

In `class/data_path.h`, set `shard_on_off` to `false` for:

| preset | index (1-based) | file |
|---|---|---|
| `POT_D` | 22 | `Pot_D_Piece_22_Breakline_0.pcd` |
| `POT_F` | 7 | `Pot_F_Piece_07_Breakline_0.pcd` |
| `POT_I` | 28, 29, 30 | `Pot_I_Piece_{28,29,30}_Breakline_0.pcd` |
| `POT_All` | the same five, at their offsets | |

`POT_J` already has index 9 disabled. Without this, `CountResult`'s `total_sherd`
counts 29 / 7 / 30 instead of 28 / 6 / 27 and no run can ever match the published
percentage. (The generated preset from Task 1 derives the same set from the
ground-truth graph automatically; this step keeps the shipped presets from being
landmines.)

- [ ] **Step 4: Fix two upstream defects in the same file**

- `POT_J`: `SHARD_NUMBER 19` but only **12** entries in `file_path` and
  `shard_on_off`; the tail is value-initialised to empty strings and `false`, and
  `GT_graph` is sized 19×19 against a 12×12 graph file. Correct it to `12`.
- `POT_C_J`: the same defect, `SHARD_NUMBER 23` against **16** entries. Correct
  it to `16`.
- `POT_CTRL` (and `POT_TEST`, same block shape): its ground-truth paths use
  `GroundTruth/Transformation/...` and `GroundTruth/Pot_A_simple_graph.txt`,
  neither of which resolves — `Dataset/SfS_pp` spells the folder `Ground Truth/`
  with a space and has no `Transformation/` subdirectory, and in the ICCV layout
  the graph lives *inside* `Transformation/`. `main.cpp` reads the graph with a
  bare `ifstream` and no existence check, so today the comment claiming this run
  is scored against real ground truth is false and the accuracy report is
  silently disabled. Point both at the `Ground Truth/` layout:

```cpp
string gt_T_path[SHARD_NUMBER] = {
	path + "Ground Truth/Pot_A_Piece_1_T.txt", ...
};
string gt_graph_path[1] = { path + "Ground Truth/Pot_A_simple_graph.txt" };
```

  Note this makes the paths relative to `SFS_DATA_ROOT`, so a `POT_CTRL` run must
  either stage the ground truth into its working root or symlink it there.

- [ ] **Step 5: Build and verify the score appears without a keypress**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp
cmake --build build -j10 2>&1 | grep -E "error:|Built target"
R="/Users/vaceslaveliseev/@dev/structure-from-sherds-pp/Dataset/SfS_pp"
REF=/Users/vaceslaveliseev/@dev/SfSpp_preprocessing/DatasetRef
ln -sfn "$R/Ground Truth" "$REF/Ground Truth"        # if not already staged
rm -rf "$REF/Result" && mkdir -p "$REF/Result"
PYTHONPATH=tools /opt/homebrew/bin/python3 -c "
from sfsbench.preset import discover, generate
generate(discover('$REF', ['A']), 'class/data_path_generated.h', gt_root='$R/Ground Truth')
"
cmake --build build -j10 --target SfSpp_gen 2>&1 | grep -E "error:|Built target"
SFS_DATA_ROOT="$REF/" SFS_AUTOSAVE=exit ./build/SfSpp_gen 2>&1 | grep -E "Sherd accuracy|Edge accuracy"
cat "$REF/Result/1. Acc.txt"
```

Expected: `********* Sherd accuracy : 8 / 8 (100%)` and a written `1. Acc.txt`.

**This is the gate for the whole plan.** If the authors' preprocessed Pot A does
not score 8/8 through the non-interactive path, stop: either the reconstruction
is wrong or the ground truth is not being loaded, and every later number would be
meaningless.

- [ ] **Step 6: Commit**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp
git add main.cpp class/data_path.h
git commit -m "Score the result on the non-interactive path

CountResult already implements the published Sherd Accuracy (relative pose,
20 deg / 50 mm, a sherd correct with at least one correct edge) but was reachable
only through the interactive 'f' key, so SFS_AUTOSAVE runs produced meshes and no
number. Also disable the fragments the paper excludes from evaluation (D 22, F 7,
I 28-30), fix POT_J/POT_C_J declaring more shards than they define, and repair
POT_CTRL's ground-truth paths, which resolved to nothing."
```

---

### Task 3: Harness — dataset staging

*(Was Task 1. Unchanged except where noted.)*

**Files:**
- Create: `tools/sfsbench/__init__.py`
- Create: `tools/sfsbench/stage.py`
- Create: `tools/sfsbench/tests/test_stage.py`

**Interfaces:**
- Produces:
  - `read_mesh(path: str) -> tuple[np.ndarray, np.ndarray]` — returns `(vertices Nx3 float64, faces Mx3 int64)`. Accepts `.ply` (binary little-endian, `float x,y,z` + `uchar red,green,blue`, `uchar`-count int faces) and `.obj`.
  - `vertex_normals(xyz: np.ndarray, tri: np.ndarray) -> np.ndarray` — area-weighted, unit length, Nx3.
  - `voxel_subsample(xyz: np.ndarray, target: int) -> np.ndarray` — returns **indices** of one representative vertex per voxel, voxel size chosen so the count lands within [0.75, 1.3] × target.
  - `write_pcd(path: str, xyz: np.ndarray, nrm: np.ndarray, binary: bool = False) -> None` — PCD, fields `x y z normal_x normal_y normal_z curvature`.
  - `write_obj(path: str, xyz: np.ndarray, tri: np.ndarray) -> None`
  - `stage_dataset(src_dir: str, dst_root: str, pot_id: str, target_points: int = 2_000_000) -> list[str]` — converts every mesh in `src_dir` into `<dst_root>/Point/Pot_<pot_id>/Pot_<pot_id>_Piece_NN_Point.pcd` and `<dst_root>/Mesh/Pot_<pot_id>/Pot_<pot_id>_Piece_NN_Mesh.obj`, returns the fragment names in order.

**Why the naming matters:** `mesh_processing.cpp:338` does
`filePath.stem().string().substr(0, 14)`. `Pot_A_Piece_01_Point` truncated to 14
characters is exactly `Pot_A_Piece_01`. It does *not* throw on shorter stems, so
the real constraint is "≤ 14 characters, or unique in the first 14" — the failure
mode for longer stems is a silent collision. Note two other naming schemes
coexist: `substr(0, find("_Point"))` at `mesh_processing.cpp:1720` and
`substr(0, find("_Mesh"))` at `edgeline_extraction.cpp:2464, 2894` — the latter
builds the axis filename, so a rename breaks the schemes inconsistently.

**Always write `_Mesh.obj`**, even when the source was `_Mesh_DS.obj`:
`mesh_processing.cpp:1722` composes the mesh filename as `<stem>_Mesh.obj`
unconditionally.

- [ ] **Step 1: Write the failing test**

```python
# tools/sfsbench/tests/test_stage.py
import os, tempfile
import numpy as np
from sfsbench.stage import (read_mesh, vertex_normals, voxel_subsample,
                            write_pcd, write_obj, stage_dataset)

def _unit_square_obj(path):
    """Two triangles in the z=0 plane; every vertex normal must be +z or -z."""
    with open(path, "w") as f:
        f.write("v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\n")
        f.write("f 1 2 3\nf 1 3 4\n")

def test_read_obj_roundtrip():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "a.obj")
        _unit_square_obj(p)
        V, F = read_mesh(p)
        assert V.shape == (4, 3)
        assert F.shape == (2, 3)
        assert F.max() == 3          # zero-based

def test_vertex_normals_are_unit_and_axis_aligned():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "a.obj")
        _unit_square_obj(p)
        V, F = read_mesh(p)
        N = vertex_normals(V, F)
        assert np.allclose(np.linalg.norm(N, axis=1), 1.0)
        assert np.allclose(np.abs(N[:, 2]), 1.0)

def test_voxel_subsample_hits_target_band():
    rng = np.random.default_rng(0)
    xyz = rng.random((200_000, 3)) * 100.0
    idx = voxel_subsample(xyz, 20_000)
    assert 15_000 <= len(idx) <= 26_000
    assert len(np.unique(idx)) == len(idx)

def test_write_pcd_header_is_readable_by_pcl_layout():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "a.pcd")
        xyz = np.zeros((3, 3)); nrm = np.tile([0.0, 0.0, 1.0], (3, 1))
        write_pcd(p, xyz, nrm)
        lines = open(p).read().splitlines()
        assert lines[2] == "FIELDS x y z normal_x normal_y normal_z curvature"
        assert lines[10] == "DATA ascii"
        assert len(lines) == 11 + 3

def test_stage_dataset_uses_14_char_names():
    with tempfile.TemporaryDirectory() as d:
        src = os.path.join(d, "src"); os.makedirs(src)
        _unit_square_obj(os.path.join(src, "whatever.obj"))
        names = stage_dataset(src, os.path.join(d, "ds"), "A", target_points=4)
        assert names == ["Pot_A_Piece_01"]
        assert len(names[0]) == 14
        assert os.path.exists(os.path.join(d, "ds", "Point", "Pot_A", "Pot_A_Piece_01_Point.pcd"))
        assert os.path.exists(os.path.join(d, "ds", "Mesh", "Pot_A", "Pot_A_Piece_01_Mesh.obj"))

def test_stage_dataset_normalises_the_downsampled_suffix():
    with tempfile.TemporaryDirectory() as d:
        src = os.path.join(d, "src"); os.makedirs(src)
        _unit_square_obj(os.path.join(src, "Pot_C_Piece_01_Mesh_DS.obj"))
        stage_dataset(src, os.path.join(d, "ds"), "C", target_points=4)
        assert os.path.exists(os.path.join(d, "ds", "Mesh", "Pot_C", "Pot_C_Piece_01_Mesh.obj"))
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
.venv/bin/pytest sfsbench/tests/test_stage.py -v
```

Expected: FAIL — `ModuleNotFoundError: No module named 'sfsbench.stage'`

- [ ] **Step 3: Write the implementation**

```python
# tools/sfsbench/stage.py
"""Convert raw scanned meshes into the layout the preprocessing binaries expect."""
import glob, os
import numpy as np


def _read_ply(path):
    with open(path, "rb") as f:
        hdr = b""
        while b"end_header" not in hdr:
            chunk = f.readline()
            if not chunk:
                raise ValueError(f"{path}: no end_header")
            hdr += chunk
        nv = nf = 0
        for line in hdr.split(b"\n"):
            if line.startswith(b"element vertex"):
                nv = int(line.split()[-1])
            elif line.startswith(b"element face"):
                nf = int(line.split()[-1])
        vdt = np.dtype([("x", "<f4"), ("y", "<f4"), ("z", "<f4"),
                        ("r", "u1"), ("g", "u1"), ("b", "u1")])
        V = np.frombuffer(f.read(nv * vdt.itemsize), dtype=vdt, count=nv)
        fdt = np.dtype([("n", "u1"), ("a", "<i4"), ("b", "<i4"), ("c", "<i4")])
        F = np.frombuffer(f.read(nf * fdt.itemsize), dtype=fdt, count=nf)
    if nf and not (F["n"] == 3).all():
        raise ValueError(f"{path}: non-triangular faces")
    xyz = np.stack([V["x"], V["y"], V["z"]], 1).astype(np.float64)
    tri = np.stack([F["a"], F["b"], F["c"]], 1).astype(np.int64)
    return xyz, tri


def _read_obj(path):
    V, F = [], []
    with open(path, "r", errors="replace") as f:
        for line in f:
            if line.startswith("v "):
                t = line.split()
                V.append((float(t[1]), float(t[2]), float(t[3])))
            elif line.startswith("f "):
                F.append([int(tok.split("/")[0]) - 1 for tok in line.split()[1:4]])
    return np.asarray(V, dtype=np.float64), np.asarray(F, dtype=np.int64)


def read_mesh(path):
    ext = os.path.splitext(path)[1].lower()
    if ext == ".ply":
        return _read_ply(path)
    if ext == ".obj":
        return _read_obj(path)
    raise ValueError(f"unsupported mesh format: {path}")


def vertex_normals(xyz, tri):
    fn = np.cross(xyz[tri[:, 1]] - xyz[tri[:, 0]], xyz[tri[:, 2]] - xyz[tri[:, 0]])
    vn = np.zeros_like(xyz)
    for k in range(3):
        np.add.at(vn, tri[:, k], fn)          # unnormalised => area weighting
    n = np.linalg.norm(vn, axis=1, keepdims=True)
    n[n == 0] = 1.0
    return vn / n


def voxel_subsample(xyz, target):
    if len(xyz) <= target:
        return np.arange(len(xyz))
    origin = xyz.min(0)
    extent = xyz.max(0) - origin
    v = (float(np.prod(extent)) / target) ** (1.0 / 3.0) * 0.35
    idx = np.arange(len(xyz))
    for _ in range(12):
        k = np.floor((xyz - origin) / v).astype(np.int64)
        dims = k.max(0) + 1
        h = k[:, 0] + dims[0] * (k[:, 1] + dims[1] * k[:, 2])
        _, idx = np.unique(h, return_index=True)
        n = len(idx)
        if 0.75 * target <= n <= 1.3 * target:
            break
        v *= (n / target) ** (1.0 / 3.0)
    return np.sort(idx)


def write_pcd(path, xyz, nrm, binary=False):
    """ASCII by default because that is what the reference files use.

    Measured: 1.75 M points of our own scans came to 97 MB per fragment, 390 MB
    for four, written line by line by np.savetxt. At the 2 M default a
    100-fragment collection is roughly 11 GB of scratch per run. Binary PCD is
    ~7x smaller and far faster to write, and PCL reads it; switch when the size
    becomes the bottleneck (see Task 12, Step 1).
    """
    n = len(xyz)
    data = np.hstack([xyz, nrm, np.zeros((n, 1))]).astype(np.float32)
    header = ("# .PCD v0.7 - Point Cloud Data file format\n"
              "VERSION 0.7\n"
              "FIELDS x y z normal_x normal_y normal_z curvature\n"
              "SIZE 4 4 4 4 4 4 4\n"
              "TYPE F F F F F F F\n"
              "COUNT 1 1 1 1 1 1 1\n"
              f"WIDTH {n}\n"
              "HEIGHT 1\n"
              "VIEWPOINT 0 0 0 1 0 0 0\n"
              f"POINTS {n}\n")
    if binary:
        with open(path, "wb") as f:
            f.write((header + "DATA binary\n").encode())
            f.write(data.tobytes())
        return
    with open(path, "w") as f:
        f.write(header + "DATA ascii\n")
        np.savetxt(f, data, fmt="%.6g")


def write_obj(path, xyz, tri):
    """Note: no `vt` lines.

    That is deliberate and load-bearing. VTK duplicates points when loading an
    OBJ that carries per-corner texture indices, which is why results from the
    authors' meshes come back with 3 x n_faces vertices while results from ours
    come back with n_vertices. sfsbench.score picks the layout from the vertex
    count; do not "fix" this by emitting vt without updating that.
    """
    with open(path, "w") as f:
        f.write("# staged by sfsbench\n")
        np.savetxt(f, xyz, fmt="v %.6g %.6g %.6g")
        np.savetxt(f, tri + 1, fmt="f %d %d %d")


def stage_dataset(src_dir, dst_root, pot_id, target_points=2_000_000):
    """Raw meshes -> <dst_root>/{Point,Mesh}/Pot_<id>/Pot_<id>_Piece_NN_*.

    The 14-character name is mandatory: mesh_processing.cpp truncates the file
    stem to 14 characters to recover the fragment id. The mesh is always written
    as _Mesh.obj even when the source was _Mesh_DS.obj, because
    mesh_processing.cpp composes that name unconditionally.
    """
    point_dir = os.path.join(dst_root, "Point", f"Pot_{pot_id}")
    mesh_dir = os.path.join(dst_root, "Mesh", f"Pot_{pot_id}")
    os.makedirs(point_dir, exist_ok=True)
    os.makedirs(mesh_dir, exist_ok=True)

    sources = sorted(glob.glob(os.path.join(src_dir, "*.ply")) +
                     glob.glob(os.path.join(src_dir, "*.obj")))
    names = []
    for i, src in enumerate(sources, start=1):
        name = f"Pot_{pot_id}_Piece_{i:02d}"
        assert len(name) == 14, f"name must be 14 chars, got {name!r}"
        xyz, tri = read_mesh(src)
        nrm = vertex_normals(xyz, tri)
        keep = voxel_subsample(xyz, target_points)
        write_pcd(os.path.join(point_dir, f"{name}_Point.pcd"), xyz[keep], nrm[keep])
        write_obj(os.path.join(mesh_dir, f"{name}_Mesh.obj"), xyz, tri)
        names.append(name)
    return names
```

Also create the package marker:

```python
# tools/sfsbench/__init__.py
"""Harness for the Structure-from-Sherds pipeline."""
```

- [ ] **Step 4: Run test to verify it passes**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
.venv/bin/pytest sfsbench/tests/test_stage.py -v
```

Expected: PASS, 6 passed

- [ ] **Step 5: Commit**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp
git add tools/sfsbench/__init__.py tools/sfsbench/stage.py tools/sfsbench/tests/test_stage.py
git commit -m "harness: stage raw meshes into the preprocessing layout"
```

---

### Task 4: The millimetre diagnostic (Kabsch)

*(Was Task 2. Demoted from "the gate" to "a diagnostic", gate value corrected,
face expansion made conditional.)*

This measures **absolute displacement in millimetres after one global rigid
alignment** — how far off a fragment is, which the published Sherd Accuracy does
not tell you. It is not comparable to SA in either direction and must never be
reported as the published number.

**Files:**
- Create: `tools/sfsbench/score.py`
- Create: `tools/sfsbench/tests/test_score.py`

**Interfaces:**
- Consumes: `read_mesh` from `sfsbench.stage`
- Produces:
  - `kabsch(P, Q) -> tuple[np.ndarray, np.ndarray]` — rigid transform `(R, t)` mapping `P` onto `Q`, reflection-free.
  - `source_points(Vr: np.ndarray, Vo: np.ndarray, Fo: np.ndarray) -> np.ndarray` — picks the source layout that matches the result's vertex count, raising if neither does.
  - `score_result(result_dir, source_mesh_dir, gt_dir, pot_id, stride=37) -> dict` — returns `{"placed": int, "total_scored": int, "within_20mm": int, "per_fragment": {id: mean_error_mm}, "mean_error": float, "mean_of_means": float}`.
  - `read_acc(result_dir: str) -> dict | None` — parses `Result/1. Acc.txt` into `{"sherd": (k, n), "edge": (k, n)}`; `None` when the file is absent (no ground truth, or the run predates Task 2).

**Why the layout must be decided, not assumed:** see Global Constraints. The
previous revision expanded the source unconditionally and then hid any mismatch
behind `n = min(len(src), len(res))`. That silently destroys correspondence for
every run the harness itself stages — the exact failure it warned about, with the
sign flipped. Raise instead of truncating.

- [ ] **Step 1: Write the failing test**

```python
# tools/sfsbench/tests/test_score.py
import numpy as np
import pytest
from sfsbench.score import kabsch, source_points

def test_kabsch_recovers_a_known_transform():
    rng = np.random.default_rng(1)
    P = rng.random((50, 3)) * 10
    theta = 0.7
    Rt = np.array([[np.cos(theta), -np.sin(theta), 0],
                   [np.sin(theta),  np.cos(theta), 0],
                   [0, 0, 1]])
    tt = np.array([5.0, -2.0, 1.0])
    Q = (Rt @ P.T).T + tt
    R, t = kabsch(P, Q)
    assert np.allclose(R, Rt, atol=1e-9)
    assert np.allclose(t, tt, atol=1e-9)

def test_kabsch_never_reflects():
    rng = np.random.default_rng(2)
    P = rng.random((30, 3))
    Q = P.copy(); Q[:, 0] *= -1          # a reflection, not a rotation
    R, _ = kabsch(P, Q)
    assert np.linalg.det(R) > 0

V = np.array([[0., 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]])
F = np.array([[0, 1, 2], [0, 2, 3]])

def test_source_points_expands_when_the_result_is_face_expanded():
    Vr = np.zeros((6, 3))                 # 3 x n_faces: the authors' meshes
    S = source_points(Vr, V, F)
    assert S.shape == (6, 3)
    assert np.allclose(S[0], V[0])
    assert np.allclose(S[3], V[0])
    assert np.allclose(S[4], V[2])

def test_source_points_stays_raw_when_the_result_is_not_expanded():
    Vr = np.zeros((4, 3))                 # n_vertices: meshes our stager writes
    S = source_points(Vr, V, F)
    assert S.shape == (4, 3)
    assert np.allclose(S, V)

def test_source_points_raises_on_a_mismatch():
    with pytest.raises(ValueError):
        source_points(np.zeros((5, 3)), V, F)
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
.venv/bin/pytest sfsbench/tests/test_score.py -v
```

Expected: FAIL — `ModuleNotFoundError: No module named 'sfsbench.score'`

- [ ] **Step 3: Write the implementation**

```python
# tools/sfsbench/score.py
"""Millimetre displacement of a reassembly against ground truth.

This is a DIAGNOSTIC, not the published metric. The gate is Sherd Accuracy, which
the binary computes itself (class/data_structure.cpp:920) and writes to
Result/1. Acc.txt; read_acc() parses it. The two are not comparable: SA is
per-edge and relative and counts a sherd correct with one correct edge, while
this pools vertices, fits one global transform and reports millimetres.
"""
import glob, os, re
import numpy as np

from sfsbench.stage import read_mesh

DISPLACEMENT_MM = 20.0        # a diagnostic band, not a published threshold


def kabsch(P, Q):
    """Rigid transform (R, t) with Q ~= R @ P + t. Reflection-free."""
    cp, cq = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - cp).T @ (Q - cq))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, cq - R @ cp


def source_points(Vr, Vo, Fo):
    """Match the source layout to the result's, from the data.

    VTK duplicates points on load only when the input OBJ carries per-corner vt
    indices. The authors' meshes do; meshes written by sfsbench.stage do not. So
    the layout is a property of the input file, never a rule.
    """
    if len(Vr) == 3 * len(Fo):
        return Vo[Fo.reshape(-1)]
    if len(Vr) == len(Vo):
        return Vo
    raise ValueError(f"result has {len(Vr)} vertices, which matches neither "
                     f"{len(Vo)} source vertices nor {3 * len(Fo)} face corners")


def read_acc(result_dir):
    """Parse the Sherd/Edge accuracy the binary wrote. None when absent."""
    path = os.path.join(result_dir, "1. Acc.txt")
    if not os.path.exists(path):
        return None
    text = open(path, errors="replace").read()
    nums = re.findall(r"(\d+)\s*/\s*(\d+)", text)
    if len(nums) < 2:
        return None
    return {"sherd": (int(nums[0][0]), int(nums[0][1])),
            "edge": (int(nums[1][0]), int(nums[1][1]))}


def _source_mesh(source_mesh_dir, pot_id, i):
    for suffix in ("_Mesh.obj", "_Mesh_DS.obj"):
        p = os.path.join(source_mesh_dir, f"Pot_{pot_id}_Piece_{i:02d}{suffix}")
        if os.path.exists(p):
            return p
    raise FileNotFoundError(
        f"no source mesh for Pot_{pot_id}_Piece_{i:02d} in {source_mesh_dir}")


def score_result(result_dir, source_mesh_dir, gt_dir, pot_id, stride=37):
    files = sorted(glob.glob(os.path.join(result_dir, "*_Top_1_OBJ_*.obj")),
                   key=lambda p: int(re.search(r"OBJ_(\d+)\.obj$", p).group(1)))
    ours, gts, ids = [], [], []
    for path in files:
        i = int(re.search(r"OBJ_(\d+)\.obj$", path).group(1))
        gt_file = os.path.join(gt_dir, f"Pot_{pot_id}_Piece_{i}_T.txt")
        if not os.path.exists(gt_file):
            continue                     # excluded fragment, or no ground truth
        Vo, Fo = read_mesh(_source_mesh(source_mesh_dir, pot_id, i))
        Vr, _ = read_mesh(path)
        src = source_points(Vr, Vo, Fo)[::stride]
        res = Vr[::stride]
        if len(src) != len(res):
            raise ValueError(f"{path}: {len(res)} result points vs {len(src)} source")
        T = np.loadtxt(gt_file)
        gts.append((T[:3, :3] @ src.T).T + T[:3, 3])
        ours.append(res)
        ids.append(i)

    out = {"placed": len(files), "total_scored": len(ids), "within_20mm": 0,
           "per_fragment": {}, "mean_error": float("nan"),
           "mean_of_means": float("nan")}
    if not ids:
        return out

    R, t = kabsch(np.vstack(ours), np.vstack(gts))
    all_err = []
    for i, O, G in zip(ids, ours, gts):
        e = np.linalg.norm(((R @ O.T).T + t) - G, axis=1)
        out["per_fragment"][i] = float(e.mean())
        out["within_20mm"] += int(e.mean() < DISPLACEMENT_MM)
        all_err.append(e)
    out["mean_error"] = float(np.concatenate(all_err).mean())        # pooled
    out["mean_of_means"] = float(np.mean(list(out["per_fragment"].values())))
    return out
```

- [ ] **Step 4: Run test to verify it passes**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
.venv/bin/pytest sfsbench/tests/test_score.py -v
```

Expected: PASS, 5 passed

- [ ] **Step 5: Verify against the two runs already on disk**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
R="/Users/vaceslaveliseev/@dev/structure-from-sherds-pp/Dataset/SfS_pp"
PYTHONPATH=. .venv/bin/python -c "
from sfsbench.score import score_result
R='$R'
for tag, d, src in [
    ('DatasetRef', '/Users/vaceslaveliseev/@dev/SfSpp_preprocessing/DatasetRef/Result', R+'/Mesh'),
    ('DatasetA',   '/Users/vaceslaveliseev/@dev/SfSpp_preprocessing/DatasetA/Result',   R+'/Mesh'),
]:
    r = score_result(d, src, R+'/Ground Truth', 'A')
    print(tag, r['within_20mm'], '/', r['total_scored'],
          'pooled', round(r['mean_error'], 2), 'mean-of-means', round(r['mean_of_means'], 2))
"
```

Expected, both reproduced during review from the artefacts on disk:

```
DatasetRef 8 / 8 pooled 3.74 mean-of-means 4.37
DatasetA   0 / 8 pooled 89.46 mean-of-means ~90
```

*(The previous revision expected `mean 4.4`. A correct implementation produces
3.739 pooled; 4.37 is the unweighted mean of the eight per-fragment means, since
the large, well-placed piece 1 dominates the pool. Both are now reported, so the
ambiguity cannot recur.)*

If **8/8 and 0/8** do not reproduce, stop — the diagnostic is wrong. The ±0.05
tolerance applies to the millimetre figures.

- [ ] **Step 6: Commit**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp
git add tools/sfsbench/score.py tools/sfsbench/tests/test_score.py
git commit -m "harness: millimetre displacement diagnostic, plus 1. Acc.txt reader"
```

---

### Task 5: Pipeline runner and the reference report

*(Was Task 3. Result/ is now cleared per run, stats are keyed by fragment, the
score regex is anchored, and the runner regenerates the preset and rebuilds.)*

**Files:**
- Create: `tools/sfsbench/run.py`
- Create: `tools/sfsbench/cli.py`
- Create: `tools/sfsbench/tests/test_run.py`

**Interfaces:**
- Consumes: `stage_dataset` (Task 3), `score_result` / `read_acc` (Task 4), `discover` / `generate` (Task 1)
- Produces:
  - `PipelinePaths` — dataclass `pp_repo, prep_repo, dataset_root, temp_root, pot_id, gt_root=None`
  - `run_pipeline(paths, env=None, binary="SfSpp_gen") -> dict` — regenerates the preset, rebuilds that target, then runs MeshPreprocessing, `extract_axis` per fragment, EdgeLineExtraction and the reassembly. Returns
    `{"clusters": {name: int}, "breaklines": {name: {"segments": int, "points": int, "rim": int, "base": int, "info": int}}, "axes": {name: int}, "matches_raw": int, "matches_pruned": int, "score": int, "result_dir": str, "placed": int}`
  - `outcome(stats: dict, acc: dict | None) -> str` — one of `"assembled"`, `"misassembled"`, `"no-joins"`.

**Outcome rules** (the design requires three, not two):
- `"no-joins"` when `matches_pruned == 0` or fewer than two fragments appear in the result
- `"assembled"` when accuracy exists and `acc["sherd"][0] == acc["sherd"][1]`
- `"misassembled"` otherwise

Note `"assembled"` is now decided by **Sherd Accuracy**, not by the millimetre
diagnostic.

Four corrections carried into this task:

- **`Result/` must be cleared at the start of every run.** Filenames are
  `<Y>_<M>_<D>_<H>_<Min>_Top_1_OBJ_<i>.obj`, so a second run in a different
  minute *adds* files instead of replacing them: `placed` inflates and
  `score_result`'s glob returns two files per index whose sort key cannot tell
  them apart, scoring a mixture of two runs.
- **`axes[name]` is a candidate count, not a fixed number.** `extract_axis`
  writes 1–10 lines, one per surviving candidate after 10° de-duplication and the
  `SFS_AXIS_COST_TOLERANCE` filter; real files in `DatasetA/Axes/` have 1, 1, 1,
  1, 1, 1, 2, 3. More than one is a signal worth reporting, not an error.
- **Anchor the `Score :` regex.** Two producers exist — `visualize.cpp:358`
  (`"Result : r/n, Score : ..."`) and `visualize.cpp:414` — and an unanchored
  `re.search` matches whichever comes first in the concatenated log.
- **Key the counts by fragment.** The previous revision declared dicts and
  returned bare lists, so a reader could not tell which fragment a count belonged
  to. The `[breakline]` log line carries the full path, so those parse by name
  directly; cluster counts do not, but `mesh_processing.cpp:1712` sorts the input
  files, so zipping with the sorted name list is sound — and asserted.

- [ ] **Step 1: Write the failing test**

```python
# tools/sfsbench/tests/test_run.py
from sfsbench.run import outcome, parse_stats

FULL = {"sherd": (8, 8), "edge": (15, 15)}
PART = {"sherd": (3, 8), "edge": (4, 15)}

def test_no_joins_when_nothing_survives_pruning():
    assert outcome({"matches_pruned": 0, "placed": 4}, None) == "no-joins"

def test_no_joins_when_fewer_than_two_placed():
    assert outcome({"matches_pruned": 12, "placed": 1}, None) == "no-joins"

def test_assembled_when_sherd_accuracy_is_complete():
    assert outcome({"matches_pruned": 46, "placed": 8}, FULL) == "assembled"

def test_misassembled_when_sherd_accuracy_is_partial():
    assert outcome({"matches_pruned": 46, "placed": 8}, PART) == "misassembled"

def test_misassembled_without_ground_truth():
    assert outcome({"matches_pruned": 46, "placed": 4}, None) == "misassembled"

def test_parse_stats_reads_the_binaries_stdout():
    log = "\n".join([
        "Number of clusters is equal to 12",
        "Number of clusters is equal to 9",
        "[breakline] /ds/Breaklines/Pot_A/Pot_A_Piece_01_Breakline_0.pcd: 5 segment(s), 209 points, rim 1, base 2 => info 2",
        "[breakline] /ds/Breaklines/Pot_A/Pot_A_Piece_02_Breakline_0.pcd: 3 segment(s), 173 points, rim 0, base 0 => info 0",
        "Total number : 3587",
        "Total number pruned : 459",
        "Result : 1/8, Score : 999",
        "Score : 205",
    ])
    s = parse_stats(log, ["Pot_A_Piece_01", "Pot_A_Piece_02"])
    assert s["matches_raw"] == 3587
    assert s["matches_pruned"] == 459
    assert s["score"] == 205                       # line-anchored, not "Result : ..., Score : 999"
    assert s["clusters"] == {"Pot_A_Piece_01": 12, "Pot_A_Piece_02": 9}
    assert s["breaklines"]["Pot_A_Piece_01"]["segments"] == 5
    assert s["breaklines"]["Pot_A_Piece_01"]["info"] == 2
    assert s["breaklines"]["Pot_A_Piece_02"]["points"] == 173

def test_parse_stats_refuses_to_guess_when_the_counts_disagree():
    log = "Number of clusters is equal to 12"
    s = parse_stats(log, ["Pot_A_Piece_01", "Pot_A_Piece_02"])
    assert s["clusters"] == {}                     # rather than mislabel one
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
.venv/bin/pytest sfsbench/tests/test_run.py -v
```

Expected: FAIL — `ModuleNotFoundError: No module named 'sfsbench.run'`

- [ ] **Step 3: Write the implementation**

```python
# tools/sfsbench/run.py
"""Drive the C++ binaries end to end and collect the numbers diagnosis needs."""
import os, re, shutil, subprocess
from dataclasses import dataclass

from sfsbench.preset import discover, generate
from sfsbench.score import read_acc

PP_REPO = "/Users/vaceslaveliseev/@dev/structure-from-sherds-pp"
PREP_REPO = "/Users/vaceslaveliseev/@dev/SfSpp_preprocessing"
REFERENCE = os.path.join(PP_REPO, "Dataset", "SfS_pp")
REFERENCE_GT = os.path.join(REFERENCE, "Ground Truth")   # note the space
ICCV = "/Users/vaceslaveliseev/@dev/structure-from-sherds/ICCV Data"   # 2021 data

BREAKLINE_RE = re.compile(
    r"\[breakline\] (?P<path>\S+): (?P<segments>\d+) segment\(s\), "
    r"(?P<points>\d+) points, rim (?P<rim>-?\d+), base (?P<base>-?\d+) "
    r"=> info (?P<info>-?\d+)")


@dataclass
class PipelinePaths:
    pp_repo: str
    prep_repo: str
    dataset_root: str
    temp_root: str
    pot_id: str
    gt_root: str | None = None


def parse_stats(log, names=None):
    def one(pattern, default=None):
        m = re.search(pattern, log, re.MULTILINE)
        return int(m.group(1)) if m else default

    clusters_seen = [int(x) for x in
                     re.findall(r"Number of clusters is equal to (\d+)", log)]
    clusters = {}
    if names and len(clusters_seen) == len(names):
        # mesh_processing.cpp:1712 sorts its input files, so log order is name
        # order. If the counts disagree the mapping is unknowable; say nothing
        # rather than mislabel.
        clusters = dict(zip(names, clusters_seen))

    breaklines = {}
    for m in BREAKLINE_RE.finditer(log):
        name = os.path.basename(m.group("path"))[:14]
        breaklines[name] = {k: int(m.group(k))
                            for k in ("segments", "points", "rim", "base", "info")}

    return {
        "matches_raw": one(r"Total number : (\d+)"),
        "matches_pruned": one(r"Total number pruned : (\d+)"),
        "score": one(r"^Score : (\d+)"),          # anchored: visualize.cpp:414
        "clusters": clusters,
        "clusters_seen": clusters_seen,
        "breaklines": breaklines,
    }


def outcome(stats, acc):
    if not stats.get("matches_pruned") or stats.get("placed", 0) < 2:
        return "no-joins"
    if acc and acc["sherd"][1] and acc["sherd"][0] == acc["sherd"][1]:
        return "assembled"
    return "misassembled"


def _run(cmd, cwd, env, log_path):
    merged = dict(os.environ)
    merged.update(env or {})
    with open(log_path, "w") as f:
        subprocess.run(cmd, cwd=cwd, env=merged, stdout=f, stderr=subprocess.STDOUT,
                       check=False)
    return open(log_path, errors="replace").read()


def run_pipeline(paths, env=None, binary="SfSpp_gen"):
    env = dict(env or {})
    prep_env = {k: v for k, v in env.items() if k.startswith("SFSPP_")}
    prep_env["SFSPP_DATASET_ROOT"] = paths.dataset_root
    prep_env["SFSPP_TEMP_ROOT"] = paths.temp_root
    prep_env["SFSPP_POT_ID"] = paths.pot_id
    sfs_env = {k: v for k, v in env.items() if k.startswith("SFS_")}
    sfs_env["SFS_DATA_ROOT"] = paths.dataset_root + os.sep
    sfs_env["SFS_AUTOSAVE"] = "exit"

    logs = os.path.join(paths.temp_root, "logs")
    os.makedirs(logs, exist_ok=True)
    prep_build = os.path.join(paths.prep_repo, "build")
    pp_build = os.path.join(paths.pp_repo, "build")
    pot = paths.pot_id

    result_dir = os.path.join(paths.dataset_root, "Result")
    shutil.rmtree(result_dir, ignore_errors=True)     # filenames carry a minute
    for sub in ("Axes", f"Breaklines/Pot_{pot}", f"Surfaces/Pot_{pot}",
                "Result", "Graph Log"):
        os.makedirs(os.path.join(paths.dataset_root, sub), exist_ok=True)

    mesh_log = _run(["./MeshPreprocessing"], prep_build, prep_env,
                    os.path.join(logs, "mesh.log"))

    data_dir = os.path.join(paths.temp_root, "Data", f"Pot_{pot}")
    names = sorted({f[:14] for f in os.listdir(data_dir) if f.endswith("_Surface_0.xyz")})
    axes = {}
    for name in names:
        axis_path = os.path.join(paths.dataset_root, "Axes", f"{name}_Axis.xyz")
        _run(["./extract_axis",
              os.path.join(data_dir, f"{name}_Surface_0.xyz"),
              os.path.join(data_dir, f"{name}_Surface_1.xyz"),
              axis_path, "1000"],
             pp_build, {}, os.path.join(logs, f"axis_{name}.log"))
        # 1-10 lines: one per candidate surviving 10-degree dedup and the
        # SFS_AXIS_COST_TOLERANCE filter. More than one is a signal, not an error.
        axes[name] = sum(1 for line in open(axis_path) if line.strip())

    edge_log = _run(["./EdgeLineExtraction"], prep_build, prep_env,
                    os.path.join(logs, "edge.log"))

    # Regenerate the preset for whatever is now staged, then rebuild one TU.
    generate(discover(paths.dataset_root, [pot]),
             os.path.join(paths.pp_repo, "class", "data_path_generated.h"),
             gt_root=paths.gt_root)
    _run(["cmake", "--build", "build", "-j10", "--target", binary],
         paths.pp_repo, {}, os.path.join(logs, "build.log"))

    sfs_log = _run([f"./{binary}"], pp_build, sfs_env,
                   os.path.join(logs, "reassembly.log"))

    stats = parse_stats(mesh_log + "\n" + edge_log + "\n" + sfs_log, names)
    stats["axes"] = axes
    stats["names"] = names
    stats["result_dir"] = result_dir
    stats["placed"] = len([f for f in os.listdir(result_dir) if "_Top_1_OBJ_" in f])
    stats["acc"] = read_acc(result_dir)
    return stats
```

```python
# tools/sfsbench/cli.py
"""sfsbench — one command from raw meshes to a scored reassembly."""
import argparse, os, sys

from sfsbench.run import PipelinePaths, run_pipeline, outcome, PP_REPO, PREP_REPO
from sfsbench.score import score_result
from sfsbench.stage import stage_dataset

# Published Sherd Accuracy from the SfS++ project page, back-solved over the
# fragments the authors actually evaluate (the all-zero rows of
# Pot_<X>_simple_graph.txt are excluded, and disabled via shard_on_off).
# All ten are hard targets; see the plan-review note for the derivation.
REFERENCE_TARGETS = {
    "A": (8, 8), "B": (9, 9), "C": (4, 4), "D": (23, 28), "E": (31, 31),
    "F": (6, 6), "G": (7, 7), "H": (10, 11), "I": (25, 27), "J": (8, 11),
}


def main(argv=None):
    ap = argparse.ArgumentParser(prog="sfsbench")
    ap.add_argument("--src", required=True, help="folder of raw .ply/.obj fragments")
    ap.add_argument("--pot-id", required=True, help="single letter, e.g. A")
    ap.add_argument("--work", required=True, help="working directory for this run")
    ap.add_argument("--gt", default=None, help="ground-truth folder ('Ground Truth')")
    ap.add_argument("--source-mesh-dir", default=None,
                    help="folder holding the original Pot_X_Piece_NN_Mesh.obj used for scoring")
    ap.add_argument("--target-points", type=int, default=2_000_000)
    ap.add_argument("--binary", default="SfSpp_gen")
    ap.add_argument("--env", action="append", default=[], metavar="KEY=VALUE")
    args = ap.parse_args(argv)

    dataset_root = os.path.join(args.work, "Dataset")
    temp_root = os.path.join(args.work, "Temp")
    os.makedirs(temp_root, exist_ok=True)

    stage_dataset(args.src, dataset_root, args.pot_id, args.target_points)

    env = {}
    for item in args.env:
        k, _, v = item.partition("=")
        env[k] = v

    paths = PipelinePaths(PP_REPO, PREP_REPO, dataset_root, temp_root,
                          args.pot_id, args.gt)
    stats = run_pipeline(paths, env, args.binary)
    acc = stats["acc"]

    scored = {}
    if args.gt:
        scored = score_result(
            stats["result_dir"],
            args.source_mesh_dir or os.path.join(dataset_root, "Mesh", f"Pot_{args.pot_id}"),
            args.gt, args.pot_id)

    verdict = outcome(stats, acc)
    print(f"pot {args.pot_id}: {verdict}")
    print(f"  fragments placed   : {stats['placed']}")
    print(f"  matches raw/pruned : {stats['matches_raw']} / {stats['matches_pruned']}")
    print(f"  score              : {stats['score']}")
    print(f"  clusters           : {stats['clusters'] or stats['clusters_seen']}")
    print(f"  breakline info     : "
          f"{ {n: b['info'] for n, b in stats['breaklines'].items()} }")
    print(f"  axis candidates    : {stats['axes']}")
    if acc:
        target = REFERENCE_TARGETS.get(args.pot_id)
        line = (f"  sherd accuracy     : {acc['sherd'][0]} / {acc['sherd'][1]}")
        if target:
            line += f"   (target {target[0]}/{target[1]})"
        print(line)
        print(f"  edge accuracy      : {acc['edge'][0]} / {acc['edge'][1]}")
    else:
        print("  sherd accuracy     : not scored (no ground truth loaded)")
    if scored.get("total_scored"):
        print(f"  displacement       : pooled {scored['mean_error']:.2f} mm, "
              f"mean-of-means {scored['mean_of_means']:.2f} mm, "
              f"{scored['within_20mm']} / {scored['total_scored']} under 20 mm "
              f"(diagnostic, not the published metric)")
    return 0 if verdict == "assembled" else 1


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: Run test to verify it passes**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
.venv/bin/pytest sfsbench/tests/test_run.py -v
```

Expected: PASS, 7 passed

- [ ] **Step 5: Reproduce the known 0/8 baseline through the harness**

First isolate Pot A's meshes — the reference `Mesh/` folder is flat and holds all
164 files, so pointing `--src` at it would stage ten pots as one collection:

```bash
R="/Users/vaceslaveliseev/@dev/structure-from-sherds-pp/Dataset/SfS_pp"
mkdir -p /tmp/potA-src && cp "$R"/Mesh/Pot_A_Piece_0*_Mesh.obj /tmp/potA-src/
ls /tmp/potA-src | wc -l          # must print 8
```

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
PYTHONPATH=. .venv/bin/python -m sfsbench.cli \
  --src /tmp/potA-src --pot-id A --work /tmp/sfsbench-A \
  --gt "$R/Ground Truth" --source-mesh-dir "$R/Mesh"
```

Expected: `pot A: misassembled`, `sherd accuracy : 0 / 8`, matching the
diagnosis. The harness is correct when it reproduces the *known wrong* answer —
that proves it is driving the same pipeline, not a different one.

- [ ] **Step 6: Commit**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp
git add tools/sfsbench/run.py tools/sfsbench/cli.py tools/sfsbench/tests/test_run.py
git commit -m "harness: end-to-end runner with a three-way outcome"
```

---

### Task 6: Diagnosis — which defect is fatal

*(Was Task 4. Four runs instead of three; `$I` fixed; the fragment list is
derived, not hard-coded.)*

**Files:**
- Create: `tools/sfsbench/ablate.py`
- Create: `docs/superpowers/notes/2026-08-21-ablation-results.md`

**Interfaces:**
- Consumes: `run_pipeline`, `read_acc`
- Produces:
  - `substitute(dataset_root, pot_id, component, reference_root) -> list[str]` — symlinks the authors' files over ours for exactly one component (`"axes"`, `"breaklines"`, `"surfaces"`), returning the fragment names it acted on.
  - `substitute_rim_flags(dataset_root, pot_id, reference_root) -> dict[str, int]` — rewrites only the third field of the second header line of each of our `_Breakline_0.pcd` with the reference value, leaving geometry and segment ranges untouched.

This task produces **a decision, not a feature**. Its deliverable is the notes
file recording which run repaired the assembly.

The runs:

| Run | Surfaces | Axes | Breaklines | Answers |
|---|---|---|---|---|
| A | ours | **theirs** | ours | is the axis error fatal |
| B1 | ours | ours | ours + **their rim flags** | is the missing rim flag fatal |
| B2 | ours | ours | **theirs, whole file** | is anything else about our breaklines fatal |
| C | **theirs** | theirs (re-derived) | ours (or re-derived — say which) | is segmentation fatal |

Two things this table says that the previous revision's did not, and both change
how the result must be read: **Run C's axes are not ours** (substituting surfaces
forces `extract_axis` to be re-run on them, so C is confounded with A), and
**breaklines are downstream of surfaces** (`EdgeLineExtraction` consumes surfaces
and the axis), so unless `EdgeLineExtraction` is re-run in Run C the breaklines
there are still derived from our surfaces — a fourth, unnamed mixed state.

B1 exists because the previous Run B could not answer the question it was labelled
with. Substituting whole breakline files replaces geometry, segment ranges *and*
the flag at once — and measurement shows the geometry is already fine: segment
counts are **identical** for all eight fragments (5,3,4,4,4,2,3,3) and point
counts differ by at most 3.2 %. The only substantive difference is the third
header field. So overwriting just that field isolates the flag completely.

Reference values, read from `Dataset/SfS_pp/Breaklines/Pot_A_Piece_NN_Breakline_0.pcd`:

```
piece   1  2  3  4  5  6  7  8
theirs  0  1  1  1  0  1  0  1      (5 rims, no base flag)
ours    2  0  0  1  0  0  0  0      (1 rim, piece 04; base detected on piece 01)
```

`info` encodes `1` = rim, `2` = base, `3` = both, negative = base present but not
"sane" (`data_structure.cpp:218-225`).

- [ ] **Step 1: Write the substitution helpers**

```python
# tools/sfsbench/ablate.py
"""Swap one preprocessing product for the authors' own, to isolate the defect."""
import glob, os, re

COMPONENTS = {
    "axes":       ("Axes",          ["{n}_Axis.xyz"]),
    "breaklines": ("Breaklines",    ["{n}_Breakline_0.pcd"]),
    "surfaces":   ("Surfaces",      ["{n}_Surface_0.xyz", "{n}_Surface_1.xyz"]),
}


def _staged_names(dataset_root, pot_id):
    """Derive the fragment list from what is staged, not from range(1, 9)."""
    pattern = os.path.join(dataset_root, "Breaklines", f"Pot_{pot_id}",
                           f"Pot_{pot_id}_Piece_*_Breakline_0.pcd")
    names = sorted(os.path.basename(p)[:14] for p in glob.glob(pattern))
    if not names:
        raise FileNotFoundError(f"nothing staged under {pattern}")
    return names


def substitute(dataset_root, pot_id, component, reference_root):
    """Symlink the authors' files over ours for exactly one component.

    reference_root is Dataset/SfS_pp, whose Axes/, Breaklines/ and Surfaces/ are
    flat: they hold Pot_A_Piece_01_* directly, with no per-pot subdirectory,
    unlike the staged layout we produce (flat only for Axes/).
    """
    if component not in COMPONENTS:
        raise ValueError(f"unknown component {component!r}")
    top, templates = COMPONENTS[component]
    dst_dir = os.path.join(dataset_root, top) if component == "axes" else \
        os.path.join(dataset_root, top, f"Pot_{pot_id}")
    os.makedirs(dst_dir, exist_ok=True)
    src_dir = os.path.join(reference_root, top)

    names = _staged_names(dataset_root, pot_id)
    for n in names:
        for template in templates:
            f = template.format(n=n)
            src, dst = os.path.join(src_dir, f), os.path.join(dst_dir, f)
            if not os.path.exists(src):
                raise FileNotFoundError(src)
            if os.path.islink(dst) or os.path.exists(dst):
                os.remove(dst)
            os.symlink(src, dst)
    return names


def substitute_rim_flags(dataset_root, pot_id, reference_root):
    """Overwrite only the info byte of our breakline headers.

    Segment counts already match the reference exactly and point counts are
    within 3.2 %, so this isolates the rim/base flag from breakline geometry —
    which substituting whole files cannot do.
    """
    out = {}
    for n in _staged_names(dataset_root, pot_id):
        ours = os.path.join(dataset_root, "Breaklines", f"Pot_{pot_id}",
                            f"{n}_Breakline_0.pcd")
        theirs = os.path.join(reference_root, "Breaklines", f"{n}_Breakline_0.pcd")
        ref_line = open(theirs, errors="replace").read().splitlines()[1]
        ref_info = int(ref_line.split()[3])
        lines = open(ours, errors="replace").read().splitlines(keepends=True)
        head = lines[1].rstrip("\n").split()
        head[3] = str(ref_info)
        lines[1] = " ".join(head) + "\n"
        with open(ours, "w") as f:
            f.writelines(lines)
        out[n] = ref_info
    return out
```

- [ ] **Step 2: Run ablation A — their axes**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
R="/Users/vaceslaveliseev/@dev/structure-from-sherds-pp/Dataset/SfS_pp"
rm -rf /tmp/ablate-axes && cp -R /tmp/sfsbench-A /tmp/ablate-axes
PYTHONPATH=. .venv/bin/python -c "
from sfsbench.ablate import substitute
substitute('/tmp/ablate-axes/Dataset', 'A', 'axes', '$R')
"
```

*(The previous revision passed `'\$I'` here — an undefined variable that the shell
expanded to the empty string, making `reference_root` `''`.)*

Re-run only the reassembly, since preprocessing output is already staged:

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp
PYTHONPATH=tools /opt/homebrew/bin/python3 -c "
from sfsbench.preset import discover, generate
generate(discover('/tmp/ablate-axes/Dataset', ['A']),
         'class/data_path_generated.h', gt_root='$R/Ground Truth')
"
cmake --build build -j10 --target SfSpp_gen 2>&1 | grep -E "error:|Built target"
rm -rf /tmp/ablate-axes/Dataset/Result && mkdir -p /tmp/ablate-axes/Dataset/Result
cd build && SFS_DATA_ROOT=/tmp/ablate-axes/Dataset/ SFS_AUTOSAVE=exit ./SfSpp_gen 2>&1 \
  | grep -E "Sherd accuracy|Edge accuracy|Total number"
```

Record the Sherd Accuracy.

- [ ] **Step 3: Run ablation B1 — their rim flags only**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
rm -rf /tmp/ablate-rimflags && cp -RL /tmp/sfsbench-A /tmp/ablate-rimflags
PYTHONPATH=. .venv/bin/python -c "
from sfsbench.ablate import substitute_rim_flags
print(substitute_rim_flags('/tmp/ablate-rimflags/Dataset', 'A', '$R'))
"
```

Expected printout: `{... '_01': 0, '_02': 1, '_03': 1, '_04': 1, '_05': 0,
'_06': 1, '_07': 0, '_08': 1}`. Note `cp -RL` (dereference) so the copy does not
inherit symlinks from an earlier ablation. Then regenerate, rebuild, run and
record, as in Step 2.

- [ ] **Step 4: Run ablation B2 — their whole breakline files**

Same as Step 2 with `/tmp/ablate-breaklines` and
`substitute(..., 'breaklines', ...)`. Record the number.

- [ ] **Step 5: Run ablation C — their surfaces**

Same as Step 2 with `/tmp/ablate-surfaces` and `substitute(..., 'surfaces', ...)`,
then re-derive the axes, because they come from the surfaces:

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/build
for i in 01 02 03 04 05 06 07 08; do
  ./extract_axis /tmp/ablate-surfaces/Dataset/Surfaces/Pot_A/Pot_A_Piece_${i}_Surface_0.xyz \
                 /tmp/ablate-surfaces/Dataset/Surfaces/Pot_A/Pot_A_Piece_${i}_Surface_1.xyz \
                 /tmp/ablate-surfaces/Dataset/Axes/Pot_A_Piece_${i}_Axis.xyz 1000 > /dev/null 2>&1
done
```

**Decide and record whether `EdgeLineExtraction` is re-run here.** If it is not,
Run C's breaklines are still derived from *our* surfaces and the run tests only
the registration-time use of `surface_in`/`surface_out`. If it is, Run C is a
clean "their surfaces, everything downstream re-derived" run. Either is
defensible; leaving it unstated is not.

- [ ] **Step 6: Cheap extra — validate `extract_axis` across all ten pots**

The dataset ships `Ground Truth Axes/` (141 files) separately from `Axes/`.
`extract_axis` was validated at 0.35° against `Axes/` on Pot A alone. If
`Ground Truth Axes/` holds the manually corrected axes it is the better target,
and it covers ten pots for the cost of one loop:

```bash
ls "$R/Ground Truth Axes" | head -3
ls "$R/Ground Truth Axes" | wc -l
# then, per fragment, the angle between the axis we extract from THEIR surfaces
# and the ground-truth axis; report the median per pot.
```

Record the median per pot in the notes. This is evidence about the axis stage
across the whole dataset, which nothing else in the plan provides.

- [ ] **Step 7: Write the decision down**

Create `docs/superpowers/notes/2026-08-21-ablation-results.md` with the five
numbers (baseline 0/8, plus A, B1, B2, C), the axis-validation medians, and one
sentence naming the component that repaired the assembly. Tasks 7 and 8 are gated
on this result: implement only the one the ablation blames, and mark the other as
not-needed in the notes.

- [ ] **Step 8: Commit**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp
git add tools/sfsbench/ablate.py docs/superpowers/notes
git commit -m "diagnosis: isolate which preprocessing product breaks reassembly"
```

---

### Task 7: Port the reference rim criterion

*(Was Task 5. The verification vector is corrected, the insertion point is
corrected, and the commit is scoped to one hunk.)*

**Do this task only if Task 6 blamed the breaklines or the rim flags.**

**It does nothing for our own fragments** — they have no rim to detect. This task
serves the reference pots and the diagnosis; the path to `potB` runs through
Task 9 instead. Say so in the notes so nobody reads a rim improvement as progress
towards `potB`.

**Files:**
- Modify: `/Users/vaceslaveliseev/@dev/SfSpp_preprocessing/edgeline_extraction.cpp` — add `isBreaklineSegARim_Profile` among the existing static helpers (the helper block runs 2652–2829; `processFragmentData` begins at **2831**, not ~2650), and call it at the two `isBreaklineSegARim` call sites (**2927** and **3051** — the only two)

**Interfaces:**
- Produces: `static bool isBreaklineSegARim_Profile(pcl::PointCloud<pcl::PointXYZ>::Ptr seg, const std::string& axisPath)`

**Criterion, from `AxisExtraction/check_base_and_rim.m`** — a segment is a rim
when, measured against the axis:

- `std(r) < 1.0 mm` — radius near-constant along the segment
- `std(h) < 1.0 mm` — height near-constant, where `h = direction · point`
- `|mean(diff(r)[1:-1])| < 0.1` and `|mean(diff(h)[1:-1])| < 0.1`
- at least 20 points in the segment

Three semantics notes for the port, all confirmed against the MATLAB source: it
is `std`, not variance, despite the variable names `mean_r_err` / `mean_h_err`;
the threshold is on `r`, not `r²` (`compute_l2p_distance.m:12` returns a true
perpendicular distance); and `h` is the raw projection `dir · p` with no origin
subtraction, which is harmless for `std` and `diff`. The argument order is
`[direction; position]`, i.e. `vt(4:6)` is the direction. There is also a rim
tie-break at lines 123–130 that this port does not implement — note it if the
result disagrees with the reference on a fragment or two.

The MATLAB rim path is **not runnable as shipped** (`check_base_and_rim.m:57`
calls `align_point_cloud` and `read_axis.m:14` calls `load_root_dir` with a
hard-coded `D:/SFS_BB_temp/Plt_A`; neither function exists in the repository), so
there is nothing to diff the port against except the criterion itself and the
reference output files.

All four lengths are millimetres at the reference scale and must go through
`sfsEnvDouble` so Task 11 can rescale them.

- [ ] **Step 1: Write the implementation**

```cpp
// Port of the rim half of AxisExtraction/check_base_and_rim.m. The shipped C++
// heuristic (isBreaklineSegARim / _DistanceFromAxis) finds one rim on Pot A
// where the authors' own files mark five; rim pruning is what removes roughly a
// third of the candidate matches, all false.
static bool isBreaklineSegARim_Profile(pcl::PointCloud<pcl::PointXYZ>::Ptr seg,
                                       const std::string& axisPath) {
    const double r_threshold  = sfsEnvDouble("SFSPP_RIM_R_STD", 1.0);
    const double h_threshold  = sfsEnvDouble("SFSPP_RIM_H_STD", 1.0);
    const double dr_threshold = sfsEnvDouble("SFSPP_RIM_DR", 0.1);
    const double dh_threshold = sfsEnvDouble("SFSPP_RIM_DH", 0.1);
    const int    min_points   = static_cast<int>(sfsEnvDouble("SFSPP_RIM_MIN_POINTS", 20.0));

    if (!seg || static_cast<int>(seg->points.size()) < min_points) return false;

    std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> axes;
    if (!loadAxes(axisPath, axes)) return false;

    auto stdev = [](const std::vector<double>& v) {
        if (v.size() < 2) return std::numeric_limits<double>::max();
        const double m = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
        double s = 0.0;
        for (double x : v) s += (x - m) * (x - m);
        return std::sqrt(s / (v.size() - 1));       // MATLAB std, N-1
    };
    auto mean_inner_diff = [](const std::vector<double>& v) {
        if (v.size() < 4) return std::numeric_limits<double>::max();
        double s = 0.0; int n = 0;
        for (size_t i = 2; i + 1 < v.size(); ++i) { s += v[i] - v[i - 1]; ++n; }
        return n ? s / n : std::numeric_limits<double>::max();
    };

    for (const auto& ax : axes) {
        const Eigen::Vector3d org = ax.first, dir = ax.second;
        std::vector<double> r, h;
        r.reserve(seg->points.size());
        h.reserve(seg->points.size());
        for (const auto& p : seg->points) {
            const Eigen::Vector3d w(p.x - org.x(), p.y - org.y(), p.z - org.z());
            const double hi = w.dot(dir);
            h.push_back(hi);
            r.push_back((w - hi * dir).norm());
        }
        if (stdev(r) < r_threshold && stdev(h) < h_threshold &&
            std::abs(mean_inner_diff(r)) < dr_threshold &&
            std::abs(mean_inner_diff(h)) < dh_threshold) {
            return true;
        }
    }
    return false;
}
```

At both call sites, replace

```cpp
bool isRim = isBreaklineSegARim(cloud_breakLineSeg, outPath + fileNameOnly + "_SampledWithNormals.ply", segCount);
```

with

```cpp
// Prefer the reference criterion when an axis is available; the shipped
// heuristic stays as the fallback for axis-less runs.
const std::string rimAxisPath = datasetRoot() + "Axes/" + baseFileName + "_Axis.xyz";
bool isRim = std::filesystem::exists(rimAxisPath)
                 ? isBreaklineSegARim_Profile(cloud_breakLineSeg, rimAxisPath)
                 : isBreaklineSegARim(cloud_breakLineSeg,
                                      outPath + fileNameOnly + "_SampledWithNormals.ply", segCount);
```

At the second call site the surrounding code uses `outPathTemp` rather than
`outPath`; keep that name there. `baseFileName` is in scope in both functions.

- [ ] **Step 2: Build**

```bash
cd /Users/vaceslaveliseev/@dev/SfSpp_preprocessing
cmake --build build -j10 2>&1 | grep -E "error:|Built target"
```

Expected: `Built target EdgeLineExtraction`, no errors.

- [ ] **Step 3: Verify against the reference rim flags**

```bash
cd /Users/vaceslaveliseev/@dev/SfSpp_preprocessing/build
SFSPP_POT_ID=A SFSPP_DATASET_ROOT=/tmp/sfsbench-A/Dataset \
SFSPP_TEMP_ROOT=/tmp/sfsbench-A/Temp ./EdgeLineExtraction 2>&1 | grep '\[breakline\]'
R="/Users/vaceslaveliseev/@dev/structure-from-sherds-pp/Dataset/SfS_pp"
for i in 01 02 03 04 05 06 07 08; do
  echo "Piece_$i ours: $(sed -n 2p /tmp/sfsbench-A/Dataset/Breaklines/Pot_A/Pot_A_Piece_${i}_Breakline_0.pcd)  theirs: $(sed -n 2p "$R/Breaklines/Pot_A_Piece_${i}_Breakline_0.pcd")"
done
```

**The reference vector is `0, 1, 1, 1, 0, 1, 0, 1` for pieces 1–8**, i.e. five
rims and no base flag, read directly from
`Dataset/SfS_pp/Breaklines/Pot_A_Piece_NN_Breakline_0.pcd`. Expected: our third
field matches theirs for at least 6 of 8.

Two corrections behind that number:

- The previous revision expected `2, 1, 1, 1, 1, 1, 0, 0`. That vector is real,
  but it belongs to the **older ICCV collection**
  (`structure-from-sherds/ICCV Data/Breaklines/`), which also has different
  segment counts (4,3,4,4,3,2,3,3 against SfS++'s 5,3,4,4,4,2,3,3). Checking a
  SfS++ run against it would reject a correct port on pieces 1, 5, 7 and 8.
- The review that caught this proposed `2, 1, 1, 1, 0, 1, 0, 1`, keeping the
  leading `2`. That is also wrong: the SfS++ file has `0` there. The two
  collections genuinely disagree on the base flag — ICCV marks piece 1 as base,
  SfS++ marks no base on Pot A at all — which is consistent with `NO_BASE_INFO`
  being active in `main.cpp` and with the paper being base-agnostic. Both
  collections mark exactly five rims; they disagree on *which* (ICCV: 2,3,4,5,6;
  SfS++: 2,3,4,6,8).

Current state before the port, for comparison: `2, 0, 0, 1, 0, 0, 0, 0` — one
rim (piece 04) and a base on piece 01 from our own ported base detection.

- [ ] **Step 4: Re-run the full pipeline and score**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
R="/Users/vaceslaveliseev/@dev/structure-from-sherds-pp/Dataset/SfS_pp"
PYTHONPATH=. .venv/bin/python -m sfsbench.cli \
  --src /tmp/potA-src --pot-id A --work /tmp/sfsbench-A2 \
  --gt "$R/Ground Truth" --source-mesh-dir "$R/Mesh"
```

Expected: `sherd accuracy` strictly greater than 0/8. Record the number.

- [ ] **Step 5: Commit — the rim hunk only**

The `SFSPP_DATASET_ROOT` plumbing in this file was committed by Task 0, so the
working tree should now hold nothing but the rim change. Verify before committing:

```bash
cd /Users/vaceslaveliseev/@dev/SfSpp_preprocessing
git diff --stat                      # edgeline_extraction.cpp only
git diff | grep -c "^+.*SFSPP_DATASET_ROOT"   # must print 0
git add edgeline_extraction.cpp
git commit -m "Port the reference rim criterion from check_base_and_rim.m

The shipped heuristic finds one rim on Pot A where the authors' own files mark
five. Rim pruning removes roughly a third of the candidate matches, all false,
which is why our run starts from 224 matches where the reference starts from 76."
```

---

### Task 8: Segment the mesh instead of the point cloud

*(Was Task 6. Unchanged.)*

**Do this task only if Task 6 blamed the surfaces and Task 7 did not reach the target.**

**Files:**
- Create: `/Users/vaceslaveliseev/@dev/SfSpp_preprocessing/mesh_segmentation.cpp`
- Modify: `/Users/vaceslaveliseev/@dev/SfSpp_preprocessing/CMakeLists.txt` — add a `MeshSegmentation` executable alongside the existing two

**Rationale:** the authors discard the mesh and grow regions over a point cloud
using PCA normals from 5 neighbours spanning about a millimetre. Measured on our
scans, 46–68 % of neighbouring normal pairs exceed the 4.5° growing threshold, so
the surface shatters. Face normals on a triangulated mesh are exact and
connectivity is known, so a dihedral-angle cut is far more stable.

**Algorithm:**
1. Read the mesh; compute a unit normal per face.
2. Build face adjacency over shared edges (hash the sorted vertex-index pair).
3. Flood-fill faces into components, crossing an edge only when the angle between
   the two face normals is below `SFSPP_DIHEDRAL_DEG` (default 25.0).
4. Keep the two components with the largest total area — the inner and outer shell.
5. Emit them as `<name>_Surface_0.xyz` and `<name>_Surface_1.xyz`, one line per
   face centroid: `cx cy cz nx ny nz`.

**Interfaces:**
- Produces the same two files `extract_axis` and the reassembly already consume,
  so nothing downstream changes.

- [ ] **Step 1: Write the failing test**

Build a mesh whose answer is known: a box with the two large faces smooth and the
sides meeting them at 90°. Exactly two components should survive, one per large face.

```bash
cd /Users/vaceslaveliseev/@dev/SfSpp_preprocessing
mkdir -p tests
/opt/homebrew/bin/python3 - <<'PY'
import numpy as np
# a 100 x 100 x 5 slab, tessellated 20 x 20 on the large faces
n = 21
xs = np.linspace(0, 100, n); ys = np.linspace(0, 100, n)
V, F = [], []
for z in (0.0, 5.0):
    base = len(V)
    for y in ys:
        for x in xs:
            V.append((x, y, z))
    for j in range(n - 1):
        for i in range(n - 1):
            a = base + j * n + i; b = a + 1; c = a + n; d = c + 1
            F += [(a, b, d), (a, d, c)] if z == 0.0 else [(a, d, b), (a, c, d)]
with open("tests/slab.obj", "w") as f:
    for v in V: f.write("v %g %g %g\n" % v)
    for t in F: f.write("f %d %d %d\n" % (t[0]+1, t[1]+1, t[2]+1))
print("wrote tests/slab.obj:", len(V), "verts", len(F), "faces")
PY
```

- [ ] **Step 2: Run to verify it fails**

```bash
cd /Users/vaceslaveliseev/@dev/SfSpp_preprocessing
./build/MeshSegmentation tests/slab.obj /tmp/slab
```

Expected: FAIL — no such file or directory (the binary does not exist yet).

- [ ] **Step 3: Write the implementation**

```cpp
// mesh_segmentation.cpp
//
// Inner/outer shell segmentation by dihedral angle over mesh connectivity,
// replacing the point-cloud region growing in mesh_processing.cpp.
//
//   MeshSegmentation <mesh.obj> <out_prefix>
//
// writes <out_prefix>_Surface_0.xyz and <out_prefix>_Surface_1.xyz, each line
// "cx cy cz nx ny nz" for one face of that shell.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include <Eigen/Dense>

namespace {

double envDouble(const char* name, double fallback) {
    const char* v = std::getenv(name);
    if (!v || !*v) return fallback;
    try { return std::stod(v); } catch (...) { return fallback; }
}

struct Mesh {
    std::vector<Eigen::Vector3d> v;
    std::vector<std::array<int, 3>> f;
};

bool readObj(const std::string& path, Mesh& m) {
    std::ifstream in(path);
    if (!in) return false;
    std::string line;
    while (std::getline(in, line)) {
        if (line.rfind("v ", 0) == 0) {
            std::istringstream ss(line.substr(2));
            double x, y, z; ss >> x >> y >> z;
            m.v.emplace_back(x, y, z);
        } else if (line.rfind("f ", 0) == 0) {
            std::istringstream ss(line.substr(2));
            std::array<int, 3> t{};
            for (int k = 0; k < 3; ++k) {
                std::string tok; ss >> tok;
                t[k] = std::stoi(tok.substr(0, tok.find('/'))) - 1;
            }
            m.f.push_back(t);
        }
    }
    return !m.v.empty() && !m.f.empty();
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "usage: MeshSegmentation <mesh.obj> <out_prefix>\n";
        return 2;
    }
    const double max_dihedral = envDouble("SFSPP_DIHEDRAL_DEG", 25.0);
    const double cos_limit = std::cos(max_dihedral * M_PI / 180.0);

    Mesh m;
    if (!readObj(argv[1], m)) { std::cerr << "cannot read " << argv[1] << "\n"; return 1; }

    const size_t nf = m.f.size();
    std::vector<Eigen::Vector3d> normal(nf), centroid(nf);
    std::vector<double> area(nf);
    for (size_t i = 0; i < nf; ++i) {
        const Eigen::Vector3d& a = m.v[m.f[i][0]];
        const Eigen::Vector3d& b = m.v[m.f[i][1]];
        const Eigen::Vector3d& c = m.v[m.f[i][2]];
        Eigen::Vector3d n = (b - a).cross(c - a);
        area[i] = 0.5 * n.norm();
        normal[i] = (area[i] > 0) ? n.normalized() : Eigen::Vector3d::UnitZ();
        centroid[i] = (a + b + c) / 3.0;
    }

    std::map<std::pair<int, int>, std::vector<int>> edge_faces;
    for (size_t i = 0; i < nf; ++i)
        for (int k = 0; k < 3; ++k) {
            int a = m.f[i][k], b = m.f[i][(k + 1) % 3];
            edge_faces[{std::min(a, b), std::max(a, b)}].push_back(static_cast<int>(i));
        }

    std::vector<std::vector<int>> adj(nf);
    for (const auto& kv : edge_faces)
        if (kv.second.size() == 2) {
            const int x = kv.second[0], y = kv.second[1];
            if (normal[x].dot(normal[y]) >= cos_limit) {
                adj[x].push_back(y);
                adj[y].push_back(x);
            }
        }

    std::vector<int> comp(nf, -1);
    std::vector<double> comp_area;
    for (size_t s = 0; s < nf; ++s) {
        if (comp[s] != -1) continue;
        const int id = static_cast<int>(comp_area.size());
        comp_area.push_back(0.0);
        std::vector<int> stack{static_cast<int>(s)};
        comp[s] = id;
        while (!stack.empty()) {
            const int cur = stack.back(); stack.pop_back();
            comp_area[id] += area[cur];
            for (int nb : adj[cur])
                if (comp[nb] == -1) { comp[nb] = id; stack.push_back(nb); }
        }
    }

    std::vector<int> order(comp_area.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](int a, int b) { return comp_area[a] > comp_area[b]; });

    std::cout << "faces " << nf << ", components " << comp_area.size();
    if (!order.empty()) std::cout << ", largest area " << comp_area[order[0]];
    std::cout << std::endl;
    if (order.size() < 2) { std::cerr << "fewer than two shells found\n"; return 1; }

    for (int k = 0; k < 2; ++k) {
        std::ofstream out(std::string(argv[2]) + "_Surface_" + std::to_string(k) + ".xyz");
        out.precision(8);
        for (size_t i = 0; i < nf; ++i)
            if (comp[i] == order[k])
                out << centroid[i].x() << ' ' << centroid[i].y() << ' ' << centroid[i].z()
                    << ' ' << normal[i].x() << ' ' << normal[i].y() << ' ' << normal[i].z()
                    << '\n';
    }
    return 0;
}
```

Add to `CMakeLists.txt`, after the two existing `add_executable` lines:

```cmake
add_executable(MeshSegmentation mesh_segmentation.cpp)
target_link_libraries(MeshSegmentation PRIVATE ${PCL_LIBRARIES})
```

- [ ] **Step 4: Run to verify it passes**

```bash
cd /Users/vaceslaveliseev/@dev/SfSpp_preprocessing
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release > /dev/null && cmake --build build -j10 --target MeshSegmentation 2>&1 | grep -E "error:|Built target"
./build/MeshSegmentation tests/slab.obj /tmp/slab
wc -l /tmp/slab_Surface_0.xyz /tmp/slab_Surface_1.xyz
```

Expected: `Built target MeshSegmentation`, then `components` at least 2, and both
surface files holding 800 lines each (the slab's two large faces, 20×20×2 triangles).

- [ ] **Step 5: Compare against the authors' surfaces on Pot A**

```bash
cd /Users/vaceslaveliseev/@dev/SfSpp_preprocessing
R="/Users/vaceslaveliseev/@dev/structure-from-sherds-pp/Dataset/SfS_pp"
for i in 01 02 03 04 05 06 07 08; do
  ./build/MeshSegmentation "$R/Mesh/Pot_A_Piece_${i}_Mesh.obj" /tmp/segA_${i}
done
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/build
for i in 01 02 03 04 05 06 07 08; do
  ./extract_axis /tmp/segA_${i}_Surface_0.xyz /tmp/segA_${i}_Surface_1.xyz /tmp/segA_${i}_Axis.xyz 1000 2>/dev/null | tail -1
done
```

Then compare each `/tmp/segA_NN_Axis.xyz` against `$R/Axes/Pot_A_Piece_NN_Axis.xyz`
(and against `Ground Truth Axes/`, per Task 6 Step 6).

Expected: median angle below 1°, against 3.5° from the point-cloud segmentation.
If it is not better, stop and reconsider rather than wiring it in.

- [ ] **Step 6: Commit**

```bash
cd /Users/vaceslaveliseev/@dev/SfSpp_preprocessing
git add mesh_segmentation.cpp CMakeLists.txt tests/slab.obj
git commit -m "Segment the shell over mesh connectivity, not the point cloud

Point-cloud region growing uses PCA normals from 5 neighbours spanning about a
millimetre; 46-68% of neighbouring normal pairs exceed the 4.5 degree growing
threshold, so the surface shatters. Face normals are exact and connectivity is
known, so a dihedral-angle cut is stable."
```

---

### Task 9: The rimless run — is our own data even in scope

*(New. This is the only evidence that says whether rimless four-fragment
collections are in scope at all, and it costs one rebuild and three small pots.)*

Our four scans have neither rim nor base. `main.cpp:31` carries
`//#define NO_RIM_INFO`, commented out — the exact switch for that case — and
nothing in the plan used it. Meanwhile `NO_BASE_INFO` is already active, so base
information never reaches reassembly regardless.

Run it **before** spending anything on `potB`. Paper Table X's "Non-enhanced"
ablation suggests rim-less reassembly degrades but does not collapse; this
measures it on data with ground truth.

**Files:**
- Modify: `CMakeLists.txt` — a `SfSpp_norim` target
- Create: `docs/superpowers/notes/2026-08-21-rimless.md`

- [ ] **Step 1: Add the target**

Rather than editing `main.cpp`, guard the existing define so the build system can
set it:

```cpp
//#define NO_RIM_INFO      -> becomes:
#ifdef SFS_NO_RIM
#define NO_RIM_INFO
#endif
#define NO_BASE_INFO
```

and in `CMakeLists.txt`, clone the `SfSpp_gen` block as `SfSpp_norim` with
`POT_GEN SFS_NO_RIM` in its `target_compile_definitions`.

- [ ] **Step 2: Sweep the three small pots both ways**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release > /dev/null && cmake --build build -j10 2>&1 | grep -E "error:|Built target"
# for pot in A B C: stage, generate, build, run twice - SfSpp_gen and SfSpp_norim
```

Use the sweep driver from Task 10 if it is already written; otherwise a shell
loop over the Task 5 CLI with `--binary SfSpp_gen` and `--binary SfSpp_norim`.

- [ ] **Step 3: Record the verdict**

Write `docs/superpowers/notes/2026-08-21-rimless.md` with six numbers (A, B, C ×
with/without rim information) and one sentence: whether a rimless collection can
still reach the published accuracy, and therefore whether `potB` is in scope for
the pipeline as it stands.

- [ ] **Step 4: Commit**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp
git add CMakeLists.txt main.cpp docs/superpowers/notes
git commit -m "Measure reassembly without rim information

Our own fragments have no rim, so the rim criterion port does nothing for them.
NO_RIM_INFO already existed as the switch for that case and nothing used it."
```

---

### Task 10: The reference sweep (ten pots)

*(Was Task 8. Rebuilt on the generated preset — the previous version could not
have run: it globbed `_Mesh.obj` when 147 of 164 reference meshes are
`_Mesh_DS.obj`, compiled the preprocessing for one pot, and pointed ten new
binaries at a flat layout the harness never produces.)*

**Files:**
- Create: `tools/sfsbench/sweep.py`
- Modify: `docs/superpowers/notes/2026-08-21-ablation-results.md` — append the sweep table

**Interfaces:**
- Consumes: everything above
- Produces: `sweep(pots, work_root, resume=True) -> dict[str, dict]` — runs the full chain per pot, checkpointing each result to `<work_root>/<pot>/result.json` so an interruption does not lose the run.

- [ ] **Step 1: Write the sweep**

```python
# tools/sfsbench/sweep.py
"""Run every reference pot through the whole chain and tabulate the result."""
import glob, json, os, shutil

from sfsbench.cli import REFERENCE_TARGETS
from sfsbench.run import (REFERENCE, REFERENCE_GT, PP_REPO, PREP_REPO,
                          PipelinePaths, run_pipeline, outcome)
from sfsbench.score import score_result
from sfsbench.stage import stage_dataset

MESH_SUFFIXES = ("_Mesh.obj", "_Mesh_DS.obj")


def _pot_source(pot, scratch):
    """The reference Mesh folder is flat and holds all ten pots; isolate one.

    Pots A and B ship _Mesh.obj, C-J ship _Mesh_DS.obj - 147 of the 164 files.
    Both suffixes must be accepted here; stage_dataset normalises the name on
    the way out, because mesh_processing.cpp composes <stem>_Mesh.obj.
    """
    out = os.path.join(scratch, f"src_{pot}")
    shutil.rmtree(out, ignore_errors=True)
    os.makedirs(out, exist_ok=True)
    found = []
    for suffix in MESH_SUFFIXES:
        found += glob.glob(os.path.join(REFERENCE, "Mesh", f"Pot_{pot}_Piece_*{suffix}"))
    if not found:
        raise FileNotFoundError(f"no meshes for pot {pot} under {REFERENCE}/Mesh")
    for src in sorted(found):
        shutil.copy(src, out)
    return out


def sweep(pots, work_root, resume=True):
    results = {}
    for pot in pots:
        work = os.path.join(work_root, pot)
        checkpoint = os.path.join(work, "result.json")
        if resume and os.path.exists(checkpoint):
            results[pot] = json.load(open(checkpoint))
            print(f"pot {pot}: resumed from checkpoint")
            continue

        dataset_root = os.path.join(work, "Dataset")
        temp_root = os.path.join(work, "Temp")
        os.makedirs(temp_root, exist_ok=True)
        names = stage_dataset(_pot_source(pot, work_root), dataset_root, pot)
        paths = PipelinePaths(PP_REPO, PREP_REPO, dataset_root, temp_root, pot,
                              REFERENCE_GT)
        stats = run_pipeline(paths, {}, binary="SfSpp_gen")
        scored = score_result(stats["result_dir"], os.path.join(REFERENCE, "Mesh"),
                              REFERENCE_GT, pot)
        record = {
            "acc": stats["acc"],
            "outcome": outcome(stats, stats["acc"]),
            "target": REFERENCE_TARGETS.get(pot),
            "staged": len(names),
            "placed": stats["placed"],
            "displacement_mm": scored.get("mean_error"),
        }
        json.dump(record, open(checkpoint, "w"), indent=2)
        results[pot] = record
        print(f"pot {pot}: {record['outcome']} {record['acc']}")
    return results


def format_table(results):
    lines = ["| Pot | sherd acc | target | edge acc | outcome | displacement |",
             "|---|---|---|---|---|---|"]
    for pot, r in results.items():
        acc = r.get("acc") or {}
        s = "/".join(map(str, acc.get("sherd", ("-", "-"))))
        e = "/".join(map(str, acc.get("edge", ("-", "-"))))
        t = "/".join(map(str, r["target"])) if r.get("target") else "-"
        d = r.get("displacement_mm")
        lines.append(f"| {pot} | {s} | {t} | {e} | {r['outcome']} | "
                     f"{'-' if d is None else format(d, '.1f') + ' mm'} |")
    return "\n".join(lines)
```

- [ ] **Step 2: Run the sweep**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
PYTHONPATH=. .venv/bin/python -c "
from sfsbench.sweep import sweep, format_table
r = sweep(['C','A','B','F','G','J','H','D','E','I'], '/tmp/sfs-sweep')
print(format_table(r))
"
```

Targets, all ten gated:

| Pot | evaluated | target SA |
|---|---|---|
| A | 8 | 8/8 |
| B | 9 | 9/9 |
| C | 4 | 4/4 |
| D | 28 | 23/28 |
| E | 31 | 31/31 |
| F | 6 | 6/6 |
| G | 7 | 7/7 |
| H | 11 | 10/11 |
| I | 27 | 25/27 |
| J | 11 | 8/11 |

The order is deliberate: the small pots finish quickly and fail fast, while D, E
and I hold 28–31 evaluated fragments each. **Budget a working day, not "a couple
of hours"** — ICCV 2021 Table 2, on the authors' hardware, puts preprocessing
alone at 47.8 min (D) and 37.2 min (E), with reassembly at (b=3, k=5) a further
22.5 and 17.0 min; Pot I is a third pot of that size, and seven more pots follow.
The per-pot checkpoint means an interruption costs one pot, not the sweep.

Gate on **Sherd Accuracy**, not Edge Accuracy: Pot D's published EA of 83.3 % does
not back-solve over the 69 edges in its graph.

- [ ] **Step 3: Append the table to the notes and commit**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp
git add tools/sfsbench/sweep.py docs/superpowers/notes
git commit -m "harness: reference sweep over all ten pots, checkpointed per pot"
```

---

### Task 11: The length-scale envelope

*(Was Task 7, and moved after the diagnosis. The old design inferred a scale
proportional to wall thickness and checked it by "the reference pots must come
out at 1.0". Measured across all ten pots, that check fails on nine of them: the
thresholds are fixed and work from 2.80 mm to 8.20 mm median wall thickness, so
they are demonstrably not thickness-proportional.)*

**Do this task only after the diagnosis has cleared surface segmentation**, or
derive thickness from the mesh rather than from `Surface_0`/`Surface_1`.
Calibrating the thresholds by which segmentation's output is judged, on a
measurement taken from that output, is circular.

**Files:**
- Create: `tools/sfsbench/scale.py`
- Create: `tools/sfsbench/tests/test_scale.py`
- Modify: `tools/sfsbench/cli.py` — warn, and optionally set `SFS_LENGTH_SCALE`, when a collection falls outside the envelope

**Interfaces:**
- Produces:
  - `wall_thickness(surface0_path, surface1_path, samples=3000) -> float` — median nearest-neighbour distance, inner surface to outer
  - `mesh_thickness(mesh_path, samples=3000) -> float` — the same quantity taken across the triangulation, needing no segmentation
  - `envelope_check(thickness: float, extent: float) -> dict` — `{"inside": bool, "thickness_scale": float, "extent_scale": float, "reason": str}`; both scales are 1.0 whenever `inside` is true

**The envelope**, measured across all 163 reference surface pairs and 164 meshes:

```python
# Per-fragment ranges over pots A-J, every one of which reaches its published
# accuracy at SFS_LENGTH_SCALE = 1.0. Inside this box, scale is 1.0 by
# definition; outside it, the fixed thresholds have no evidence behind them.
THICKNESS_RANGE_MM = (2.4, 13.5)
EXTENT_RANGE_MM = (22.0, 211.0)

# Per-pot medians, for reference and for the unit test:
#   thickness  A 3.67  B 3.44  C 6.07  D 8.20  E 6.79
#              F 4.72  G 2.80  H 5.87  I 7.51  J 5.53      (2.93x spread)
#   extent     A 93.4  B 78.7  C 86.8  D 79.8  E 73.5
#              F 35.6  G 54.4  H 72.1  I 86.3  J 40.9      (2.63x spread)
REFERENCE_THICKNESS_MM = 3.67    # Pot A median, measured 3.6731
REFERENCE_EXTENT_MM = 93.4       # Pot A median max-extent, measured 93.42
```

Two constants corrected from the previous revision:
`REFERENCE_EXTENT_MM = 119.0` was documented as "median bbox diagonal-max across
the 8 fragments" and is neither — Pot A's median max-extent is 93.42 and its
median bbox diagonal 126.50; 119.0 is `Pot_A_Piece_01`'s max extent, the largest
fragment, at the 89th percentile of all 164 reference fragments. And the claim
that "piece 1, the base fragment, is 7.26 mm" is dropped: it is not reproducible
under any definition tried, and piece 1's distance distribution is bimodal.

Also note before wiring anything: **`SFS_LENGTH_SCALE` moves 12 of roughly 57
length-valued constants.** `sfsLen()` / `sfsArea()` wrap exactly 12 literals
across 11 lines in `feature_matching.cpp` and `ranking_system.cpp`; roughly 45
comparable millimetre constants are untouched — `RejectOutlier` distance gates
(20 mm at reconstruction.cpp:460, 524, 587, 1187, 1298; 10 mm at 660, 685),
`CauchyLoss(5.0)` on `loss_dist` and `CauchyLoss(2.0)` on `loss_axis`/`loss_rim`,
`isConverge` translation terms (2.0 / 1.0 / 0.5 mm), `FillMatchedPoints`'
`dist < 3.0`, `ranking_system.cpp:20` `t_threshold = 20.0`, the LCS feature
quantisation steps at `feature_matching.cpp:1518-1521`, `filter.cpp:248`,
`data_structure.h:271, 275`. On a collection whose scale genuinely differs,
moving 12 of 57 leaves the pipeline internally inconsistent — arguably worse than
moving none. Decide explicitly: extend `sfsLen` to the rest, or document the knob
as partial and diagnostic.

- [ ] **Step 1: Write the failing test**

```python
# tools/sfsbench/tests/test_scale.py
import os, tempfile
import numpy as np
from sfsbench.scale import (wall_thickness, envelope_check,
                            THICKNESS_RANGE_MM, EXTENT_RANGE_MM)

POT_MEDIANS = {           # measured across all ten reference pots
    "A": (3.67, 93.4), "B": (3.44, 78.7), "C": (6.07, 86.8),
    "D": (8.20, 79.8), "E": (6.79, 73.5), "F": (4.72, 35.6),
    "G": (2.80, 54.4), "H": (5.87, 72.1), "I": (7.51, 86.3),
    "J": (5.53, 40.9),
}

def _write_plane(path, z, n=60, span=50.0):
    xs = np.linspace(0, span, n)
    g = np.stack(np.meshgrid(xs, xs), -1).reshape(-1, 2)
    rows = np.hstack([g, np.full((len(g), 1), z), np.tile([0., 0, 1], (len(g), 1))])
    np.savetxt(path, rows, fmt="%.6g")

def test_wall_thickness_of_two_parallel_planes():
    with tempfile.TemporaryDirectory() as d:
        a = os.path.join(d, "s0.xyz"); b = os.path.join(d, "s1.xyz")
        _write_plane(a, 0.0); _write_plane(b, 4.0)
        assert abs(wall_thickness(a, b) - 4.0) < 0.05

def test_wall_thickness_finds_a_neighbour_beyond_the_first_ring():
    """The old cell size was span/40 with a 27-cell search, so a true nearest
    neighbour further than one ring was missed and its point silently dropped."""
    with tempfile.TemporaryDirectory() as d:
        a = os.path.join(d, "s0.xyz"); b = os.path.join(d, "s1.xyz")
        _write_plane(a, 0.0, span=200.0); _write_plane(b, 40.0, span=200.0)
        assert abs(wall_thickness(a, b) - 40.0) < 0.5

def test_every_reference_pot_is_inside_the_envelope():
    """The real correctness criterion: all ten already reach their published
    accuracy at scale 1.0, so inference must be a no-op for all ten."""
    for pot, (thickness, extent) in POT_MEDIANS.items():
        out = envelope_check(thickness, extent)
        assert out["inside"], f"{pot} fell outside: {out['reason']}"
        assert out["thickness_scale"] == 1.0
        assert out["extent_scale"] == 1.0

def test_our_own_fragments_are_outside_the_envelope():
    out = envelope_check(37.6, 508.5)     # measured on data/potB
    assert not out["inside"]
    assert out["thickness_scale"] > 5.0
```

*(The previous revision's two unit tests fed the reference constants back into
`infer_scales` and asserted the result was 1.0 — tautologies. Task 7 Step 5 was
tautological too: it passed `[119.0]` as the extent list, so `extent_scale` came
out 1.0 by construction, and it only ever ran on Pot A.)*

- [ ] **Step 2: Run test to verify it fails**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
.venv/bin/pytest sfsbench/tests/test_scale.py -v
```

Expected: FAIL — `ModuleNotFoundError: No module named 'sfsbench.scale'`

- [ ] **Step 3: Write the implementation**

```python
# tools/sfsbench/scale.py
"""Decide whether a collection is inside the envelope the fixed thresholds cover.

Measured across pots A-J: median wall thickness spans 2.80-8.20 mm (2.93x) and
median max-extent 35.6-93.4 mm (2.63x), and ALL TEN reach their published
accuracy at SFS_LENGTH_SCALE = 1.0. So the thresholds are not proportional to
wall thickness - they are fixed, and they work across that whole range. Inferring
a proportional scale would return 1.0 on Pot A by construction and 0.76-2.23 on
the other nine, perturbing nine pots that currently work.
"""
import numpy as np

THICKNESS_RANGE_MM = (2.4, 13.5)     # per-fragment min-max over A-J
EXTENT_RANGE_MM = (22.0, 211.0)
REFERENCE_THICKNESS_MM = 3.67        # Pot A median (3.6731)
REFERENCE_EXTENT_MM = 93.4           # Pot A median max-extent (93.42)


def _nearest_distances(P, Q, cell=None, samples=3000, seed=0):
    """Nearest-neighbour distances P -> Q via a voxel hash with ring expansion.

    A KD-tree would be natural but scipy is not available here. The ring grows
    until the best distance found is provably nearer than the next ring, so a
    neighbour further than one cell is neither missed nor silently dropped - the
    previous version searched only the 27 adjacent cells at cell = span/40
    (~3 mm on a 120 mm fragment), so a true neighbour at 3.7 mm could be missed
    and its sample discarded from the median.
    """
    Q = np.asarray(Q, dtype=np.float64)
    P = np.asarray(P, dtype=np.float64)
    if cell is None:
        span = float(np.max(Q.max(0) - Q.min(0)))
        cell = max(span / 40.0, 1e-6)
    origin = Q.min(0)
    k = np.floor((Q - origin) / cell).astype(np.int64)
    dims = k.max(0) + 1
    h = k[:, 0] + dims[0] * (k[:, 1] + dims[1] * k[:, 2])
    order = np.argsort(h)
    hs = h[order]
    max_ring = int(np.max(dims)) + 1

    rng = np.random.default_rng(seed)
    picks = rng.choice(len(P), min(samples, len(P)), replace=False)
    out = []
    for i in picks:
        kk = np.floor((P[i] - origin) / cell).astype(np.int64)
        best = np.inf
        ring = 0
        while ring <= max_ring:
            cand = []
            for dx in range(-ring, ring + 1):
                for dy in range(-ring, ring + 1):
                    for dz in range(-ring, ring + 1):
                        if ring and max(abs(dx), abs(dy), abs(dz)) != ring:
                            continue           # only the new shell
                        c = kk + np.array([dx, dy, dz])
                        if (c < 0).any() or (c >= dims).any():
                            continue
                        hh = c[0] + dims[0] * (c[1] + dims[1] * c[2])
                        lo = np.searchsorted(hs, hh, "left")
                        hi = np.searchsorted(hs, hh, "right")
                        if hi > lo:
                            cand.append(order[lo:hi])
            if cand:
                idx = np.concatenate(cand)
                best = min(best, float(np.linalg.norm(Q[idx] - P[i], axis=1).min()))
            if best <= ring * cell:            # nothing outside can be nearer
                break
            ring += 1
        if np.isfinite(best):
            out.append(best)
    return np.asarray(out)


def wall_thickness(surface0_path, surface1_path, samples=3000):
    A = np.loadtxt(surface0_path)[:, :3]
    B = np.loadtxt(surface1_path)[:, :3]
    d = _nearest_distances(A, B, samples=samples)
    return float(np.median(d)) if len(d) else float("nan")


def mesh_thickness(mesh_path, samples=3000):
    """Shell-to-shell distance across the triangulation.

    Needed because Surface_0/Surface_1 are the output of surface segmentation,
    which the design lists as a suspect for the total failure: calibrating on it
    would be circular. Split faces by normal against the fragment's dominant
    axis, then measure one side against the other.
    """
    raise NotImplementedError("write this if Task 6 does not clear segmentation")


def envelope_check(thickness, extent):
    """Scales are 1.0 inside the envelope, by definition. Outside, they are the
    ratio to Pot A - and a warning, because 45 of ~57 length constants ignore
    SFS_LENGTH_SCALE (see the plan)."""
    t_lo, t_hi = THICKNESS_RANGE_MM
    e_lo, e_hi = EXTENT_RANGE_MM
    reasons = []
    if not (t_lo <= thickness <= t_hi):
        reasons.append(f"thickness {thickness:.2f} outside {t_lo}-{t_hi} mm")
    if not (e_lo <= extent <= e_hi):
        reasons.append(f"extent {extent:.1f} outside {e_lo}-{e_hi} mm")
    if not reasons:
        return {"inside": True, "thickness_scale": 1.0, "extent_scale": 1.0,
                "reason": "inside the reference envelope"}
    return {"inside": False,
            "thickness_scale": thickness / REFERENCE_THICKNESS_MM,
            "extent_scale": extent / REFERENCE_EXTENT_MM,
            "reason": "; ".join(reasons)}
```

In `cli.py`, report the check and **do not** set `SFS_LENGTH_SCALE` silently:

```python
from sfsbench.scale import envelope_check, wall_thickness

check = envelope_check(thickness, extent)
print(f"  scale envelope     : {check['reason']}")
if not check["inside"]:
    print(f"  suggested scales   : thickness x{check['thickness_scale']:.2f}, "
          f"extent x{check['extent_scale']:.2f}  "
          f"(SFS_LENGTH_SCALE covers 12 of ~57 length constants -- "
          f"pass --env SFS_LENGTH_SCALE=... deliberately)")
```

- [ ] **Step 4: Run test to verify it passes**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
.venv/bin/pytest sfsbench/tests/test_scale.py -v
```

Expected: PASS, 4 passed

- [ ] **Step 5: Verify on all ten reference pots, from the files**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
R="/Users/vaceslaveliseev/@dev/structure-from-sherds-pp/Dataset/SfS_pp"
PYTHONPATH=. .venv/bin/python -c "
import glob, os, numpy as np
from sfsbench.scale import wall_thickness, envelope_check
from sfsbench.stage import read_mesh
R='$R'
for pot in 'ABCDEFGHIJ':
    pairs=sorted(glob.glob(f'{R}/Surfaces/Pot_{pot}_Piece_*_Surface_0.xyz'))
    t=[wall_thickness(p, p.replace('_Surface_0','_Surface_1')) for p in pairs]
    ext=[]
    for m in sorted(glob.glob(f'{R}/Mesh/Pot_{pot}_Piece_*_Mesh*.obj')):
        V,_=read_mesh(m); ext.append(float(np.max(V.max(0)-V.min(0))))
    c=envelope_check(float(np.median(t)), float(np.median(ext)))
    print(pot, round(float(np.median(t)),2), round(float(np.median(ext)),1), c['inside'], c['reason'])
"
```

Expected: **`True` for all ten**, with medians close to the table above. Anything
else means the envelope is drawn wrong — fix the envelope, not the pots.
*(This replaces the previous Step 5, which ran on Pot A alone and passed
`[119.0]` as the extent, making the result 1.0 by construction.)*

- [ ] **Step 6: Commit**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp
git add tools/sfsbench/scale.py tools/sfsbench/tests/test_scale.py tools/sfsbench/cli.py
git commit -m "harness: check a collection against the reference scale envelope

Measured across pots A-J, median wall thickness spans 2.80-8.20 mm and all ten
reach their published accuracy at SFS_LENGTH_SCALE = 1.0, so the thresholds are
fixed rather than thickness-proportional. Inference applies only outside that
envelope."
```

---

### Task 12: Run our own fragments with everything in place

*(Was Task 9. A physical measurement is now a precondition, the pairwise distance
computation is rewritten, and `--pot-id A` is explained.)*

**Files:**
- Create: `docs/superpowers/notes/2026-08-21-potb-result.md`

- [ ] **Step 1: Establish the physical scale — before anything is run**

The four scans are 371–680 units in max extent with a 37.6-unit wall, against
reference fragments whose median max extent is 74.7 and whose thickest wall is
13.4 mm. The PLY headers say `Created in RealityCapture` (three files) and `by
Geomagic Studio` (one) — photogrammetry output, **unscaled unless a scale
constraint was supplied at reconstruction**. So the coordinates may not be
millimetres.

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
PYTHONPATH=. .venv/bin/python -c "
import glob, numpy as np
from sfsbench.stage import read_mesh
for f in sorted(glob.glob('../data/potB/*.ply')):
    V,_ = read_mesh(f)
    e = V.max(0)-V.min(0)
    print(f.split('/')[-1], len(V), 'extent', np.round(e,1), 'max', round(float(e.max()),1))
"
```

Then **measure one sherd with a caliper** and compare. Record the ratio in the
notes.

Without this, a "no joins found" verdict is uninterpretable: it could equally
mean "the correspondence window is 5 units on a 500-unit object". Do not skip it
and do not substitute a guess. If the scans turn out not to be in millimetres,
rescale the meshes at staging time — that is far safer than moving 12 of 57
thresholds (Task 11).

Two more preconditions worth checking here:

- Read the Task 9 verdict. If a rimless reference pot cannot reach its target,
  say so before running and treat the outcome accordingly.
- Do **not** use `data/test_fragments` for any geometric statistic. It is the
  same four scans decimated ~12×, but the decimated copies carry stray points
  30–55 units off-surface in z (z-extent 145.15 vs 89.76 for FY234021; 131.39 vs
  91.64 for FY234094; 103.41 vs 71.47 for FY234104). Clean them first, or use
  `data/potB`.

- [ ] **Step 2: Run the full-resolution fragments**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
PYTHONPATH=. .venv/bin/python -m sfsbench.cli \
  --src ../data/potB --pot-id A --work /tmp/sfsbench-potB
```

`--pot-id A` is not arbitrary and is not about the vessel: it names the staged
directory (`Point/Pot_A/`, `Mesh/Pot_A/`, …) and the fragment filenames
(`Pot_A_Piece_NN`), which must be 14 characters for
`mesh_processing.cpp:338`'s `substr(0, 14)`. With the generated preset from
Task 1 any letter works, since the paths come from the directory listing rather
than a hard-wired preset — but keeping `A` keeps the staged layout comparable
with every other run in this plan.

There is no ground truth, so no Sherd Accuracy: `read_acc` returns `None` and the
verdict comes from `outcome` plus Step 3.

- [ ] **Step 3: Measure the joins**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
PYTHONPATH=. .venv/bin/python -c "
import glob, re, numpy as np
from sfsbench.stage import read_mesh
from sfsbench.scale import _nearest_distances
fs = sorted(glob.glob('/tmp/sfsbench-potB/Dataset/Result/*_Top_1_OBJ_*.obj'),
            key=lambda p: int(re.search(r'OBJ_(\d+)', p).group(1)))
M = [read_mesh(f)[0] for f in fs]
for i in range(len(M)):
    for j in range(i + 1, len(M)):
        d = _nearest_distances(M[i], M[j], samples=4000)
        print(f'{i+1}-{j+1}: min {d.min():8.2f}  p05 {np.percentile(d,5):8.2f}')
"
```

A pair whose minimum is under about 1 mm (or the equivalent in whatever the
units turn out to be, per Step 1) is a claimed join; everything else is separate.

*(The previous revision computed this as `M[i][:,None,:] - M[j][None,::7,:]` at
stride 17/7. On the `potB` result — 615 160 vertices per fragment — that
allocates a 36186 × 5170 × 3 float64 array, **4.5 GB per pair**, six pairs. It
thrashes or raises `MemoryError`. The voxel-grid helper Task 11 needs anyway does
the same job in bounded memory.)*

- [ ] **Step 4: Write the result down**

Record in `docs/superpowers/notes/2026-08-21-potb-result.md`: the physical scale
established in Step 1, the harness verdict, the pairwise minimum distances, the
envelope check, and whether the rimless run (Task 9) said this configuration was
in scope.

State plainly whether the outcome is `no-joins`. Given four fragments of a vessel
over half a metre across, that is a legitimate answer and not a failure of the
software — **provided** Step 1 established the units and the reference sweep
reached its targets. Absent either, the verdict is "not yet interpretable", and
that is what the note must say.

One caveat to record if joins *are* claimed: the harness assumes single-vessel
input. It has no test for "these fragments came from two different pots", though
the dataset and `NUM_MIXED_SHERD` support that case and the published mixed
experiments (Dishes 97.6 %, C+J 80.0 %, E+I 96.6 %, ABCDE 92.5 %, All 87.3 %)
are the harder half of the result. Mixed collections are out of scope here on
cost grounds — Table IX puts the All-mixed run at 42.7 h at (k,b) = (10,5) — but
"no joins" must never be reported when the real answer is "these are two pots".

- [ ] **Step 5: Commit**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp
git add docs/superpowers/notes
git commit -m "Record the potB result under the finished pipeline"
```
