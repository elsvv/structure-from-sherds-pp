# Arbitrary Fragment Reassembly — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reassemble a folder of raw scanned pot fragments without per-vessel parameter tuning, verified against the authors' published results for pots A–E.

**Architecture:** A Python harness (`sfsbench`) drives the existing C++ binaries end to end and scores the result against ground truth. It is built first because the diagnosis that decides all later work is its first use. Diagnosis then identifies which preprocessing defect is fatal; that defect is fixed; finally the millimetre thresholds are derived from the data instead of being hand-set.

**Tech Stack:** Python 3.13 (numpy only), C++17, PCL 1.15, CGAL 6.2, Ceres 2.2, CMake, Homebrew.

## Global Constraints

- Three repositories, all on branch `macos-port`:
  - `/Users/vaceslaveliseev/@dev/structure-from-sherds-pp` — reassembly (`build/SfSpp`, `build/SfSpp_ctrl`, `build/extract_axis`). Remote: `origin` = user's fork. **Harness lives here.**
  - `/Users/vaceslaveliseev/@dev/SfSpp_preprocessing` — preprocessing (`build/MeshPreprocessing`, `build/EdgeLineExtraction`). **No remote yet** — commit locally only.
  - `/Users/vaceslaveliseev/@dev/structure-from-sherds` — ICCV 2021 repo, holds `ICCV Data/` (reference dataset + ground truth). **Read-only in this plan.**
- Python is `/opt/homebrew/bin/python3`. Only `numpy` may be imported. `scipy`, `trimesh`, `open3d` are NOT installed and must NOT be added.
- Rebuild after any C++ change: `cmake --build build -j10` in that repository.
- Reference dataset root: `/Users/vaceslaveliseev/@dev/structure-from-sherds/ICCV Data`
- Ground truth: `ICCV Data/GroundTruth/Transformation/Pot_<X>_Piece_<n>_T.txt` — note the piece number here is **not** zero-padded, while mesh filenames are (`Pot_A_Piece_01_Mesh.obj`).
- Success threshold per fragment, from the paper: mean placement error < 20 mm after one global rigid alignment of the whole assembly.
- Environment variables already implemented and available:
  - Preprocessing: `SFSPP_DATASET_ROOT`, `SFSPP_TEMP_ROOT`, `SFSPP_SAMPLING_RADIUS`, `SFSPP_NORMAL_NEIGHBORS`, `SFSPP_RG_NEIGHBORS`, `SFSPP_SMOOTHNESS_DEG`, `SFSPP_MIN_CLUSTER`, `SFSPP_CURVATURE`, `SFSPP_BASE_RADIUS`, `SFSPP_BASE_NORMAL_DEG`, `SFSPP_BASE_ANGLE_DIFF_DEG`, `SFSPP_BASE_ZBIN`, `SFSPP_BASE_MIN_COUNT`, `SFSPP_BASE_MIN_POINTS`, `SFSPP_BASE_STD_RATIO`
  - Reassembly: `SFS_DATA_ROOT`, `SFS_AUTOSAVE`, `SFS_LENGTH_SCALE`, `SFS_AXIS_COST_TOLERANCE`
- Pipeline order is fixed and non-obvious: **MeshPreprocessing → extract_axis → EdgeLineExtraction → SfSpp**. `EdgeLineExtraction` reads the axis while classifying rim segments, so the axis must exist before it runs. Upstream's README has the axis last; that is wrong.
- Result OBJs written by the reassembly are **face-expanded**: vertex `3k+j` is corner `j` of face `k` of the source mesh. Any comparison against a source mesh must expand the source the same way (`V[F.reshape(-1)]`), or correspondence is lost and every measurement is meaningless.
- Commit after every task. Never use `git add -A` in the pp repository — `data/` and `ICCV Data/` are large and gitignored, but scratch output is not.

---

### Task 1: Harness skeleton — dataset staging

**Files:**
- Create: `tools/sfsbench/__init__.py`
- Create: `tools/sfsbench/stage.py`
- Create: `tools/sfsbench/tests/test_stage.py`

**Interfaces:**
- Produces:
  - `read_mesh(path: str) -> tuple[np.ndarray, np.ndarray]` — returns `(vertices Nx3 float64, faces Mx3 int64)`. Accepts `.ply` (binary little-endian, `float x,y,z` + `uchar red,green,blue`, `uchar`-count int faces) and `.obj`.
  - `vertex_normals(xyz: np.ndarray, tri: np.ndarray) -> np.ndarray` — area-weighted, unit length, Nx3.
  - `voxel_subsample(xyz: np.ndarray, target: int) -> np.ndarray` — returns **indices** of one representative vertex per voxel, voxel size chosen so the count lands within [0.75, 1.3] × target.
  - `write_pcd(path: str, xyz: np.ndarray, nrm: np.ndarray) -> None` — ASCII PCD, fields `x y z normal_x normal_y normal_z curvature`.
  - `write_obj(path: str, xyz: np.ndarray, tri: np.ndarray) -> None`
  - `stage_dataset(src_dir: str, dst_root: str, pot_id: str, target_points: int = 2_000_000) -> list[str]` — converts every mesh in `src_dir` into `<dst_root>/Point/Pot_<pot_id>/Pot_<pot_id>_Piece_NN_Point.pcd` and `<dst_root>/Mesh/Pot_<pot_id>/Pot_<pot_id>_Piece_NN_Mesh.obj`, returns the fragment names in order.

**Why the naming matters:** `mesh_processing.cpp:318` does `filePath.stem().string().substr(0, 14)`. `Pot_A_Piece_01_Point` truncated to 14 characters is exactly `Pot_A_Piece_01`. Any other naming silently breaks every downstream path.

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
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
/opt/homebrew/bin/python3 -m pytest sfsbench/tests/test_stage.py -v
```

Expected: FAIL — `ModuleNotFoundError: No module named 'sfsbench.stage'`

If `pytest` is missing, install it: `/opt/homebrew/bin/python3 -m pip install --user pytest`

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


def write_pcd(path, xyz, nrm):
    n = len(xyz)
    data = np.hstack([xyz, nrm, np.zeros((n, 1))]).astype(np.float32)
    with open(path, "w") as f:
        f.write("# .PCD v0.7 - Point Cloud Data file format\n")
        f.write("VERSION 0.7\n")
        f.write("FIELDS x y z normal_x normal_y normal_z curvature\n")
        f.write("SIZE 4 4 4 4 4 4 4\n")
        f.write("TYPE F F F F F F F\n")
        f.write("COUNT 1 1 1 1 1 1 1\n")
        f.write(f"WIDTH {n}\n")
        f.write("HEIGHT 1\n")
        f.write("VIEWPOINT 0 0 0 1 0 0 0\n")
        f.write(f"POINTS {n}\n")
        f.write("DATA ascii\n")
        np.savetxt(f, data, fmt="%.6g")


def write_obj(path, xyz, tri):
    with open(path, "w") as f:
        f.write("# staged by sfsbench\n")
        np.savetxt(f, xyz, fmt="v %.6g %.6g %.6g")
        np.savetxt(f, tri + 1, fmt="f %d %d %d")


def stage_dataset(src_dir, dst_root, pot_id, target_points=2_000_000):
    """Raw meshes -> <dst_root>/{Point,Mesh}/Pot_<id>/Pot_<id>_Piece_NN_*.

    The 14-character name is mandatory: mesh_processing.cpp truncates the file
    stem to 14 characters to recover the fragment id.
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
/opt/homebrew/bin/python3 -m pytest sfsbench/tests/test_stage.py -v
```

Expected: PASS, 5 passed

- [ ] **Step 5: Commit**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp
git add tools/sfsbench
git commit -m "harness: stage raw meshes into the preprocessing layout"
```

---

### Task 2: Scoring against ground truth

**Files:**
- Create: `tools/sfsbench/score.py`
- Create: `tools/sfsbench/tests/test_score.py`

**Interfaces:**
- Consumes: `read_mesh` from `sfsbench.stage`
- Produces:
  - `kabsch(P: np.ndarray, Q: np.ndarray) -> tuple[np.ndarray, np.ndarray]` — rigid transform `(R, t)` mapping `P` onto `Q`, reflection-free.
  - `expand_faces(V: np.ndarray, F: np.ndarray) -> np.ndarray` — `V[F.reshape(-1)]`, the layout the reassembly writes.
  - `score_result(result_dir: str, source_mesh_dir: str, gt_dir: str, pot_id: str, stride: int = 37) -> dict` — returns
    `{"placed": int, "total_scored": int, "correct": int, "per_fragment": {id: mean_error_mm}, "mean_error": float}`.
    `correct` counts fragments whose mean error is < 20 mm.

**Why `expand_faces`:** the reassembly writes each face's three corners as three separate vertices, so a result OBJ has `3 × n_faces` vertices, not `n_vertices`. Comparing raw vertex arrays silently destroys correspondence — during diagnosis this produced a 55 mm "error" on a perfect match.

- [ ] **Step 1: Write the failing test**

```python
# tools/sfsbench/tests/test_score.py
import numpy as np
from sfsbench.score import kabsch, expand_faces

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

def test_expand_faces_matches_reassembly_layout():
    V = np.array([[0., 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]])
    F = np.array([[0, 1, 2], [0, 2, 3]])
    E = expand_faces(V, F)
    assert E.shape == (6, 3)
    assert np.allclose(E[0], V[0])
    assert np.allclose(E[3], V[0])
    assert np.allclose(E[4], V[2])
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
/opt/homebrew/bin/python3 -m pytest sfsbench/tests/test_score.py -v
```

Expected: FAIL — `ModuleNotFoundError: No module named 'sfsbench.score'`

- [ ] **Step 3: Write the implementation**

```python
# tools/sfsbench/score.py
"""Score a reassembly against manually restored ground truth."""
import glob, os, re
import numpy as np

from sfsbench.stage import read_mesh

SUCCESS_MM = 20.0        # the paper's per-fragment success threshold


def kabsch(P, Q):
    """Rigid transform (R, t) with Q ~= R @ P + t. Reflection-free."""
    cp, cq = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - cp).T @ (Q - cq))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, cq - R @ cp


def expand_faces(V, F):
    """Corner-per-vertex layout, as written by SaveResult."""
    return V[F.reshape(-1)]


def score_result(result_dir, source_mesh_dir, gt_dir, pot_id, stride=37):
    files = sorted(glob.glob(os.path.join(result_dir, "*_Top_1_OBJ_*.obj")),
                   key=lambda p: int(re.search(r"OBJ_(\d+)\.obj$", p).group(1)))
    ours, gts, ids = [], [], []
    for path in files:
        i = int(re.search(r"OBJ_(\d+)\.obj$", path).group(1))
        gt_file = os.path.join(gt_dir, f"Pot_{pot_id}_Piece_{i}_T.txt")
        if not os.path.exists(gt_file):
            continue                     # no ground truth for this collection
        Vo, Fo = read_mesh(os.path.join(source_mesh_dir,
                                        f"Pot_{pot_id}_Piece_{i:02d}_Mesh.obj"))
        Vr, _ = read_mesh(path)
        src = expand_faces(Vo, Fo)[::stride]
        res = Vr[::stride]
        n = min(len(src), len(res))
        T = np.loadtxt(gt_file)
        gts.append((T[:3, :3] @ src[:n].T).T + T[:3, 3])
        ours.append(res[:n])
        ids.append(i)

    out = {"placed": len(files), "total_scored": len(ids),
           "correct": 0, "per_fragment": {}, "mean_error": float("nan")}
    if not ids:
        return out

    R, t = kabsch(np.vstack(ours), np.vstack(gts))
    all_err = []
    for i, O, G in zip(ids, ours, gts):
        e = np.linalg.norm(((R @ O.T).T + t) - G, axis=1)
        out["per_fragment"][i] = float(e.mean())
        out["correct"] += int(e.mean() < SUCCESS_MM)
        all_err.append(e)
    out["mean_error"] = float(np.concatenate(all_err).mean())
    return out
```

- [ ] **Step 4: Run test to verify it passes**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
/opt/homebrew/bin/python3 -m pytest sfsbench/tests/test_score.py -v
```

Expected: PASS, 3 passed

- [ ] **Step 5: Verify against a known-good run**

`score_result` must reproduce the reference result measured during diagnosis:
the authors' preprocessed Pot A gives 8/8 with a mean error of about 4.4 mm.

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
/opt/homebrew/bin/python3 -c "
from sfsbench.score import score_result
I='/Users/vaceslaveliseev/@dev/structure-from-sherds/ICCV Data'
r=score_result('/Users/vaceslaveliseev/@dev/SfSpp_preprocessing/DatasetRef/Result',
               I+'/Mesh', I+'/GroundTruth/Transformation', 'A')
print(r['correct'], '/', r['total_scored'], 'mean', round(r['mean_error'],2))
"
```

Expected output: `8 / 8 mean 4.4` (± 0.1)

If this does not reproduce, stop — the scorer is wrong and every later measurement would be worthless.

- [ ] **Step 6: Commit**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp
git add tools/sfsbench/score.py tools/sfsbench/tests/test_score.py
git commit -m "harness: score a reassembly against ground truth"
```

---

### Task 3: Pipeline runner and the reference report

**Files:**
- Create: `tools/sfsbench/run.py`
- Create: `tools/sfsbench/cli.py`
- Create: `tools/sfsbench/tests/test_run.py`

**Interfaces:**
- Consumes: `stage_dataset` (Task 1), `score_result` (Task 2)
- Produces:
  - `PipelinePaths` — dataclass with fields `pp_repo: str`, `prep_repo: str`, `dataset_root: str`, `temp_root: str`, `pot_id: str`
  - `run_pipeline(paths: PipelinePaths, env: dict[str, str] | None = None, binary: str = "SfSpp") -> dict` — runs MeshPreprocessing, extract_axis for every fragment, EdgeLineExtraction, then the named reassembly binary. Returns
    `{"clusters": {name: int}, "breakline_points": {name: int}, "breakline_segments": {name: int}, "axes": {name: int}, "matches_raw": int, "matches_pruned": int, "score": int, "result_dir": str}`
  - `outcome(stats: dict, scored: dict) -> str` — one of `"assembled"`, `"misassembled"`, `"no-joins"`.

**Outcome rules** (the spec requires three, not two):
- `"no-joins"` when `matches_pruned == 0` or fewer than two fragments appear in the result
- `"assembled"` when ground truth exists and `scored["correct"] == scored["total_scored"]`
- `"misassembled"` otherwise

- [ ] **Step 1: Write the failing test**

```python
# tools/sfsbench/tests/test_run.py
from sfsbench.run import outcome, parse_stats

def test_no_joins_when_nothing_survives_pruning():
    assert outcome({"matches_pruned": 0, "placed": 4}, {}) == "no-joins"

def test_no_joins_when_fewer_than_two_placed():
    assert outcome({"matches_pruned": 12, "placed": 1}, {}) == "no-joins"

def test_assembled_when_every_scored_fragment_is_correct():
    assert outcome({"matches_pruned": 46, "placed": 8},
                   {"correct": 8, "total_scored": 8}) == "assembled"

def test_misassembled_when_some_fragment_is_wrong():
    assert outcome({"matches_pruned": 46, "placed": 8},
                   {"correct": 3, "total_scored": 8}) == "misassembled"

def test_misassembled_without_ground_truth():
    assert outcome({"matches_pruned": 46, "placed": 4},
                   {"correct": 0, "total_scored": 0}) == "misassembled"

def test_parse_stats_reads_the_binaries_stdout():
    log = "\n".join([
        "Number of clusters is equal to 12",
        "Total number : 3587",
        "Total number pruned : 459",
        "Score : 205",
    ])
    s = parse_stats(log)
    assert s["matches_raw"] == 3587
    assert s["matches_pruned"] == 459
    assert s["score"] == 205
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
/opt/homebrew/bin/python3 -m pytest sfsbench/tests/test_run.py -v
```

Expected: FAIL — `ModuleNotFoundError: No module named 'sfsbench.run'`

- [ ] **Step 3: Write the implementation**

```python
# tools/sfsbench/run.py
"""Drive the C++ binaries end to end and collect the numbers diagnosis needs."""
import os, re, subprocess
from dataclasses import dataclass

PP_REPO = "/Users/vaceslaveliseev/@dev/structure-from-sherds-pp"
PREP_REPO = "/Users/vaceslaveliseev/@dev/SfSpp_preprocessing"
ICCV = "/Users/vaceslaveliseev/@dev/structure-from-sherds/ICCV Data"


@dataclass
class PipelinePaths:
    pp_repo: str
    prep_repo: str
    dataset_root: str
    temp_root: str
    pot_id: str


def parse_stats(log):
    def one(pattern, cast=int, default=None):
        m = re.search(pattern, log)
        return cast(m.group(1)) if m else default
    return {
        "matches_raw": one(r"Total number : (\d+)"),
        "matches_pruned": one(r"Total number pruned : (\d+)"),
        "score": one(r"Score : (\d+)"),
        "clusters": [int(x) for x in re.findall(r"Number of clusters is equal to (\d+)", log)],
        "breakline_points": [int(x) for x in re.findall(r"(\d+) points, rim", log)],
        "breakline_segments": [int(x) for x in re.findall(r"(\d+) segment\(s\)", log)],
    }


def outcome(stats, scored):
    if not stats.get("matches_pruned") or stats.get("placed", 0) < 2:
        return "no-joins"
    if scored.get("total_scored", 0) and scored["correct"] == scored["total_scored"]:
        return "assembled"
    return "misassembled"


def _run(cmd, cwd, env, log_path):
    merged = dict(os.environ)
    merged.update(env or {})
    with open(log_path, "w") as f:
        subprocess.run(cmd, cwd=cwd, env=merged, stdout=f, stderr=subprocess.STDOUT,
                       check=False)
    return open(log_path, errors="replace").read()


def run_pipeline(paths, env=None, binary="SfSpp"):
    env = dict(env or {})
    prep_env = {k: v for k, v in env.items() if k.startswith("SFSPP_")}
    prep_env["SFSPP_DATASET_ROOT"] = paths.dataset_root
    prep_env["SFSPP_TEMP_ROOT"] = paths.temp_root
    sfs_env = {k: v for k, v in env.items() if k.startswith("SFS_")}
    sfs_env["SFS_DATA_ROOT"] = paths.dataset_root + os.sep
    sfs_env["SFS_AUTOSAVE"] = "exit"

    logs = os.path.join(paths.temp_root, "logs")
    os.makedirs(logs, exist_ok=True)
    prep_build = os.path.join(paths.prep_repo, "build")
    pp_build = os.path.join(paths.pp_repo, "build")
    pot = paths.pot_id

    for sub in ("Axes", f"Breaklines/Pot_{pot}", f"Surfaces/Pot_{pot}", "Result", "Graph Log"):
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
        axes[name] = sum(1 for line in open(axis_path) if line.strip())

    edge_log = _run(["./EdgeLineExtraction"], prep_build, prep_env,
                    os.path.join(logs, "edge.log"))
    sfs_log = _run([f"./{binary}"], pp_build, sfs_env,
                   os.path.join(logs, "reassembly.log"))

    stats = parse_stats(mesh_log + "\n" + edge_log + "\n" + sfs_log)
    stats["axes"] = axes
    result_dir = os.path.join(paths.dataset_root, "Result")
    stats["result_dir"] = result_dir
    stats["placed"] = len([f for f in os.listdir(result_dir) if "_Top_1_OBJ_" in f])
    return stats
```

```python
# tools/sfsbench/cli.py
"""sfsbench — one command from raw meshes to a scored reassembly."""
import argparse, os, sys

from sfsbench.run import PipelinePaths, run_pipeline, outcome, PP_REPO, PREP_REPO, ICCV
from sfsbench.score import score_result
from sfsbench.stage import stage_dataset

REFERENCE_TARGETS = {"A": 8, "B": 9, "C": 4, "D": 22, "E": 19}


def main(argv=None):
    ap = argparse.ArgumentParser(prog="sfsbench")
    ap.add_argument("--src", required=True, help="folder of raw .ply/.obj fragments")
    ap.add_argument("--pot-id", required=True, help="single letter, e.g. A")
    ap.add_argument("--work", required=True, help="working directory for this run")
    ap.add_argument("--gt", default=None, help="ground-truth Transformation folder")
    ap.add_argument("--source-mesh-dir", default=None,
                    help="folder holding the original Pot_X_Piece_NN_Mesh.obj used for scoring")
    ap.add_argument("--target-points", type=int, default=2_000_000)
    ap.add_argument("--binary", default="SfSpp")
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

    paths = PipelinePaths(PP_REPO, PREP_REPO, dataset_root, temp_root, args.pot_id)
    stats = run_pipeline(paths, env, args.binary)

    scored = {}
    if args.gt:
        scored = score_result(stats["result_dir"],
                              args.source_mesh_dir or os.path.join(dataset_root, "Mesh", f"Pot_{args.pot_id}"),
                              args.gt, args.pot_id)

    verdict = outcome(stats, scored)
    print(f"pot {args.pot_id}: {verdict}")
    print(f"  fragments placed   : {stats['placed']}")
    print(f"  matches raw/pruned : {stats['matches_raw']} / {stats['matches_pruned']}")
    print(f"  score              : {stats['score']}")
    print(f"  clusters           : {stats['clusters']}")
    print(f"  breakline points   : {stats['breakline_points']}")
    print(f"  axes per fragment  : {list(stats['axes'].values())}")
    if scored.get("total_scored"):
        target = REFERENCE_TARGETS.get(args.pot_id)
        line = f"  correct            : {scored['correct']} / {scored['total_scored']}"
        if target:
            line += f"   (target {target})"
        print(line)
        print(f"  mean error         : {scored['mean_error']:.2f} mm")
    return 0 if verdict == "assembled" else 1


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: Run test to verify it passes**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
/opt/homebrew/bin/python3 -m pytest sfsbench/tests/test_run.py -v
```

Expected: PASS, 6 passed

- [ ] **Step 5: Reproduce the known 0/8 baseline through the harness**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
I="/Users/vaceslaveliseev/@dev/structure-from-sherds/ICCV Data"
PYTHONPATH=. /opt/homebrew/bin/python3 -m sfsbench.cli \
  --src "$I/Mesh" --pot-id A --work /tmp/sfsbench-A \
  --gt "$I/GroundTruth/Transformation" --source-mesh-dir "$I/Mesh" \
  --binary SfSpp_ctrl
```

Expected: `pot A: misassembled`, `correct : 0 / 8`, matching the diagnosis result.
The harness is correct when it reproduces the *known wrong* answer — that proves
it is driving the same pipeline, not a different one.

Note `--src "$I/Mesh"` picks up all 90 meshes in that folder, not only Pot A.
If the staging step produces more than 8 fragments, copy the eight
`Pot_A_Piece_0*_Mesh.obj` files into a scratch folder first and point `--src` there.

- [ ] **Step 6: Commit**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp
git add tools/sfsbench/run.py tools/sfsbench/cli.py tools/sfsbench/tests/test_run.py
git commit -m "harness: end-to-end runner with a three-way outcome"
```

---

### Task 4: Diagnosis — which defect is fatal

**Files:**
- Create: `tools/sfsbench/ablate.py`
- Create: `docs/superpowers/notes/2026-08-21-ablation-results.md`

**Interfaces:**
- Consumes: `run_pipeline`, `score_result`
- Produces: `substitute(dataset_root: str, pot_id: str, component: str, reference_root: str) -> None` — replaces one component of a staged dataset with the authors' files by symlink. `component` is one of `"axes"`, `"breaklines"`, `"surfaces"`.

This task produces **a decision, not a feature**. Its deliverable is the notes
file recording which run repaired the assembly.

- [ ] **Step 1: Write the substitution helper**

```python
# tools/sfsbench/ablate.py
"""Swap one preprocessing product for the authors' own, to isolate the defect."""
import os

def substitute(dataset_root, pot_id, component, reference_root):
    """Symlink the authors' files over ours for exactly one component.

    reference_root is the flat ICCV Data folder: Axes/, Breaklines/, Surfaces/
    hold Pot_A_Piece_01_* directly, with no per-pot subdirectory.
    """
    plan = {
        "axes":       ("Axes",                       "{n}_Axis.xyz"),
        "breaklines": (f"Breaklines/Pot_{pot_id}",   "{n}_Breakline_0.pcd"),
        "surfaces":   (f"Surfaces/Pot_{pot_id}",     "{n}_Surface_0.xyz"),
    }
    if component not in plan:
        raise ValueError(f"unknown component {component!r}")
    subdir, template = plan[component]
    dst_dir = os.path.join(dataset_root, subdir)
    os.makedirs(dst_dir, exist_ok=True)

    names = [f"Pot_{pot_id}_Piece_{i:02d}" for i in range(1, 9)]
    files = []
    for n in names:
        files.append(template.format(n=n))
        if component == "surfaces":
            files.append(f"{n}_Surface_1.xyz")

    src_dir = os.path.join(reference_root, subdir.split("/")[0])
    for f in files:
        src = os.path.join(src_dir, f)
        dst = os.path.join(dst_dir, f)
        if not os.path.exists(src):
            raise FileNotFoundError(src)
        if os.path.islink(dst) or os.path.exists(dst):
            os.remove(dst)
        os.symlink(src, dst)
```

- [ ] **Step 2: Run ablation A — our surfaces and breaklines, their axes**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
I="/Users/vaceslaveliseev/@dev/structure-from-sherds/ICCV Data"
cp -R /tmp/sfsbench-A /tmp/ablate-axes
PYTHONPATH=. /opt/homebrew/bin/python3 -c "
from sfsbench.ablate import substitute
substitute('/tmp/ablate-axes/Dataset', 'A', 'axes', '$I')
"
```

Then re-run only the reassembly (preprocessing output is already staged):

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/build
rm -rf /tmp/ablate-axes/Dataset/Result && mkdir -p /tmp/ablate-axes/Dataset/Result
SFS_DATA_ROOT=/tmp/ablate-axes/Dataset/ SFS_AUTOSAVE=exit ./SfSpp_ctrl 2>&1 | grep -E "Score|Total number"
cd ../tools
PYTHONPATH=. /opt/homebrew/bin/python3 -c "
from sfsbench.score import score_result
I='$I'
r=score_result('/tmp/ablate-axes/Dataset/Result', I+'/Mesh',
               I+'/GroundTruth/Transformation', 'A')
print('axes substituted ->', r['correct'], '/', r['total_scored'])
"
```

Record the number.

- [ ] **Step 3: Run ablation B — their breaklines**

Same as Step 2 with `/tmp/ablate-breaklines` and `substitute(..., 'breaklines', ...)`.
Record the number.

- [ ] **Step 4: Run ablation C — their surfaces**

Same as Step 2 with `/tmp/ablate-surfaces` and `substitute(..., 'surfaces', ...)`.
Note that `extract_axis` must be re-run after substituting surfaces, because the
axes are derived from them:

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/build
for i in 01 02 03 04 05 06 07 08; do
  ./extract_axis /tmp/ablate-surfaces/Dataset/Surfaces/Pot_A/Pot_A_Piece_${i}_Surface_0.xyz \
                 /tmp/ablate-surfaces/Dataset/Surfaces/Pot_A/Pot_A_Piece_${i}_Surface_1.xyz \
                 /tmp/ablate-surfaces/Dataset/Axes/Pot_A_Piece_${i}_Axis.xyz 1000 > /dev/null 2>&1
done
```

Record the number.

- [ ] **Step 5: Write the decision down**

Create `docs/superpowers/notes/2026-08-21-ablation-results.md` with a table of the
four numbers (baseline 0/8, plus A, B, C) and one sentence naming the component
that repaired the assembly. Tasks 5 and 6 are gated on this result: implement
only the one the ablation blames, and mark the other as not-needed in the notes.

- [ ] **Step 6: Commit**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp
git add tools/sfsbench/ablate.py docs/superpowers/notes
git commit -m "diagnosis: isolate which preprocessing product breaks reassembly"
```

---

### Task 5: Port the reference rim criterion

**Do this task only if Task 4 blamed the breaklines.**

**Files:**
- Modify: `/Users/vaceslaveliseev/@dev/SfSpp_preprocessing/edgeline_extraction.cpp` — add `isBreaklineSegARim_Profile` next to the existing helpers (they sit just above `void processFragmentData`, around line 2650), and call it at the two `isBreaklineSegARim` call sites (currently lines 2927 and 3051)

**Interfaces:**
- Produces: `static bool isBreaklineSegARim_Profile(pcl::PointCloud<pcl::PointXYZ>::Ptr seg, const std::string& axisPath)`

**Criterion, from `AxisExtraction/check_base_and_rim.m`** — a segment is a rim when,
measured against the axis:

- `std(r) < 1.0 mm` — radius near-constant along the segment
- `std(h) < 1.0 mm` — height near-constant, where `h = direction · point`
- `|mean(diff(r)[1:-1])| < 0.1` and `|mean(diff(h)[1:-1])| < 0.1`
- at least 20 points in the segment

All four lengths are millimetres at the reference scale and must go through
`sfsEnvDouble` so Task 7 can rescale them.

- [ ] **Step 1: Write the implementation**

```cpp
// Port of the rim half of AxisExtraction/check_base_and_rim.m. The shipped C++
// heuristic (isBreaklineSegARim / _DistanceFromAxis) finds no rim at all on
// Pot A, where the authors' own files mark five of six fragments; rim pruning is
// what removes roughly a third of the candidate matches, all false.
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

Run `EdgeLineExtraction` on the Pot A working directory from Task 3 and compare
the `info` field of each produced header against the authors' files.

```bash
cd /Users/vaceslaveliseev/@dev/SfSpp_preprocessing/build
SFSPP_DATASET_ROOT=/tmp/sfsbench-A/Dataset SFSPP_TEMP_ROOT=/tmp/sfsbench-A/Temp ./EdgeLineExtraction 2>&1 | grep '\[breakline\]'
I="/Users/vaceslaveliseev/@dev/structure-from-sherds/ICCV Data"
for i in 01 02 03 04 05 06 07 08; do
  echo "Piece_$i ours: $(sed -n 2p /tmp/sfsbench-A/Dataset/Breaklines/Pot_A/Pot_A_Piece_${i}_Breakline_0.pcd)  theirs: $(sed -n 2p "$I/Breaklines/Pot_A_Piece_${i}_Breakline_0.pcd")"
done
```

Expected: the third field of our `#` line matches theirs for at least 6 of 8
fragments. The reference is `2, 1, 1, 1, 1, 1, 0, 0` for pieces 1–8.

- [ ] **Step 4: Re-run the full pipeline and score**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
I="/Users/vaceslaveliseev/@dev/structure-from-sherds/ICCV Data"
PYTHONPATH=. /opt/homebrew/bin/python3 -m sfsbench.cli \
  --src /tmp/potA-src --pot-id A --work /tmp/sfsbench-A2 \
  --gt "$I/GroundTruth/Transformation" --source-mesh-dir "$I/Mesh" \
  --binary SfSpp_ctrl
```

(`/tmp/potA-src` is the eight-mesh scratch folder from Task 3, Step 5.)

Expected: `correct` strictly greater than 0. Record the number.

- [ ] **Step 5: Commit**

```bash
cd /Users/vaceslaveliseev/@dev/SfSpp_preprocessing
git add edgeline_extraction.cpp
git commit -m "Port the reference rim criterion from check_base_and_rim.m

The shipped heuristic finds no rim on Pot A where the authors' own files mark
five of six. Rim pruning removes roughly a third of the candidate matches, all
false, which is why our run starts from 224 matches where the reference starts
from 76."
```

---

### Task 6: Segment the mesh instead of the point cloud

**Do this task only if Task 4 blamed the surfaces and Task 5 did not reach the target.**

**Files:**
- Create: `/Users/vaceslaveliseev/@dev/SfSpp_preprocessing/mesh_segmentation.cpp`
- Modify: `/Users/vaceslaveliseev/@dev/SfSpp_preprocessing/CMakeLists.txt` — add a `MeshSegmentation` executable alongside the existing two

**Rationale:** the authors discard the mesh and grow regions over a point cloud
using PCA normals from 5 neighbours spanning about a millimetre. Measured on our
scans, 46–68% of neighbouring normal pairs exceed the 4.5° growing threshold, so
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
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
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
    std::string tag;
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

`mesh_segmentation.cpp` also needs `#include <array>` and `#include <sstream>`;
add them with the other standard headers.

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
I="/Users/vaceslaveliseev/@dev/structure-from-sherds/ICCV Data"
for i in 01 02 03 04 05 06 07 08; do
  ./build/MeshSegmentation "$I/Mesh/Pot_A_Piece_${i}_Mesh.obj" /tmp/segA_${i}
done
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/build
for i in 01 02 03 04 05 06 07 08; do
  ./extract_axis /tmp/segA_${i}_Surface_0.xyz /tmp/segA_${i}_Surface_1.xyz /tmp/segA_${i}_Axis.xyz 1000 2>/dev/null | tail -1
done
```

Then compare each `/tmp/segA_NN_Axis.xyz` against `$I/Axes/Pot_A_Piece_NN_Axis.xyz`
using the angle measurement from Task 2's scorer.

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

### Task 7: Infer the length scale from the data

**Files:**
- Create: `tools/sfsbench/scale.py`
- Create: `tools/sfsbench/tests/test_scale.py`
- Modify: `tools/sfsbench/cli.py` — set `SFS_LENGTH_SCALE` and the `SFSPP_*` lengths from the inferred values unless the caller overrode them

**Interfaces:**
- Consumes: nothing from earlier tasks except `read_mesh`
- Produces:
  - `wall_thickness(surface0_path: str, surface1_path: str, samples: int = 3000) -> float` — median nearest-neighbour distance from the inner surface to the outer, in mm
  - `infer_scales(surface_pairs: list[tuple[str, str]], extents: list[float]) -> dict` — returns
    `{"thickness_scale": float, "extent_scale": float}`, each relative to the reference constants below

**Reference constants, measured on the authors' Pot A:**

```python
REFERENCE_THICKNESS_MM = 3.69     # median across the 8 fragments
REFERENCE_EXTENT_MM = 119.0       # median bbox diagonal-max across the 8 fragments
```

Two scales, not one, because the thresholds measure different things: those that
live across the fracture (correspondence window, inlier threshold, overlap area)
follow wall thickness, while those along the profile (profile bin size, edge-line
resampling) follow fragment extent. On our own fragments the two ratios are 2.2
and roughly 4–5 — they genuinely disagree, so one knob cannot serve both.

- [ ] **Step 1: Write the failing test**

```python
# tools/sfsbench/tests/test_scale.py
import os, tempfile
import numpy as np
from sfsbench.scale import wall_thickness, infer_scales, REFERENCE_THICKNESS_MM

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

def test_reference_thickness_gives_unit_scale():
    # a collection matching the reference must score exactly 1.0
    out = infer_scales([], [], _thickness_override=REFERENCE_THICKNESS_MM,
                       _extent_override=119.0)
    assert abs(out["thickness_scale"] - 1.0) < 1e-9
    assert abs(out["extent_scale"] - 1.0) < 1e-9

def test_thicker_wall_scales_up():
    out = infer_scales([], [], _thickness_override=2 * REFERENCE_THICKNESS_MM,
                       _extent_override=119.0)
    assert abs(out["thickness_scale"] - 2.0) < 1e-9
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
/opt/homebrew/bin/python3 -m pytest sfsbench/tests/test_scale.py -v
```

Expected: FAIL — `ModuleNotFoundError: No module named 'sfsbench.scale'`

- [ ] **Step 3: Write the implementation**

```python
# tools/sfsbench/scale.py
"""Derive the pipeline's millimetre thresholds from the fragments themselves."""
import numpy as np

# Measured on the authors' Pot A, which every published threshold was tuned for.
REFERENCE_THICKNESS_MM = 3.69
REFERENCE_EXTENT_MM = 119.0


def _nearest_distances(P, Q, cell=None, samples=3000, seed=0):
    """Median-friendly nearest-neighbour distances P -> Q via a voxel hash.

    A KD-tree would be natural but scipy is not available in this environment.
    """
    if cell is None:
        span = float(np.max(Q.max(0) - Q.min(0)))
        cell = max(span / 40.0, 1e-6)
    origin = Q.min(0)
    k = np.floor((Q - origin) / cell).astype(np.int64)
    dims = k.max(0) + 1
    h = k[:, 0] + dims[0] * (k[:, 1] + dims[1] * k[:, 2])
    order = np.argsort(h)
    hs = h[order]

    rng = np.random.default_rng(seed)
    picks = rng.choice(len(P), min(samples, len(P)), replace=False)
    out = []
    for i in picks:
        kk = np.floor((P[i] - origin) / cell).astype(np.int64)
        cand = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
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
            out.append(float(np.linalg.norm(Q[idx] - P[i], axis=1).min()))
    return np.asarray(out)


def wall_thickness(surface0_path, surface1_path, samples=3000):
    A = np.loadtxt(surface0_path)[:, :3]
    B = np.loadtxt(surface1_path)[:, :3]
    d = _nearest_distances(A, B, samples=samples)
    return float(np.median(d)) if len(d) else float("nan")


def infer_scales(surface_pairs, extents,
                 _thickness_override=None, _extent_override=None):
    """Two scale factors relative to the reference pot.

    surface_pairs: [(surface_0.xyz, surface_1.xyz), ...] for the collection
    extents:       largest bounding-box side per fragment, in mm
    """
    if _thickness_override is not None:
        thickness = _thickness_override
    else:
        vals = [wall_thickness(a, b) for a, b in surface_pairs]
        vals = [v for v in vals if np.isfinite(v)]
        thickness = float(np.median(vals)) if vals else REFERENCE_THICKNESS_MM

    extent = _extent_override if _extent_override is not None else (
        float(np.median(extents)) if len(extents) else REFERENCE_EXTENT_MM)

    return {"thickness_scale": thickness / REFERENCE_THICKNESS_MM,
            "extent_scale": extent / REFERENCE_EXTENT_MM}
```

Then in `cli.py`, after `run_pipeline` has produced surfaces but before the
reassembly, set the environment from the inferred scales. Restructure `main` so
staging and preprocessing happen first, then:

```python
from sfsbench.scale import infer_scales

pairs = [(os.path.join(temp_root, "Data", f"Pot_{args.pot_id}", f"{n}_Surface_0.xyz"),
          os.path.join(temp_root, "Data", f"Pot_{args.pot_id}", f"{n}_Surface_1.xyz"))
         for n in names]
scales = infer_scales(pairs, extents)
env.setdefault("SFS_LENGTH_SCALE", f"{scales['thickness_scale']:.4f}")
print(f"  inferred scales    : thickness x{scales['thickness_scale']:.2f}, "
      f"extent x{scales['extent_scale']:.2f}")
```

`extents` is the largest bounding-box side per fragment. `stage_dataset`
returns only names, so compute them here from the staged meshes:

```python
from sfsbench.stage import read_mesh

extents = []
for n in names:
    V, _ = read_mesh(os.path.join(dataset_root, "Mesh", f"Pot_{args.pot_id}", f"{n}_Mesh.obj"))
    extents.append(float(np.max(V.max(0) - V.min(0))))
```

This needs `import numpy as np` at the top of `cli.py`.

- [ ] **Step 4: Run test to verify it passes**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
/opt/homebrew/bin/python3 -m pytest sfsbench/tests/test_scale.py -v
```

Expected: PASS, 3 passed

- [ ] **Step 5: Verify the inferred scale is a no-op on the reference pot**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
I="/Users/vaceslaveliseev/@dev/structure-from-sherds/ICCV Data"
PYTHONPATH=. /opt/homebrew/bin/python3 -c "
from sfsbench.scale import infer_scales
I='$I'
pairs=[(f'{I}/Surfaces/Pot_A_Piece_{i:02d}_Surface_0.xyz',
        f'{I}/Surfaces/Pot_A_Piece_{i:02d}_Surface_1.xyz') for i in range(1,9)]
print(infer_scales(pairs, [119.0]))
"
```

Expected: `thickness_scale` within 0.95–1.05. Anything else means the change
would alter behaviour on the pots that already work — fix before continuing.

- [ ] **Step 6: Commit**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp
git add tools/sfsbench/scale.py tools/sfsbench/tests/test_scale.py tools/sfsbench/cli.py
git commit -m "harness: infer millimetre thresholds from wall thickness and extent"
```

---

### Task 8: The reference sweep

**Files:**
- Create: `tools/sfsbench/sweep.py`
- Modify: `docs/superpowers/notes/2026-08-21-ablation-results.md` — append the sweep table

**Interfaces:**
- Consumes: everything above
- Produces: `sweep(pots: list[str], work_root: str) -> dict[str, dict]` — runs the full chain for each reference pot and returns its scored result

The pots are selected in `class/data_path.h` at compile time, so each needs its
own binary. Add one target per pot in `CMakeLists.txt` following the existing
`SfSpp_ctrl` pattern, defining `POT_A` … `POT_E` respectively, and guard each
preset the way `POT_TEST` is already guarded so no two are ever active together.

- [ ] **Step 1: Add the per-pot binaries**

In `class/data_path.h`, extend the existing guard so every preset is mutually
exclusive:

```cpp
#if !defined(POT_CTRL) && !defined(POT_A) && !defined(POT_B) && \
    !defined(POT_C) && !defined(POT_D) && !defined(POT_E)
#define POT_TEST			// our own scans, see BUILD-macOS.md
#endif
```

In `CMakeLists.txt`, after the `SfSpp_ctrl` block:

```cmake
foreach(pot A B C D E)
  add_executable(SfSpp_${pot} ${PROJ_SRC})
  target_include_directories(SfSpp_${pot} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/class ${PCL_INCLUDE_DIRS})
  target_link_libraries(SfSpp_${pot} PRIVATE Ceres::ceres ${PCL_LIBRARIES})
  target_compile_definitions(SfSpp_${pot} PRIVATE ${PCL_DEFINITIONS} POT_${pot})
  if(OpenMP_CXX_FOUND)
    target_link_libraries(SfSpp_${pot} PRIVATE OpenMP::OpenMP_CXX)
  endif()
  if(NOT MSVC)
    target_compile_options(SfSpp_${pot} PRIVATE -Wno-deprecated-declarations)
  endif()
endforeach()
```

- [ ] **Step 2: Build and verify all five exist**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release > /dev/null && cmake --build build -j10 2>&1 | grep -E "error:|Built target"
ls -1 build/SfSpp_[A-E]
```

Expected: five binaries, no errors.

- [ ] **Step 3: Write the sweep**

```python
# tools/sfsbench/sweep.py
"""Run every reference pot through the whole chain and tabulate the result."""
import os

from sfsbench.cli import REFERENCE_TARGETS
from sfsbench.run import ICCV, PP_REPO, PREP_REPO, PipelinePaths, run_pipeline, outcome
from sfsbench.score import score_result
from sfsbench.stage import stage_dataset

SHERD_COUNT = {"A": 8, "B": 9, "C": 4, "D": 28, "E": 31}


def _pot_source(pot, scratch):
    """The reference Mesh folder is flat and holds every pot; isolate one."""
    import glob, shutil
    out = os.path.join(scratch, f"src_{pot}")
    os.makedirs(out, exist_ok=True)
    for src in sorted(glob.glob(os.path.join(ICCV, "Mesh", f"Pot_{pot}_Piece_*_Mesh.obj"))):
        shutil.copy(src, out)
    return out


def sweep(pots, work_root):
    results = {}
    for pot in pots:
        work = os.path.join(work_root, pot)
        dataset_root = os.path.join(work, "Dataset")
        temp_root = os.path.join(work, "Temp")
        os.makedirs(temp_root, exist_ok=True)
        stage_dataset(_pot_source(pot, work_root), dataset_root, pot)
        paths = PipelinePaths(PP_REPO, PREP_REPO, dataset_root, temp_root, pot)
        stats = run_pipeline(paths, {}, binary=f"SfSpp_{pot}")
        scored = score_result(stats["result_dir"], os.path.join(ICCV, "Mesh"),
                              os.path.join(ICCV, "GroundTruth", "Transformation"), pot)
        results[pot] = {"stats": stats, "scored": scored,
                        "outcome": outcome(stats, scored),
                        "target": REFERENCE_TARGETS.get(pot),
                        "total": SHERD_COUNT[pot]}
    return results


def format_table(results):
    lines = ["| Pot | correct | of | target | outcome |", "|---|---|---|---|---|"]
    for pot, r in results.items():
        lines.append(f"| {pot} | {r['scored'].get('correct', 0)} | {r['total']} | "
                     f"{r['target']} | {r['outcome']} |")
    return "\n".join(lines)
```

- [ ] **Step 4: Run the sweep**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
PYTHONPATH=. /opt/homebrew/bin/python3 -c "
from sfsbench.sweep import sweep, format_table
r = sweep(['A','B','C','D','E'], '/tmp/sfs-sweep')
print(format_table(r))
"
```

Expected: pots A, B and C at their full targets. D and E at or above the paper's
`b=3, k=5` numbers. This is the definition of done from the spec.

Pots D and E take tens of minutes each; run them last and expect the sweep to
occupy the better part of an hour.

- [ ] **Step 5: Append the table to the notes and commit**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp
git add tools/sfsbench/sweep.py CMakeLists.txt class/data_path.h docs/superpowers/notes
git commit -m "harness: reference sweep over pots A-E"
```

---

### Task 9: Run our own fragments with everything in place

**Files:**
- Create: `docs/superpowers/notes/2026-08-21-potb-result.md`

- [ ] **Step 1: Run the full-resolution fragments**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
PYTHONPATH=. /opt/homebrew/bin/python3 -m sfsbench.cli \
  --src ../data/potB --pot-id A --work /tmp/sfsbench-potB
```

There is no ground truth, so the verdict comes from `outcome` plus the render.

- [ ] **Step 2: Measure the joins**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp/tools
PYTHONPATH=. /opt/homebrew/bin/python3 -c "
import glob, re, numpy as np
from sfsbench.stage import read_mesh
fs=sorted(glob.glob('/tmp/sfsbench-potB/Dataset/Result/*_Top_1_OBJ_*.obj'),
          key=lambda p:int(re.search(r'OBJ_(\d+)',p).group(1)))
M=[read_mesh(f)[0][::17] for f in fs]
for i in range(len(M)):
    for j in range(i+1,len(M)):
        d=np.min(np.linalg.norm(M[i][:,None,:]-M[j][None,::7,:],axis=2))
        print(f'{i+1}-{j+1}: {d:8.1f} mm')
"
```

A pair under about 1 mm is a claimed join; everything else is separate.

- [ ] **Step 3: Write the result down**

Record in `docs/superpowers/notes/2026-08-21-potb-result.md`: the harness verdict,
the pairwise distances, the inferred scales, and — if joins are claimed — whether
the edge-on renders show a continuous wall or two shells crossing at an angle.

State plainly whether the outcome is `no-joins`. Given four fragments of a vessel
over half a metre across, that is a legitimate answer and not a failure of the
software.

- [ ] **Step 4: Commit**

```bash
cd /Users/vaceslaveliseev/@dev/structure-from-sherds-pp
git add docs/superpowers/notes
git commit -m "Record the potB result under the finished pipeline"
```
