// ---------------------------------------------------------------------------
// extract_axis - axis of symmetry for one fragment, from its two segmented
// surfaces.
//
// The published preprocessing pipeline (github.com/DominicoRyu/SfSpp_preprocessing)
// ships this stage as MATLAB only, and its `refine_axis.m` is missing from the
// repository altogether, so the chain cannot be completed without MATLAB. This
// is a port of AxisExtraction/run_potsac.m built on the C++ PotSAC that already
// lives in this repository (class/axis_estimation.cpp), so no MATLAB is needed.
//
// It reproduces the *modified* PotSAC of the paper (Sec. 3.2): rather than
// keeping only the RANSAC winner, it keeps the ten best candidates, refines
// them, discards near-duplicates (axes within 10 degrees) and those whose cost
// is more than 10% above the best, and emits every surviving axis. Fragments
// whose axis is ambiguous - small ones, or ones cut nearly parallel to the axis
// - then carry several hypotheses into matching instead of one wrong answer.
//
//   extract_axis <Surface_0.xyz> <Surface_1.xyz> <out_Axis.xyz> [iters]
//
// Output: one line per axis, "px py pz nx ny nz" - the point first, then the
// direction, which is the column order BreakLine::ReadAxis expects.
// ---------------------------------------------------------------------------
#include "data_structure.h"
#include "axis_estimation.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

namespace {

constexpr int    kSampleSize        = 6;      // Pottmann needs >= 6 point-normal pairs
constexpr int    kNumCandidates     = 10;     // run_potsac.m keeps the top 10
constexpr double kDuplicateCosAngle = 0.9848; // cos(10 deg)

// run_potsac.m also discards any refined axis costing more than 10% above the
// best. The paper (Sec. 3.2) describes pruning by angle alone, and the tighter
// filter is what collapses the deliberately-ambiguous fragments back to a
// single (sometimes wrong) axis - Pot A sherd 7 in the paper's Fig. 6 is
// exactly such a case. Tunable, and reported so the choice is visible.
double costTolerance() {
    if (const char* v = std::getenv("SFS_AXIS_COST_TOLERANCE")) {
        try { return std::stod(v); } catch (...) {}
    }
    return 1.10;
}
constexpr int    kRefineIters       = 300;
constexpr int    kNumThreads        = 16;
constexpr int    kStride            = 10;     // run_potsac.m works on C(:,1:10:end)

struct Axis {
    Eigen::Vector3d point;
    Eigen::Vector3d direction;
    double          cost = 0.0;
};

// Index space is the concatenation [inner surface, outer surface], matching
// ComputePottmannAxis in class/axis_estimation.cpp.
void pointNormalAt(const Geom& g, int idx, Eigen::Vector3d& p, Eigen::Vector3d& n) {
    const int num_inner = static_cast<int>(g.sur_in_.point_.cols());
    if (idx < num_inner) {
        p = g.sur_in_.point_.col(idx);
        n = g.sur_in_.normal_.col(idx);
    } else {
        p = g.sur_out_.point_.col(idx - num_inner);
        n = g.sur_out_.normal_.col(idx - num_inner);
    }
}

int totalPoints(const Geom& g) {
    return static_cast<int>(g.sur_in_.point_.cols() + g.sur_out_.point_.cols());
}

// Port of compute_pottmann_axis.m, with the conditioning compute_axis_of_symmetry.m
// applies around it: the sample is centred and scaled before the linear solve and
// the result is mapped back. Without it the normal equations are badly scaled -
// our fragments live at coordinates in the hundreds of millimetres.
bool pottmannAxis(const Geom& g,
                  const std::vector<int>& indices,
                  Eigen::Vector3d& axis_point,
                  Eigen::Vector3d& axis_direction) {
    const int n = static_cast<int>(indices.size());
    if (n < 3) return false;

    Eigen::Matrix3Xd pts(3, n), nrm(3, n);
    for (int i = 0; i < n; ++i) {
        Eigen::Vector3d p, v;
        pointNormalAt(g, indices[i], p, v);
        pts.col(i) = p;
        nrm.col(i) = v.normalized();
    }

    const Eigen::Vector3d centre = pts.rowwise().mean();
    pts.colwise() -= centre;
    const double scale =
        (pts.cwiseAbs().sum() + nrm.cwiseAbs().sum()) / static_cast<double>(6 * n);
    if (!(scale > 0.0) || !std::isfinite(scale)) return false;
    pts /= scale;

    // J = [moment; direction] per sample, as in euc2plucker
    Eigen::Matrix3Xd jt_normal(3, n), jt_point(3, n);
    for (int i = 0; i < n; ++i) {
        jt_normal.col(i) = pts.col(i).cross(nrm.col(i));
        jt_point.col(i)  = nrm.col(i);
    }

    const Eigen::Matrix3d M11 = jt_normal * jt_normal.transpose();
    const Eigen::Matrix3d M22 = jt_point * jt_point.transpose();
    const Eigen::Matrix3d M12 = jt_normal * jt_point.transpose();

    Eigen::FullPivLU<Eigen::Matrix3d> lu(M22);
    if (!lu.isInvertible()) return false;
    const Eigen::Matrix3d M22_inv = lu.inverse();
    const Eigen::Matrix3d M = M11 - M12 * M22_inv * M12.transpose();
    if (!M.allFinite()) return false;

    // Smallest singular vector; Eigen orders singular values descending.
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(M, Eigen::ComputeFullV);
    Eigen::Vector3d v = svd.matrixV().col(2);
    if (v.norm() < 1e-12) return false;
    v.normalize();

    const Eigen::Vector3d v_bar = -(M22_inv * (M12.transpose() * v));
    Eigen::Vector3d t = v.cross(v_bar);

    t = t * scale + centre;
    t -= v.dot(t) * v;             // keep the offset orthogonal to the direction
    if (!t.allFinite()) return false;

    axis_point     = t;
    axis_direction = v;
    return true;
}

// The residual the C++ pipeline optimises (class/axis_estimation.h). Note it is
// the author's own variant of Cao's error: a norm rather than a squared norm,
// flagged "different from PotSAC" in that header.
double caoResidual(const Geom& g, int idx,
                   const Eigen::Vector3d& point, const Eigen::Vector3d& direction) {
    Eigen::Vector3d p, n;
    pointNormalAt(g, idx, p, n);
    BiaxialCaoError err(p, n);
    double residual = 0.0;
    err(point.data(), direction.data(), &residual);
    return std::isfinite(residual) ? residual : std::numeric_limits<double>::max();
}

// Geman-McClure, as compute_axis_of_symmetry.m scores RANSAC hypotheses.
double gemanMcClureCost(const Geom& g, const Eigen::Vector3d& p, const Eigen::Vector3d& d) {
    const int n = totalPoints(g);
    double cost = 0.0;
    for (int i = 0; i < n; ++i) {
        const double r = caoResidual(g, i, p, d);
        const double s = r * r;
        cost += 2.0 * s / (s + 1.0);
    }
    return cost;
}

// Huber, as run_potsac.m ranks the refined candidates.
double huberCost(const Geom& g, const Eigen::Vector3d& p, const Eigen::Vector3d& d) {
    const int n = totalPoints(g);
    double cost = 0.0;
    for (int i = 0; i < n; ++i) {
        const double r = caoResidual(g, i, p, d);
        const double s = r * r;
        cost += (s <= 1.0) ? s : (2.0 * std::sqrt(s) - 1.0);
    }
    return cost;
}

// Every kStride-th point of each surface. Keeping the two surfaces separate
// matters: ComputePottmannAxis and RefineAxis both index them independently.
Geom subsample(const Geom& g, int stride) {
    Geom out;
    auto take = [stride](const MatrixXd& src_p, const MatrixXd& src_n,
                         MatrixXd& dst_p, MatrixXd& dst_n) {
        const int n = static_cast<int>(src_p.cols());
        const int m = (n + stride - 1) / stride;
        dst_p.resize(3, m);
        dst_n.resize(3, m);
        for (int i = 0, j = 0; i < n && j < m; i += stride, ++j) {
            dst_p.col(j) = src_p.col(i);
            dst_n.col(j) = src_n.col(i);
        }
    };
    take(g.sur_in_.point_, g.sur_in_.normal_, out.sur_in_.point_, out.sur_in_.normal_);
    take(g.sur_out_.point_, g.sur_out_.normal_, out.sur_out_.point_, out.sur_out_.normal_);
    return out;
}

// Drop axes within 10 degrees of an earlier, cheaper one, and optionally
// anything more than `costTolerance()` above the best remaining cost.
//
// run_potsac.m applies the cost filter twice: once on the coarsely refined
// candidates and again after refining the survivors on every point. The second
// pass is unsound, because the final refinement reorders the costs - a
// candidate that was cheapest on the subsample need not stay cheapest - so it
// re-derives its budget from a list whose ranking has changed and discards an
// axis that had already been selected. Pot A sherd 7 is precisely that case:
// the correct axis wins the coarse round, loses the final one by a factor of
// 1.136, and is thrown away, leaving only the wrong axis. The paper (Sec. 3.2)
// describes pruning by angle alone, so the final pass here does that.
std::vector<Axis> pruneAxes(std::vector<Axis> axes, bool apply_cost_filter) {
    if (axes.empty()) return axes;
    std::stable_sort(axes.begin(), axes.end(),
                     [](const Axis& a, const Axis& b) { return a.cost < b.cost; });

    std::vector<Axis> unique;
    for (const Axis& a : axes) {
        const bool duplicate =
            std::any_of(unique.begin(), unique.end(), [&a](const Axis& kept) {
                return std::abs(kept.direction.dot(a.direction)) > kDuplicateCosAngle;
            });
        if (!duplicate) unique.push_back(a);
    }

    if (apply_cost_filter) {
        const double budget = unique.front().cost * costTolerance();
        unique.erase(std::remove_if(unique.begin(), unique.end(),
                                    [budget](const Axis& a) { return a.cost > budget; }),
                     unique.end());
    }
    return unique;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "usage: extract_axis <Surface_0.xyz> <Surface_1.xyz> "
                     "<out_Axis.xyz> [ransac_iters]\n";
        return 2;
    }
    const std::string inner_path = argv[1];
    const std::string outer_path = argv[2];
    const std::string out_path   = argv[3];
    const int ransac_iters = (argc > 4) ? std::atoi(argv[4]) : 1000;

    Geom geom;
    LoadSurface(&geom, inner_path, outer_path);
    const int num_points = totalPoints(geom);
    if (num_points < kSampleSize) {
        std::cerr << "error: only " << num_points << " surface points, need at least "
                  << kSampleSize << "\n";
        return 1;
    }
    std::cout << "surfaces: " << geom.sur_in_.point_.cols() << " inner + "
              << geom.sur_out_.point_.cols() << " outer = " << num_points << " points\n";

    const Geom coarse = subsample(geom, kStride);
    const int coarse_points = totalPoints(coarse);
    std::cout << "scoring on every " << kStride << "th point (" << coarse_points << ")\n";

    // --- RANSAC over the coarse cloud -------------------------------------
    std::mt19937 gen(20250206u);   // fixed seed: reruns must reproduce
    std::uniform_int_distribution<int> pick_point(0, coarse_points - 1);

    std::vector<Axis> candidates;
    candidates.reserve(ransac_iters);
    for (int iter = 0; iter < ransac_iters; ++iter) {
        std::vector<int> sample;
        sample.reserve(kSampleSize);
        while (static_cast<int>(sample.size()) < kSampleSize) {
            const int idx = pick_point(gen);
            if (std::find(sample.begin(), sample.end(), idx) == sample.end())
                sample.push_back(idx);
        }
        Axis a;
        if (!pottmannAxis(coarse, sample, a.point, a.direction)) continue;
        a.cost = gemanMcClureCost(coarse, a.point, a.direction);
        if (std::isfinite(a.cost)) candidates.push_back(a);
    }
    if (candidates.empty()) {
        std::cerr << "error: no valid axis hypothesis in " << ransac_iters << " iterations\n";
        return 1;
    }
    std::stable_sort(candidates.begin(), candidates.end(),
                     [](const Axis& a, const Axis& b) { return a.cost < b.cost; });
    candidates.resize(std::min<size_t>(candidates.size(), kNumCandidates));
    std::cout << "RANSAC: " << candidates.size() << " candidates, best GM cost "
              << candidates.front().cost << "\n";

    // --- refine on the coarse cloud, prune --------------------------------
    for (Axis& a : candidates) {
        RefineAxis(&const_cast<Geom&>(coarse), a.point, a.direction,
                   kRefineIters, kNumThreads, 1.0, true);
        a.direction.normalize();
        a.cost = huberCost(coarse, a.point, a.direction);
    }
    std::vector<Axis> axes = pruneAxes(std::move(candidates), /*apply_cost_filter=*/true);
    std::cout << "after coarse refinement: " << axes.size() << " distinct axes (cost tolerance "
              << costTolerance() << ")\n";
    for (size_t i = 0; i < axes.size(); ++i)
        std::cout << "    candidate " << i << " huber cost " << axes[i].cost
                  << "  dir [" << axes[i].direction.transpose() << "]\n";

    // --- refine survivors on every point, prune again ---------------------
    for (Axis& a : axes) {
        RefineAxis(&geom, a.point, a.direction, kRefineIters, kNumThreads, 1.0, true);
        a.direction.normalize();
        a.cost = huberCost(geom, a.point, a.direction);
    }
    axes = pruneAxes(std::move(axes), /*apply_cost_filter=*/false);

    std::ofstream out(out_path);
    if (!out) {
        std::cerr << "error: cannot write " << out_path << "\n";
        return 1;
    }
    out.precision(9);
    for (const Axis& a : axes) {
        out << a.point.x()     << ' ' << a.point.y()     << ' ' << a.point.z()     << ' '
            << a.direction.x() << ' ' << a.direction.y() << ' ' << a.direction.z() << '\n';
    }
    std::cout << "wrote " << axes.size() << " axis/axes to " << out_path << "\n";
    for (size_t i = 0; i < axes.size(); ++i) {
        std::cout << "  axis " << i << ": dir [" << axes[i].direction.transpose()
                  << "]  point [" << axes[i].point.transpose()
                  << "]  huber cost " << axes[i].cost << "\n";
    }
    return 0;
}
