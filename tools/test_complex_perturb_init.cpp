// tools/test_complex_perturb_init.cpp
//
// Runtime construction test for SuperpositionTwoProductInitMps.
// NOT added to CMakeLists.txt. Compile manually on cluster (Intel icpc):
//   icpc -std=c++17 -I${HOME}/.local/include -I../src \
//       -DUSE_COMPLEX=1 -DSYM_LEVEL=1 \
//       -o /tmp/test_perturb test_complex_perturb_init.cpp \
//       (link via CMakeFiles/dmrg.dir/link.txt with the object swapped)
//
// Test data: N=4, a={0,1,0,1}, b={0,0,1,1}
// Swap at sites 1,2: a[1]=1(dn)->0(up), a[2]=0(up)->1(dn)
//
// Run modes:
//   (no args)      — single fixed test: N=4, eps=0.1
//   eps=<value>    — same, with given eps (must be finite and > 0)
//   swap=s1,s2     — choose which sites to swap (must be <4, s1<s2, and
//                    carry opposite spins in the a configuration)
//
// Tolerances (per-quantity, NOT scaled with eps):
//   norm_tol      = 1e-10   absolute, normalization
//   ratio_rtol    = 1e-10   relative, ||r|/eps - 1| <= ratio_rtol
//   phase_tol     = 1e-10   absolute, radians
//   extra_amp_tol = 1e-10   absolute, stray-configuration amplitude
//
// Physical-basis note: label mapping (0=up,1=dn,2=hole) was verified on the
// cluster (qlmps tJ_sites.h: spin_up=0, spin_down=1, empty=2).
//
// Numerically reliable epsilon range: 1e-150 <= eps <= 1e300. Outside this
// range the checks may hit double-representation limits rather than exposing
// construction errors: sqrt(2)*hypot(1,eps) overflows for eps > ~1.27e308,
// and the ratio phase degrades through subnormal noise for eps below ~1e-311.
// All default and planned stress points (0.05, 0.1, 1.0, 1e200) are inside.

#include "qlmps/qlmps.h"
#include "qlten/qlten.h"
#include "tJ_hilbert_space.h"
#include "complex_perturb_init.h"

#include <iostream>
#include <cmath>
#include <complex>
#include <cassert>
#include <vector>
#include <string>
#include <cstdlib>
#include <cerrno>

using namespace qlmps;
using namespace qlten;

static int failures = 0;

void CHECK(bool cond, const char *msg) {
  if (!cond) {
    std::cerr << "FAIL: " << msg << std::endl;
    failures++;
  } else {
    std::cout << "PASS: " << msg << std::endl;
  }
}

void CHECK_CLOSE(double a, double b, double tol, const char *msg) {
  // Reject non-finite inputs BEFORE comparing: NaN comparisons are always
  // false and would otherwise produce a fake PASS.
  if (!std::isfinite(a) || !std::isfinite(b) ||
      !std::isfinite(tol) || tol < 0.0) {
    std::cerr << "FAIL: " << msg << " (non-finite or invalid: a=" << a
              << ", b=" << b << ", tol=" << tol << ")" << std::endl;
    failures++;
    return;
  }
  if (std::abs(a - b) > tol) {
    std::cerr << "FAIL: " << msg << " (got " << a << ", expected " << b
              << ", diff " << std::abs(a - b) << ")" << std::endl;
    failures++;
  } else {
    std::cout << "PASS: " << msg << std::endl;
  }
}

// Compute <config|mps> by contracting the MPS with a product state config.
// config[i] is the physical label at site i.
TenElemT ComputeOverlap(
    FiniteMPS<TenElemT, QNT> &mps,
    const std::vector<size_t> &config,
    size_t N
) {
  // Load all tensors into memory if needed
  for (size_t i = 0; i < N; ++i) {
    mps.LoadTen(kMpsPath, i);
  }

  // Site 0: mps[0] has shape (1, d, D0). Extract physical label config[0].
  auto shape0 = mps[0].GetShape();
  size_t D0 = shape0[2];
  std::vector<TenElemT> vec(D0, TenElemT(0.0));
  for (size_t r = 0; r < D0; ++r) {
    vec[r] = mps[0]({0, config[0], r});
  }

  // Sites 1..N-1: contract with next tensor (right boundary is trivial at end).
  for (size_t i = 1; i < N; ++i) {
    auto shape = mps[i].GetShape();
    size_t Dl = shape[0];
    size_t Dr = shape[2];
    std::vector<TenElemT> new_vec(Dr, TenElemT(0.0));
    for (size_t r = 0; r < Dr; ++r) {
      for (size_t l = 0; l < Dl; ++l) {
        new_vec[r] += vec[l] * mps[i]({l, config[i], r});
      }
    }
    vec = std::move(new_vec);
  }

  assert(vec.size() == 1);
  return vec[0];
}

int main(int argc, char *argv[]) {
  const size_t N = 4;
  double eps = 0.1;
  size_t swap_s1 = 1, swap_s2 = 2;

  // ---- Input validation: eps must be finite and strictly positive ----
  // ---- swap must be two distinct in-range sites with opposite spins ----
  for (int i = 1; i < argc; ++i) {
    std::string arg(argv[i]);
    if (arg.rfind("eps=", 0) == 0) {
      std::string v = arg.substr(4);
      if (v.empty()) {
        std::cerr << "[error] eps value is empty" << std::endl;
        return 2;
      }
      char *end = nullptr;
      errno = 0;
      double parsed = std::strtod(v.c_str(), &end);
      if (end == v.c_str() || *end != '\0') {
        std::cerr << "[error] cannot parse eps: '" << v << "'" << std::endl;
        return 2;
      }
      if (errno == ERANGE) {
        std::cerr << "[error] eps out of range: '" << v << "'" << std::endl;
        return 2;
      }
      if (!std::isfinite(parsed) || parsed <= 0.0) {
        std::cerr << "[error] eps must be finite and strictly positive, got: '"
                  << v << "'" << std::endl;
        return 2;
      }
      eps = parsed;
    } else if (arg.rfind("swap=", 0) == 0) {
      std::string pair = arg.substr(5);
      auto comma = pair.find(',');
      if (comma == std::string::npos) {
        std::cerr << "[error] swap expects 's1,s2': '" << pair << "'" << std::endl;
        return 2;
      }
      std::string v1 = pair.substr(0, comma);
      std::string v2 = pair.substr(comma + 1);
      char *end1 = nullptr, *end2 = nullptr;
      errno = 0;
      unsigned long p1 = std::strtoul(v1.c_str(), &end1, 10);
      if (end1 == v1.c_str() || *end1 != '\0' || errno == ERANGE) {
        std::cerr << "[error] cannot parse swap site 1: '" << v1 << "'" << std::endl;
        return 2;
      }
      errno = 0;
      unsigned long p2 = std::strtoul(v2.c_str(), &end2, 10);
      if (end2 == v2.c_str() || *end2 != '\0' || errno == ERANGE) {
        std::cerr << "[error] cannot parse swap site 2: '" << v2 << "'" << std::endl;
        return 2;
      }
      if (p1 >= N || p2 >= N) {
        std::cerr << "[error] swap sites must be < N=4, got (" << p1 << ","
                  << p2 << ")" << std::endl;
        return 2;
      }
      if (p1 == p2 || p1 > p2) {
        std::cerr << "[error] swap sites must satisfy s1 < s2 and s1 != s2, got ("
                  << p1 << "," << p2 << ")" << std::endl;
        return 2;
      }
      swap_s1 = p1;
      swap_s2 = p2;
    } else {
      std::cerr << "Unknown arg: " << arg << std::endl;
      return 2;
    }
  }

  // ---- Per-quantity tolerances (NOT scaled with eps) ----
  const double norm_tol = 1e-10;
  const double ratio_rtol = 1e-10;
  const double phase_tol = 1e-10;
  const double extra_amp_tol = 1e-10;

  // Setup site vector
  const SiteVec<TenElemT, QNT> sites = tj_site.GenUniformSites<TenElemT>(N);
  FiniteMPS<TenElemT, QNT> mps(sites);

  // Test configurations: swap (s1,s2) on a={0,1,0,1}.
  std::vector<size_t> a_labs = {0, 1, 0, 1};  // up, dn, up, dn
  if (a_labs[swap_s1] == a_labs[swap_s2]) {
    std::cerr << "[error] swap sites (" << swap_s1 << "," << swap_s2
              << ") do not carry opposite spins in a={0,1,0,1}" << std::endl;
    return 2;
  }
  std::vector<size_t> b_labs = a_labs;
  std::swap(b_labs[swap_s1], b_labs[swap_s2]);

  // Normalized amplitudes: |psi> = (|a> + eps_c|b>) / hypot(1, eps)
  const double norm_factor = std::hypot(1.0, eps);
  const double amp_a_val = 1.0 / norm_factor;
  const double amp_b_re = eps / (std::sqrt(2.0) * norm_factor);
  const double amp_b_im = amp_b_re;
  CHECK(std::isfinite(amp_a_val) && std::isfinite(amp_b_re),
        "amplitudes finite (hypot handles huge eps)");

  TenElemT amp_a(amp_a_val);
  TenElemT amp_b(amp_b_re, amp_b_im);

  std::vector<std::vector<size_t>> stat_labs_set = {a_labs, b_labs};

  // --- Construct superposition MPS (4-arg signature with amp_a, amp_b) ---
  SuperpositionTwoProductInitMps(mps, stat_labs_set, amp_a, amp_b);
  mps.Dump(kMpsPath, false);

  std::cout << "=== Section 7.1: Sector Structure ===" << std::endl;

  // Expected right-bond structure of cut k (between sites k-1 and k, 1-based)
  // is driven by the swap pair (s1, s2) and the observed QLTEN sector
  // organization (cluster run, swap=(1,2)):
  //   k <= s1     : the two lanes carry identical QN, QLTEN merges them
  //                 -> 1 sector, dim 2; coordinates 0 and 1 share the QN.
  //   s1 < k <= s2: lanes differ -> 2 sectors, dim 1 each, different QNs.
  //   k > s2      : one remaining sector, dim 1 -> coordinate 1 is out of
  //                 range and must NOT be accessed.
  // Every coordinate access below is gated by the computed total bond
  // dimension; out-of-range coordinates are never touched.
  for (size_t i = 0; i < N - 1; ++i) {
    const size_t k = i + 1;
    auto idx = mps[i].GetIndexes()[2];  // right virtual index
    const size_t num_scts = idx.GetQNSctNum();

    // Total bond dimension = sum of sector dims (summed explicitly so no
    // out-of-range sector API is ever called).
    size_t total_dim = 0;
    for (size_t s = 0; s < num_scts; ++s) {
      total_dim += idx.GetQNSct(s).dim();
    }

    if (k <= swap_s1) {
      std::string m = "Cut " + std::to_string(k) + " (k<=s1): 1 merged sector";
      CHECK(num_scts == 1, m.c_str());
      m = "Cut " + std::to_string(k) + " (k<=s1): bond dim 2";
      CHECK(total_dim == 2, m.c_str());
      if (num_scts == 1) {
        m = "Cut " + std::to_string(k) + " (k<=s1): sector dim 2";
        CHECK(idx.GetQNSct(0).dim() == 2, m.c_str());
      }
      if (total_dim >= 2) {
        auto qn0 = idx.GetQNSctFromActualCoor(0).GetQn();
        auto qn1 = idx.GetQNSctFromActualCoor(1).GetQn();
        m = "Cut " + std::to_string(k) + " (k<=s1): coor 0/1 same QN";
        CHECK(qn0 == qn1, m.c_str());
      }
    } else if (k <= swap_s2) {
      std::string m = "Cut " + std::to_string(k) + " (s1<k<=s2): 2 sectors";
      CHECK(num_scts == 2, m.c_str());
      m = "Cut " + std::to_string(k) + " (s1<k<=s2): bond dim 2";
      CHECK(total_dim == 2, m.c_str());
      if (num_scts == 2) {
        m = "Cut " + std::to_string(k) + " (s1<k<=s2): each sector dim 1";
        CHECK(idx.GetQNSct(0).dim() == 1 && idx.GetQNSct(1).dim() == 1,
              m.c_str());
      }
      if (total_dim >= 2) {
        auto qn0 = idx.GetQNSctFromActualCoor(0).GetQn();
        auto qn1 = idx.GetQNSctFromActualCoor(1).GetQn();
        m = "Cut " + std::to_string(k) + " (s1<k<=s2): coor 0/1 differ in QN";
        CHECK(!(qn0 == qn1), m.c_str());
      }
    } else {
      // k > s2: single sector dim 1; coordinate 1 is out of range here.
      std::string m = "Cut " + std::to_string(k) + " (k>s2): 1 sector";
      CHECK(num_scts == 1, m.c_str());
      m = "Cut " + std::to_string(k) + " (k>s2): bond dim 1";
      CHECK(total_dim == 1, m.c_str());
      if (num_scts == 1) {
        m = "Cut " + std::to_string(k) + " (k>s2): sector dim 1";
        CHECK(idx.GetQNSct(0).dim() == 1, m.c_str());
      }
      if (total_dim >= 1) {
        auto qn0 = idx.GetQNSctFromActualCoor(0).GetQn();
        m = "Cut " + std::to_string(k) + " (k>s2): coor 0 maps to sector QN";
        CHECK(qn0 == idx.GetQNSct(0).GetQn(), m.c_str());
      }
    }
  }

  std::cout << "\n=== Section 7.2: Normalization ===" << std::endl;

  TenElemT overlap_a = ComputeOverlap(mps, a_labs, N);
  TenElemT overlap_b = ComputeOverlap(mps, b_labs, N);

  bool overlaps_finite = std::isfinite(std::abs(overlap_a)) &&
                         std::isfinite(std::abs(overlap_b));
  CHECK(overlaps_finite, "<a|psi> and <b|psi> both finite");
  if (overlaps_finite) {
    double norm_sq = std::norm(overlap_a) + std::norm(overlap_b);
    CHECK_CLOSE(norm_sq, 1.0, norm_tol, "<a|psi>^2 + <b|psi>^2 = 1");
  }

  std::cout << "\n=== Section 7.3: Complex Amplitude Ratio ===" << std::endl;

  if (overlaps_finite) {
    TenElemT ratio = overlap_b / overlap_a;
    double ratio_mag = std::abs(ratio);
    double ratio_phase = std::arg(ratio);

    // True relative check: | |r|/eps - 1 | <= ratio_rtol. eps > 0 is ensured
    // by input validation. (A max(1,eps)-scaled absolute form would degrade
    // to an absolute tolerance for eps < 1 and could let a wrong small
    // amplitude ratio pass.)
    double ratio_relative = ratio_mag / eps;
    CHECK_CLOSE(ratio_relative, 1.0, ratio_rtol,
                "relative amplitude ratio |r|/eps = 1");
    CHECK_CLOSE(ratio_phase, M_PI / 4.0, phase_tol,
                "arg(<b|psi>/<a|psi>) = +pi/4");

    // For tiny eps the inverse ratio can exceed double range; its argument is
    // only meaningful when the value itself is finite. Non-finite is a FAIL,
    // never a silent pass.
    TenElemT inv_ratio = overlap_a / overlap_b;
    bool inv_finite = std::isfinite(std::abs(inv_ratio));
    CHECK(inv_finite, "inverse ratio <a|psi>/<b|psi> is finite");
    if (inv_finite) {
      double inv_phase = std::arg(inv_ratio);
      CHECK_CLOSE(inv_phase, -M_PI / 4.0, phase_tol,
                  "arg(<a|psi>/<b|psi>) = -pi/4");
    }
  } else {
    std::cerr << "SKIPPED ratio checks: overlaps not finite" << std::endl;
  }

  std::cout << "\n=== Section 7.4: No Extra Configurations ===" << std::endl;
  std::cout << "entering 81-configuration scan (any abort below pins the config)" << std::endl;

  // Check all 3^4 = 81 configurations (0,1,2 per site)
  size_t extra_count = 0;
  size_t scanned = 0;
  for (size_t c0 = 0; c0 < 3; ++c0)
  for (size_t c1 = 0; c1 < 3; ++c1)
  for (size_t c2 = 0; c2 < 3; ++c2)
  for (size_t c3 = 0; c3 < 3; ++c3) {
    std::vector<size_t> cfg = {c0, c1, c2, c3};
    std::cout << "config {" << c0 << "," << c1 << "," << c2 << "," << c3
              << "}: ";
    if (cfg == a_labs || cfg == b_labs) {
      std::cout << "expected lane, skip" << std::endl;
      continue;
    }
    TenElemT ov = ComputeOverlap(mps, cfg, N);
    scanned++;
    if (!std::isfinite(std::abs(ov))) {
      // Non-finite must be counted as FAIL, never as zero amplitude.
      std::cerr << "  FAIL: non-finite overlap for config {"
                << c0 << "," << c1 << "," << c2 << "," << c3 << "}: " << ov
                << std::endl;
      failures++;
    } else if (std::abs(ov) > extra_amp_tol) {
      std::cerr << "  Non-zero overlap with config {"
                << c0 << "," << c1 << "," << c2 << "," << c3
                << "}: " << ov << std::endl;
      extra_count++;
    } else {
      std::cout << "zero (ok)" << std::endl;
    }
  }
  std::cout << "scan complete (" << scanned << " non-lane configs checked)" << std::endl;
  CHECK(extra_count == 0, "No extra configuration amplitudes");

  std::cout << "\n=== Section 7.5: Per-Tensor QN Conservation ===" << std::endl;

  // Total QN of the validated physical configuration |a> (label mapping
  // 0=up, 1=dn, 2=hole verified on the cluster). Both lanes share it because
  // the spin swap preserves particle number and Sz. Zero QN is defined by
  // subtraction, matching the production constructor (complex_perturb_init.h).
  auto sites_info = mps.GetSitesInfo();
  const auto &pb_sites = sites_info.sites;
  QNT total_qn = pb_sites[0].GetQNSctFromActualCoor(a_labs[0]).GetQn();
  for (size_t i = 1; i < N; ++i) {
    total_qn += pb_sites[i].GetQNSctFromActualCoor(a_labs[i]).GetQn();
  }
  const QNT zero_qn = total_qn - total_qn;

  // Expected per-site divergence, from the production constructor's QN flow
  // (head right QN = total - p_qn; middle right QN = left - p_qn; tail right
  // boundary = zero QN):
  //   Div(T_0) = Q_total,   Div(T_i) = 0 for 1 <= i < N.
  // CalcDiv is defined as -sum(IN) + sum(OUT); after Centralize(0) the
  // tensor directions are (IN,OUT,IN) at the head, (OUT,OUT,IN) in the
  // middle, (OUT,OUT,OUT) at the tail. Every valid block is checked
  // independently, not merely the first one.
  const int failures_before_75 = failures;

  for (size_t site = 0; site < N; ++site) {
    std::cout << "Site " << site << ":" << std::endl;

    // --- A. Tensor structure ---
    // A default-constructed tensor has no indices (Rank() == 0), so a rank of
    // 3 also proves the tensor was initialized.
    const auto &ten = mps[site];
    const size_t rank = ten.Rank();
    std::cout << "  rank=" << rank << std::endl;

    // Expected directions after Centralize(0), cluster-verified:
    //   head  (site 0)  : (IN,  OUT, IN)
    //   middle sites    : (OUT, OUT, IN)
    //   tail  (site N-1): (OUT, OUT, OUT)
    const auto exp_dir0 = (site == 0) ? qlten::TenIndexDirType::IN
                                      : qlten::TenIndexDirType::OUT;
    const auto exp_dir1 = qlten::TenIndexDirType::OUT;
    const auto exp_dir2 = (site + 1 == N) ? qlten::TenIndexDirType::OUT
                                          : qlten::TenIndexDirType::IN;
    const auto dir_str = [](const qlten::TenIndexDirType d) -> std::string {
      return (d == qlten::TenIndexDirType::IN) ? "IN" : "OUT";
    };
    std::cout << "  directions=(" << dir_str(exp_dir0) << "," << dir_str(exp_dir1)
              << "," << dir_str(exp_dir2) << ") [checked per index below]"
              << std::endl;
    std::cout << "  expected_divergence=" << ((site == 0) ? "total_qn" : "zero_qn")
              << std::endl;
    std::string m = "Site " + std::to_string(site) +
                    ": rank == 3 (initialized, not default-constructed)";
    CHECK(rank == 3, m.c_str());
    if (rank != 3) {
      std::cerr << "  skipping per-block checks: invalid rank" << std::endl;
      continue;
    }

    const auto &indexes = ten.GetIndexes();
    m = "Site " + std::to_string(site) + ": index 0 direction " + dir_str(exp_dir0);
    CHECK(indexes[0].GetDir() == exp_dir0, m.c_str());
    m = "Site " + std::to_string(site) + ": index 1 direction " + dir_str(exp_dir1);
    CHECK(indexes[1].GetDir() == exp_dir1, m.c_str());
    m = "Site " + std::to_string(site) + ": index 2 direction " + dir_str(exp_dir2);
    CHECK(indexes[2].GetDir() == exp_dir2, m.c_str());

    const auto &blk_map = ten.GetBlkSparDataTen().GetBlkIdxDataBlkMap();
    std::cout << "  block_count=" << blk_map.size() << std::endl;
    m = "Site " + std::to_string(site) + ": block map non-empty";
    CHECK(!blk_map.empty(), m.c_str());
    m = "Site " + std::to_string(site) + ": block map size == GetQNBlkNum";
    CHECK(blk_map.size() == ten.GetQNBlkNum(), m.c_str());
    if (blk_map.empty()) {
      std::cerr << "  skipping per-block checks: empty block map" << std::endl;
      continue;
    }

    // --- B/C/D. Block coordinates, per-block divergence, expected value ---
    bool first_valid = false;
    QNT reference_div = zero_qn;
    size_t valid_blocks = 0;
    for (const auto &kv : blk_map) {
      const size_t blk_key = kv.first;
      const auto &data_blk = kv.second;

      const bool coors_ok = (data_blk.blk_coors.size() == rank);
      m = "Site " + std::to_string(site) + " block " + std::to_string(blk_key) +
          ": blk_coors size == rank";
      CHECK(coors_ok, m.c_str());
      m = "Site " + std::to_string(site) + " block " + std::to_string(blk_key) +
          ": shape size == rank";
      CHECK(data_blk.shape.size() == rank, m.c_str());

      bool bounds_ok = coors_ok;
      for (size_t axis = 0; axis < rank && bounds_ok; ++axis) {
        if (data_blk.blk_coors[axis] >= indexes[axis].GetQNSctNum()) {
          bounds_ok = false;
        }
      }
      m = "Site " + std::to_string(site) + " block " + std::to_string(blk_key) +
          ": blk_coors within sector counts";
      CHECK(bounds_ok, m.c_str());

      if (!coors_ok || !bounds_ok) {
        std::cerr << "  skipping CalcDiv for invalid block " << blk_key
                  << std::endl;
        continue;
      }

      const QNT block_div = qlten::CalcDiv(indexes, data_blk.blk_coors);
      valid_blocks++;

      if (!first_valid) {
        reference_div = block_div;
        first_valid = true;
      } else {
        m = "Site " + std::to_string(site) + " block " + std::to_string(blk_key) +
            ": divergence == reference divergence";
        CHECK(block_div == reference_div, m.c_str());
      }

      const QNT expected_div = (site == 0) ? total_qn : zero_qn;
      m = "Site " + std::to_string(site) + " block " + std::to_string(blk_key) +
          ": divergence == expected site divergence";
      CHECK(block_div == expected_div, m.c_str());
    }

    m = "Site " + std::to_string(site) + ": at least one valid block";
    CHECK(valid_blocks > 0, m.c_str());
  }

  CHECK(failures == failures_before_75, "Per-tensor QN conservation");

  std::cout << "\n=== Summary ===" << std::endl;
  std::cout << "params: N=4, eps=" << eps << ", swap=(" << swap_s1 << ","
            << swap_s2 << "), amp_a=" << amp_a_val
            << ", amp_b=(" << amp_b_re << "," << amp_b_im << ")" << std::endl;
  if (failures == 0) {
    std::cout << "All executed tests PASSED." << std::endl;
  } else {
    std::cout << failures << " test(s) FAILED." << std::endl;
  }

  return failures;
}
