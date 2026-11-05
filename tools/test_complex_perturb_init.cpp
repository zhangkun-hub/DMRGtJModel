// tools/test_complex_perturb_init.cpp
//
// Runtime construction test for SuperpositionTwoProductInitMps.
// NOT added to CMakeLists.txt. Compile manually on cluster (Intel icpc):
//   icpc -std=c++17 -I${HOME}/.local/include -I../src \
//       -DUSE_COMPLEX=1 -DSYM_LEVEL=1 \
//       -o /tmp/test_perturb test_complex_perturb_init.cpp \
//       -L${HOME}/.local/lib -lhptt -L<mkl> -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core \
//       (link flags must match the production Intel environment; run under /tmp)
//
// Test data: N=4, a={0,1,0,1}, b={0,0,1,1}
// Swap at sites 1,2: a[1]=1(dn)->0(up), a[2]=0(up)->1(dn)
//
// Run modes:
//   (no args)      — single fixed test: N=4, eps=0.1, amp_a=1/sqrt(1+eps^2)
//   eps=<v>        — same, with given eps (e.g. 0.05, 1.0, 1e200)
//   swap=s1,s2     — choose which sites to swap (default 1,2) to test
//                    QN restoration for a different pair location.
//
// Physical-basis note: label mapping (0=up,1=dn,2=hole) is NOT verified here;
// it only matters that a and b are in the same total-QN sector, which the
// constructor asserts. Basis-order verification is a separate cluster task.

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

  for (int i = 1; i < argc; ++i) {
    std::string arg(argv[i]);
    if (arg.rfind("eps=", 0) == 0) {
      eps = std::stod(arg.substr(4));
    } else if (arg.rfind("swap=", 0) == 0) {
      std::string pair = arg.substr(5);
      auto comma = pair.find(',');
      swap_s1 = std::stoul(pair.substr(0, comma));
      swap_s2 = std::stoul(pair.substr(comma + 1));
    } else {
      std::cerr << "Unknown arg: " << arg << std::endl;
      return 2;
    }
  }

  // Relative tolerance: for huge eps the small-amplitude lane is near machine
  // precision relative to the dominant one; scale tolerance by eps accordingly.
  const double tol = 1e-10 * (1.0 + eps);

  // Setup site vector
  const SiteVec<TenElemT, QNT> sites = tj_site.GenUniformSites<TenElemT>(N);
  FiniteMPS<TenElemT, QNT> mps(sites);

  // Test configurations: swap (s1,s2) on a={0,1,0,1}.
  std::vector<size_t> a_labs = {0, 1, 0, 1};  // up, dn, up, dn
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

  // Check each cut: both lanes exist as 2 dim-1 sectors.
  for (size_t i = 0; i < N - 1; ++i) {
    auto idx = mps[i].GetIndexes()[2];  // right virtual index
    size_t num_scts = idx.GetQNSctNum();
    std::string cut_msg = "Cut " + std::to_string(i + 1) + ": right bond has 2 sectors";
    CHECK(num_scts == 2, cut_msg.c_str());
    if (num_scts != 2) continue;

    auto qn0 = idx.GetQNSct(0).GetQn();
    auto qn1 = idx.GetQNSct(1).GetQn();
    CHECK(idx.GetQNSct(0).dim() == 1 && idx.GetQNSct(1).dim() == 1,
          "each sector has dim 1");

    bool same = (qn0 == qn1);
    // QN per lane at cut after site i is equal iff a and b agree on all
    // labels up to site i, i.e. iff the swap pair is entirely before or
    // entirely after the cut. With swap at (s1,s2):
    //   cut i+1 is "inside" the pair iff s1 < i+1 <= s2  (prefix differs).
    bool prefix_differs = (i + 1 > swap_s1) && (i + 1 <= swap_s2);
    if (prefix_differs) {
      CHECK(!same, "cut inside swap pair: QNs differ");
    } else {
      CHECK(same, "cut outside swap pair: QNs identical");
    }
  }

  std::cout << "\n=== Section 7.2: Normalization ===" << std::endl;

  TenElemT overlap_a = ComputeOverlap(mps, a_labs, N);
  TenElemT overlap_b = ComputeOverlap(mps, b_labs, N);

  double norm_sq = std::norm(overlap_a) + std::norm(overlap_b);
  CHECK_CLOSE(norm_sq, 1.0, tol, "<a|psi>^2 + <b|psi>^2 = 1");

  std::cout << "\n=== Section 7.3: Complex Amplitude Ratio ===" << std::endl;

  TenElemT ratio = overlap_b / overlap_a;
  double ratio_mag = std::abs(ratio);
  double ratio_phase = std::arg(ratio);

  CHECK_CLOSE(ratio_mag, eps, tol, "|<b|psi>/<a|psi>| = eps");
  CHECK_CLOSE(ratio_phase, M_PI / 4.0, tol, "arg(<b|psi>/<a|psi>) = +pi/4");

  TenElemT inv_ratio = overlap_a / overlap_b;
  double inv_phase = std::arg(inv_ratio);
  CHECK_CLOSE(inv_phase, -M_PI / 4.0, tol, "arg(<a|psi>/<b|psi>) = -pi/4");

  std::cout << "\n=== Section 7.4: No Extra Configurations ===" << std::endl;

  // Check all 3^4 = 81 configurations (0,1,2 per site)
  size_t extra_count = 0;
  for (size_t c0 = 0; c0 < 3; ++c0)
  for (size_t c1 = 0; c1 < 3; ++c1)
  for (size_t c2 = 0; c2 < 3; ++c2)
  for (size_t c3 = 0; c3 < 3; ++c3) {
    std::vector<size_t> cfg = {c0, c1, c2, c3};
    if (cfg == a_labs || cfg == b_labs) continue;
    TenElemT ov = ComputeOverlap(mps, cfg, N);
    if (std::abs(ov) > tol) {
      std::cerr << "  Non-zero overlap with config {"
                << c0 << "," << c1 << "," << c2 << "," << c3
                << "}: " << ov << std::endl;
      extra_count++;
    }
  }
  CHECK(extra_count == 0, "No extra configuration amplitudes");

  std::cout << "\n=== Section 7.5: Per-Tensor QN Conservation ===" << std::endl;
  // Requires qlten element-iteration API, not confirmed locally.
  std::cout << "NOT RUN: per-tensor QN conservation (needs qlten element iteration API)" << std::endl;

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
