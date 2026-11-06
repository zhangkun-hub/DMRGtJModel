// tools/verify_complex_mps.cpp
//
// Independent, read-only verification tool for DMRG MPS data.
//
// Purpose (Stage 7.3): load an existing MPS from disk (no DMRG run, no MPI),
// verify its integrity, reconstruct the two product states |a> (Neel + holes)
// and |b> (swap of |a> at two sites) with EXACTLY the same algorithm as the
// production initial state (src/dmrg.cpp), and compute:
//   c_a = <a|psi>,  c_b = <b|psi>,  r = c_b / c_a,
//   norm^2 = <psi|psi>.
//
// Compile on the cluster like tools/test_complex_perturb_init.cpp (object
// swap against CMakeFiles/dmrg.dir/link.txt of the SAME build kind):
//   complex build: use the compile command of src/dmrg.cpp from
//                  compile_commands.json (contains -DUSE_COMPLEX=1)
//   real build:    use the compile command from a -DRealCode=ON build dir
//                  (contains -DUSE_REAL=1)
// The source is shared: the scalar type follows tJ_hilbert_space.h, which
// switches TenElemT on USE_REAL. A complex-built binary cannot read a real
// MPS and vice versa (qlten reports an element-type mismatch and aborts).
// The JSON field "scalar_type" reports the COMPILE-TIME element type of this
// program; it is not read from the MPS file.
//
// Numeric range: intermediate rescaling keeps the contraction vector/matrix
// within [1e-150, 1e150] per step, but the final absolute amplitudes are
// recovered by multiplying 10^scale_log; if that factor or the result
// leaves the double range, the computation reports failure (JSON null)
// instead of emitting a wrong number.
//
// Usage:
//   verify_complex_mps params.json --swap=57,58 [--mps-dir=mps] \
//       [--eps=0.05] [--output=verify_complex_mps.json]
//   --swap=s1,s2  (required) sites swapped to build |b>; s1<s2, both
//                 non-hole, carrying opposite spins in |a>.
//   --eps=<v>     (optional) if given, also report the theoretical
//                 INITIAL-state ratio eps*(1+i)/sqrt(2) for reference.
//                 The final (post-DMRG) ratio is a measurement result and is
//                 NOT required to match it.
//
// Read-only guarantees: this tool never writes to --mps-dir, never runs
// DMRG, never modifies params.json. Output is a single JSON file.
//
// Verification boundary (see Stage 7.1-7.2 audit): the production dmrg.cpp
// dumps the initial MPS to ./mps and FiniteDMRG then OVERWRITES it, so the
// initial MPS is not preserved unless the user copies ./mps before the
// optimization finishes. This tool cannot distinguish initial vs optimized
// MPS; it only reports measured values. Do not treat post-DMRG values as an
// initial-state check.

#include "qlmps/qlmps.h"
#include "qlten/qlten.h"
#include "tJ_hilbert_space.h"
#include "params_case.h"
#include "myutil.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <string>
#include <cstdlib>
#include <cerrno>
#include <algorithm>

using namespace qlmps;
using namespace qlten;

// Hermitian conjugate of a scalar. The bra side of <psi|psi> must use the
// complex conjugate, not a plain transpose. Under USE_REAL this is identity.
#ifdef USE_REAL
inline TenElemT ConjElem(const TenElemT &z) { return z; }
#else
inline TenElemT ConjElem(const TenElemT &z) {
  return TenElemT(z.real(), -z.imag());
}
#endif

inline bool ElemFinite(const TenElemT &z) {
  return std::isfinite(std::abs(z));
}

// Real/imag accessors that compile under both USE_REAL and complex builds.
inline double RealPart(const TenElemT &z) {
#ifdef USE_REAL
  return z;
#else
  return z.real();
#endif
}

inline double ImagPart(const TenElemT &z) {
#ifdef USE_REAL
  return 0.0;
#else
  return z.imag();
#endif
}

// Rescale the contraction vector when its magnitude leaves the safe window
// [1e-150, 1e150]. scale_log accumulates the log10 factors to apply to the
// final result (so the reported value stays exact even after rescaling).
void RescaleVec(std::vector<TenElemT> &vec, double &scale_log) {
  double max_abs = 0.0;
  for (const auto &x : vec) {
    max_abs = std::max(max_abs, std::abs(x));
  }
  if (max_abs > 1e150) {
    for (auto &x : vec) { x /= 1e150; }
    scale_log += 150.0;
  } else if (max_abs < 1e-150 && max_abs > 0.0) {
    for (auto &x : vec) { x *= 1e150; }
    scale_log -= 150.0;
  }
}

// Number of mps_ten*.qlten files in the given directory (parameterized copy
// of the pattern in src/myutil.cpp GetNumofMps(), which is fixed to ./mps).
size_t CountMpsFiles(const std::string &mps_path) {
  size_t n = 0;
  for (n = 0; n < 100000; ++n) {
    std::string file = mps_path + "/" + kMpsTenBaseName + std::to_string(n) +
                       "." + kQLTenFileSuffix;
    std::ifstream ifs(file, std::ifstream::binary);
    if (ifs.good()) {
      ifs.close();
    } else {
      break;
    }
  }
  return n;
}

// <config|psi> by direct site-by-site contraction (same algorithm as the
// validated ComputeOverlap in tools/test_complex_perturb_init.cpp, with the
// mps path parameterized and intermediate rescaling added). Complexity is
// O(N * D^2 * d): polynomial in the bond dimension, never enumerates the
// Hilbert space. The result does not depend on the canonical center.
// Returns false (and leaves result unchanged) on any non-finite step.
bool ComputeOverlap(
    FiniteMPS<TenElemT, QNT> &mps,
    const std::vector<size_t> &config,
    size_t N,
    const std::string &mps_path,
    TenElemT &result,
    double &scale_log
) {
  for (size_t i = 0; i < N; ++i) {
    mps.LoadTen(mps_path, i);
  }

  auto shape0 = mps[0].GetShape();
  size_t D0 = shape0[2];
  std::vector<TenElemT> vec(D0, TenElemT(0.0));
  for (size_t r = 0; r < D0; ++r) {
    vec[r] = mps[0]({0, config[0], r});
  }

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
    RescaleVec(vec, scale_log);
  }

  if (vec.size() != 1) return false;
  result = vec[0];
  double mult = std::pow(10.0, scale_log);
  if (!std::isfinite(mult)) return false;
  result *= mult;
  return ElemFinite(result);
}

// <psi|psi> by full transfer contraction. The environment after site i is
// the D x D matrix (flattened row-major):
//   E_i(r, r') = sum_{l,l',p} E_{i-1}(l, l') * conj(A_i[l,p,r]) * A_i[l',p,r'].
// The bra and ket bond indices stay INDEPENDENT until the trivial
// boundaries; a single-bond vector would force the environment diagonal and
// silently drop the cross terms of a superposition state, so it is not
// acceptable. The bra side uses ConjElem (Hermitian conjugate), not a plain
// transpose, and no assumption about the canonical form is made.
// Complexity O(N * D^4 * d) (d = 3). Returns false on any non-finite
// intermediate.
bool ComputeNormSq(
    FiniteMPS<TenElemT, QNT> &mps,
    size_t N,
    const std::string &mps_path,
    TenElemT &norm_sq,
    double &scale_log
) {
  const size_t d = 3;  // t-J local basis: up, dn, hole
  for (size_t i = 0; i < N; ++i) {
    mps.LoadTen(mps_path, i);
  }

  // Head: E_0(r,r') = sum_p conj(A0[0,p,r]) * A0[0,p,r']
  // (the left virtual boundary is trivial, coordinate 0).
  auto shape0 = mps[0].GetShape();
  size_t D0 = shape0[2];
  std::vector<TenElemT> E(D0 * D0, TenElemT(0.0));
  for (size_t r = 0; r < D0; ++r) {
    for (size_t rp = 0; rp < D0; ++rp) {
      TenElemT acc(0.0);
      for (size_t p = 0; p < d; ++p) {
        acc += ConjElem(mps[0]({0, p, r})) * mps[0]({0, p, rp});
      }
      E[r * D0 + rp] = acc;
    }
  }

  // Middle/tail: E is (Dl x Dl) flattened; after site i it becomes
  // (Dr x Dr). Zero entries of E are skipped (exact comparison with 0.0,
  // valid for both real and complex element types).
  for (size_t i = 1; i < N; ++i) {
    auto shape = mps[i].GetShape();
    size_t Dl = shape[0];
    size_t Dr = shape[2];
    std::vector<TenElemT> newE(Dr * Dr, TenElemT(0.0));
    for (size_t r = 0; r < Dr; ++r) {
      for (size_t rp = 0; rp < Dr; ++rp) {
        TenElemT acc(0.0);
        for (size_t l = 0; l < Dl; ++l) {
          for (size_t lp = 0; lp < Dl; ++lp) {
            const TenElemT el = E[l * Dl + lp];
            if (el == TenElemT(0.0)) continue;
            for (size_t p = 0; p < d; ++p) {
              acc += el * ConjElem(mps[i]({l, p, r})) * mps[i]({lp, p, rp});
            }
          }
        }
        newE[r * Dr + rp] = acc;
      }
    }
    E = std::move(newE);
    RescaleVec(E, scale_log);
  }

  // Tail: the right virtual boundary is trivial (dim 1), so E must be 1x1.
  if (E.size() != 1) return false;
  norm_sq = E[0];
  double mult = std::pow(10.0, scale_log);
  if (!std::isfinite(mult)) return false;
  norm_sq *= mult;
  return ElemFinite(norm_sq);
}

// ---- JSON output helpers (all values finite or null; never NaN/Inf) ----
void WriteReal(std::ofstream &ofs, double v) {
  if (std::isfinite(v)) {
    ofs << std::setprecision(17) << v;
  } else {
    ofs << "null";
  }
}

void WriteComplexObj(
    std::ofstream &ofs,
    const std::string &indent,
    const TenElemT &z,
    bool valid
) {
  ofs << indent << "{\n";
  if (valid) {
    double re = RealPart(z), im = ImagPart(z);
    ofs << indent << "  \"real\": ";
    WriteReal(ofs, re);
    ofs << ",\n" << indent << "  \"imag\": ";
    WriteReal(ofs, im);
    ofs << ",\n" << indent << "  \"abs\": ";
    WriteReal(ofs, std::abs(z));
    ofs << ",\n" << indent << "  \"phase\": ";
    WriteReal(ofs, std::arg(z));
    ofs << "\n";
  } else {
    ofs << indent << "  \"real\": null,\n";
    ofs << indent << "  \"imag\": null,\n";
    ofs << indent << "  \"abs\": null,\n";
    ofs << indent << "  \"phase\": null\n";
  }
  ofs << indent << "}";
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cerr << "usage: " << argv[0]
              << " params.json --swap=s1,s2 [--mps-dir=mps]"
                 " [--eps=<v>] [--output=verify_complex_mps.json]"
              << std::endl;
    return 2;
  }

  // ---- Option parsing (params.json is argv[1]) ----
  std::string mps_dir = "mps";
  std::string output_path = "verify_complex_mps.json";
  bool has_swap = false;
  size_t swap_s1 = 0, swap_s2 = 0;
  bool has_eps = false;
  double eps = 0.0;

  for (int i = 2; i < argc; ++i) {
    std::string arg(argv[i]);
    if (arg.rfind("--mps-dir=", 0) == 0) {
      mps_dir = arg.substr(10);
      if (mps_dir.empty()) {
        std::cerr << "[error] empty --mps-dir" << std::endl;
        return 2;
      }
    } else if (arg.rfind("--output=", 0) == 0) {
      output_path = arg.substr(9);
      if (output_path.empty()) {
        std::cerr << "[error] empty --output" << std::endl;
        return 2;
      }
    } else if (arg.rfind("--eps=", 0) == 0) {
      std::string v = arg.substr(6);
      char *end = nullptr;
      errno = 0;
      double parsed = std::strtod(v.c_str(), &end);
      if (end == v.c_str() || *end != '\0' || errno == ERANGE ||
          !std::isfinite(parsed) || parsed <= 0.0) {
        std::cerr << "[error] invalid --eps value: '" << v << "'" << std::endl;
        return 2;
      }
      eps = parsed;
      has_eps = true;
    } else if (arg.rfind("--swap=", 0) == 0) {
      std::string pair = arg.substr(7);
      auto comma = pair.find(',');
      if (comma == std::string::npos) {
        std::cerr << "[error] --swap expects 's1,s2': '" << pair << "'" << std::endl;
        return 2;
      }
      const std::string v1 = pair.substr(0, comma);
      const std::string v2 = pair.substr(comma + 1);
      char *end1 = nullptr, *end2 = nullptr;
      errno = 0;
      unsigned long p1 = std::strtoul(v1.c_str(), &end1, 10);
      if (end1 == v1.c_str() || *end1 != '\0' || errno == ERANGE) {
        std::cerr << "[error] cannot parse --swap site 1: '" << v1 << "'" << std::endl;
        return 2;
      }
      errno = 0;
      unsigned long p2 = std::strtoul(v2.c_str(), &end2, 10);
      if (end2 == v2.c_str() || *end2 != '\0' || errno == ERANGE) {
        std::cerr << "[error] cannot parse --swap site 2: '" << v2 << "'" << std::endl;
        return 2;
      }
      swap_s1 = p1;
      swap_s2 = p2;
      has_swap = true;
    } else {
      std::cerr << "[error] unknown argument: " << arg << std::endl;
      return 2;
    }
  }

  if (!has_swap) {
    std::cerr << "[error] --swap=s1,s2 is required (explicit swap avoids"
                 " depending on RNG / link ordering)" << std::endl;
    return 2;
  }

  // ---- Geometry and product-state reconstruction ----
  // Exact replica of the production initial-state algorithm (src/dmrg.cpp
  // :167-205), so |a> matches the DMRG initial state by construction.
  CaseParams params(argv[1]);
  const size_t N = ComputeTotalSites(params);
  const size_t Lx = params.Lx, Ly = params.Ly;

  if (swap_s1 >= N || swap_s2 >= N || swap_s1 == swap_s2) {
    std::cerr << "[error] swap sites must satisfy 0 <= s1 < s2 < N=" << N
              << std::endl;
    return 2;
  }

  std::vector<size_t> a_labs(N);
  std::vector<size_t> A_sublattice, B_sublattice;
  for (size_t i = 0; i < N; ++i) {
    const size_t y = i % Ly;
    const size_t x = i / Ly;
    const bool a_sublattice = ((x + y) % 2 == 0);
    if (a_sublattice) {
      A_sublattice.push_back(i);
      a_labs[i] = 0;
    } else {
      B_sublattice.push_back(i);
      a_labs[i] = 1;
    }
  }

  if (params.Numhole > 1) {
    size_t nhalf = params.Numhole / 2;
    size_t stepA = A_sublattice.size() / nhalf;
    size_t stepB = B_sublattice.size() / nhalf;
    size_t offsetA = stepA / 2;
    size_t offsetB = 0;
    for (size_t k = 0; k < nhalf; ++k) {
      size_t idxA = (offsetA + k * stepA) % A_sublattice.size();
      size_t idxB = (offsetB + k * stepB) % B_sublattice.size();
      a_labs[A_sublattice[idxA]] = 2;
      a_labs[B_sublattice[idxB]] = 2;
    }
    if (2 * nhalf < params.Numhole) {
      a_labs[B_sublattice[B_sublattice.size() - 1]] = 2;
    }
  } else if (params.Numhole == 1) {
    a_labs[B_sublattice[B_sublattice.size() / 2]] = 2;
  }

  // Swap validity: both sites non-hole, opposite spins (same rule as the
  // production valid_swap_links filter, src/dmrg.cpp:257-261).
  if (a_labs[swap_s1] == 2 || a_labs[swap_s2] == 2 ||
      a_labs[swap_s1] == a_labs[swap_s2]) {
    std::cerr << "[error] swap sites (" << swap_s1 << "," << swap_s2
              << ") must be non-hole with opposite spins; labels are ("
              << a_labs[swap_s1] << "," << a_labs[swap_s2] << ")" << std::endl;
    return 2;
  }
  std::vector<size_t> b_labs = a_labs;
  std::swap(b_labs[swap_s1], b_labs[swap_s2]);

  // Label-count conservation (same check as src/dmrg.cpp:284-300).
  size_t count_before[3] = {0, 0, 0};
  size_t count_after[3] = {0, 0, 0};
  for (size_t idx = 0; idx < N; ++idx) {
    count_before[a_labs[idx]]++;
    count_after[b_labs[idx]]++;
  }
  bool counts_conserved =
      (count_before[0] == count_after[0]) &&
      (count_before[1] == count_after[1]) &&
      (count_before[2] == count_after[2]);

  std::vector<size_t> hole_sites;
  for (size_t i = 0; i < N; ++i) {
    if (a_labs[i] == 2) hole_sites.push_back(i);
  }

  // Total QN check from the physical basis (label mapping 0=up, 1=dn, 2=hole
  // verified on the cluster, qlmps tJ_sites.h).
  const SiteVec<TenElemT, QNT> sites = tj_site.GenUniformSites<TenElemT>(N);
  const auto &pb_sites = sites.sites;
  QNT total_qn_a = pb_sites[0].GetQNSctFromActualCoor(a_labs[0]).GetQn();
  for (size_t i = 1; i < N; ++i) {
    total_qn_a += pb_sites[i].GetQNSctFromActualCoor(a_labs[i]).GetQn();
  }
  QNT total_qn_b = pb_sites[0].GetQNSctFromActualCoor(b_labs[0]).GetQn();
  for (size_t i = 1; i < N; ++i) {
    total_qn_b += pb_sites[i].GetQNSctFromActualCoor(b_labs[i]).GetQn();
  }
  const bool qn_conserved = (total_qn_a == total_qn_b);

  // ---- MPS loading (read-only) ----
  FiniteMPS<TenElemT, QNT> mps(sites);
  bool load_ok = true;
  std::string load_error = "";
  const size_t file_count = CountMpsFiles(mps_dir);
  if (file_count != N) {
    load_ok = false;
    load_error = "mps file count " + std::to_string(file_count) +
                 " != expected site count " + std::to_string(N) +
                 " in directory '" + mps_dir + "'";
  } else {
    // LoadTen: a scalar-type mismatch (e.g. real binary reading a complex
    // MPS) is detected by qlten, which prints the reason and aborts the
    // process; it cannot be caught here. No implicit conversion exists.
    for (size_t i = 0; i < N; ++i) {
      mps.LoadTen(mps_dir, i);
    }
    for (size_t i = 0; i < N && load_ok; ++i) {
      if (mps[i].Rank() != 3) {
        load_ok = false;
        load_error = "tensor " + std::to_string(i) + " has rank " +
                     std::to_string(mps[i].Rank()) + " != 3";
      }
    }
  }

  // ---- Overlap, ratio, norm (only when loading succeeded) ----
  bool ov_a_ok = false, ov_b_ok = false, norm_ok = false;
  TenElemT overlap_a(0.0), overlap_b(0.0);
  TenElemT norm_sq(0.0);
  double scale_log_a = 0.0, scale_log_b = 0.0, scale_log_norm = 0.0;

  if (load_ok) {
    ov_a_ok = ComputeOverlap(mps, a_labs, N, mps_dir, overlap_a, scale_log_a);
    ov_b_ok = ComputeOverlap(mps, b_labs, N, mps_dir, overlap_b, scale_log_b);
    norm_ok = ComputeNormSq(mps, N, mps_dir, norm_sq, scale_log_norm);
  }

  double norm_abs = 0.0;
  bool norm_positive = false;
  double norm_sq_re = 0.0, norm_sq_im = 0.0;
  if (norm_ok) {
    norm_sq_re = RealPart(norm_sq);
    norm_sq_im = ImagPart(norm_sq);
    // Positivity and hermiticity sanity: <psi|psi> must be real and >= 0.
    const double imag_scale = std::max(1.0, std::abs(norm_sq_re));
    if (std::abs(norm_sq_im) > 1e-8 * imag_scale) {
      std::cout << "[warning] <psi|psi> has non-negligible imaginary part: "
                << norm_sq_im << " (real part " << norm_sq_re << ")"
                << std::endl;
    }
    if (norm_sq_re >= 0.0) {
      norm_positive = true;
      norm_abs = std::sqrt(norm_sq_re);
    }
  }

  // Ratio stability. Threshold rationale: |c_a| is bounded by ||psi||, so
  // |c_a|/||psi|| <= 1; the relative error of c_b/c_a grows like machine
  // epsilon / (|c_a|/||psi||). To keep the relative error of the ratio
  // well below 1e-8 (double noise floor is ~1e-16), require the relative
  // amplitude |c_a|/||psi|| >= 1e-8. This is an ENGINEERING GUARD against
  // amplifying noise, not a rigorous error bound; the threshold is reported
  // in the JSON so every call site is auditable.
  const double rel_amp_threshold = 1e-8;
  double ratio_threshold = 0.0;
  bool ratio_reliable = false;
  std::string ratio_reason = "";
  TenElemT ratio(0.0);
  if (ov_a_ok && ov_b_ok) {
    double scale = (norm_positive && norm_abs > 0.0) ? norm_abs : 1.0;
    ratio_threshold = rel_amp_threshold * scale;
    double amp_a = std::abs(overlap_a);
    if (!std::isfinite(amp_a)) {
      ratio_reason = "|c_a| non-finite";
    } else if (amp_a < ratio_threshold) {
      ratio_reliable = false;
      ratio_reason = "|c_a| = " + std::to_string(amp_a) +
                     " below threshold " + std::to_string(ratio_threshold) +
                     " (relative to ||psi||)";
    } else {
      ratio = overlap_b / overlap_a;
      if (ElemFinite(ratio)) {
        ratio_reliable = true;
      } else {
        ratio_reason = "c_b/c_a non-finite";
      }
    }
  } else if (!ov_a_ok) {
    ratio_reason = "overlap_a computation failed";
  } else {
    ratio_reason = "overlap_b computation failed";
  }

  // ---- JSON output ----
  std::ofstream ofs(output_path);
  if (!ofs.good()) {
    std::cerr << "[error] cannot open output file: " << output_path << std::endl;
    return 2;
  }
  const std::string ind1 = "  ";
  const std::string ind2 = "    ";
  const std::string ind3 = "      ";

#ifdef USE_REAL
  const char *scalar_type = "real";
#else
  const char *scalar_type = "complex";
#endif

  ofs << "{\n";
  ofs << ind1 << "\"mps_load_success\": " << (load_ok ? "true" : "false") << ",\n";
  ofs << ind1 << "\"mps_load_error\": ";
  if (load_ok) { ofs << "null"; } else { ofs << "\"" << load_error << "\""; }
  ofs << ",\n";
  ofs << ind1 << "\"scalar_type\": \"" << scalar_type << "\",\n";
  ofs << ind1 << "\"site_count\": " << N << ",\n";
  ofs << ind1 << "\"mps_file_count\": " << file_count << ",\n";
  ofs << ind1 << "\"mps_dir\": \"" << mps_dir << "\",\n";
  ofs << ind1 << "\"params_geometry\": \"" << params.Geometry << "\",\n";
  ofs << ind1 << "\"params_lx\": " << Lx << ",\n";
  ofs << ind1 << "\"params_ly\": " << Ly << ",\n";
  ofs << ind1 << "\"params_numhole\": " << params.Numhole << ",\n";

  ofs << ind1 << "\"product_state\": {\n";
  ofs << ind2 << "\"hole_sites\": [";
  for (size_t i = 0; i < hole_sites.size(); ++i) {
    ofs << hole_sites[i] << (i + 1 < hole_sites.size() ? ", " : "");
  }
  ofs << "],\n";
  ofs << ind2 << "\"swap_sites\": [" << swap_s1 << ", " << swap_s2 << "],\n";
  ofs << ind2 << "\"swap_original_labels\": [" << a_labs[swap_s1] << ", "
      << a_labs[swap_s2] << "],\n";
  ofs << ind2 << "\"label_counts_conserved\": "
      << (counts_conserved ? "true" : "false") << ",\n";
  ofs << ind2 << "\"total_qn_conserved\": "
      << (qn_conserved ? "true" : "false") << "\n";
  ofs << ind1 << "},\n";

  ofs << ind1 << "\"norm_squared\": ";
  if (norm_ok) {
    ofs << "{\"real\": ";
    WriteReal(ofs, norm_sq_re);
    ofs << ", \"imag\": ";
    WriteReal(ofs, norm_sq_im);
    ofs << "}";
  } else {
    ofs << "null";
  }
  ofs << ",\n";
  ofs << ind1 << "\"norm\": ";
  if (norm_positive) { WriteReal(ofs, norm_abs); } else { ofs << "null"; }
  ofs << ",\n";
  ofs << ind1 << "\"normalization_error\": ";
  if (norm_positive) { WriteReal(ofs, std::abs(norm_abs - 1.0)); }
  else { ofs << "null"; }
  ofs << ",\n";
  // Element-level iteration API for a full all-finite sweep has not been
  // confirmed on the cluster qlten; per the Stage 7.3 spec, report null
  // instead of pretending the check ran.
  ofs << ind1 << "\"all_finite\": null,\n";
  ofs << ind1 << "\"all_finite_note\": "
      << "\"element-level iteration API not confirmed; set to null\",\n";
  ofs << ind1 << "\"rescale_log\": {\n";
  ofs << ind2 << "\"overlap_a\": ";
  WriteReal(ofs, scale_log_a);
  ofs << ",\n" << ind2 << "\"overlap_b\": ";
  WriteReal(ofs, scale_log_b);
  ofs << ",\n" << ind2 << "\"norm\": ";
  WriteReal(ofs, scale_log_norm);
  ofs << "\n" << ind1 << "},\n";

  ofs << ind1 << "\"overlap_a\": ";
  if (ov_a_ok) {
    ofs << "\n";
    WriteComplexObj(ofs, ind2, overlap_a, true);
  } else {
    ofs << "null";
  }
  ofs << ",\n";
  ofs << ind1 << "\"overlap_b\": ";
  if (ov_b_ok) {
    ofs << "\n";
    WriteComplexObj(ofs, ind2, overlap_b, true);
  } else {
    ofs << "null";
  }
  ofs << ",\n";

  ofs << ind1 << "\"relative_ratio\": {\n";
  ofs << ind2 << "\"real\": ";
  if (ratio_reliable) { WriteReal(ofs, RealPart(ratio)); } else { ofs << "null"; }
  ofs << ",\n" << ind2 << "\"imag\": ";
  if (ratio_reliable) { WriteReal(ofs, ImagPart(ratio)); } else { ofs << "null"; }
  ofs << ",\n" << ind2 << "\"abs\": ";
  if (ratio_reliable) { WriteReal(ofs, std::abs(ratio)); } else { ofs << "null"; }
  ofs << ",\n" << ind2 << "\"phase\": ";
  if (ratio_reliable) { WriteReal(ofs, std::arg(ratio)); } else { ofs << "null"; }
  ofs << ",\n" << ind2 << "\"reliable\": "
      << (ratio_reliable ? "true" : "false") << ",\n";
  ofs << ind2 << "\"unreliable_reason\": ";
  if (ratio_reason.empty()) { ofs << "null"; }
  else { ofs << "\"" << ratio_reason << "\""; }
  ofs << ",\n" << ind2 << "\"threshold\": ";
  WriteReal(ofs, ratio_threshold);
  ofs << "\n" << ind1 << "},\n";

  ofs << ind1 << "\"expected_init_ratio\": ";
  if (has_eps) {
    // Theoretical INITIAL-state ratio eps*(1+i)/sqrt(2); for reference only.
    // The post-DMRG ratio is a measurement and need not match.
    double e_re = eps / std::sqrt(2.0);
#ifdef USE_REAL
    TenElemT expected(e_re);  // real build: imaginary part is meaningless
#else
    TenElemT expected(e_re, e_re);
#endif
    ofs << "\n";
    WriteComplexObj(ofs, ind2, expected, true);
  } else {
    ofs << "null";
  }
  ofs << "\n";

  ofs << "}\n";
  ofs.close();

  std::cout << "verification report written to " << output_path << std::endl;
  return load_ok ? 0 : 1;
}
