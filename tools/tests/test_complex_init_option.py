"""
Static source-level checks for the complex initial state perturbation feature.
Verifies that dmrg.cpp, myutil, and complex_perturb_init.h contain the expected
code patterns without requiring compilation.
"""

import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def read_file(relpath):
    path = os.path.join(REPO_ROOT, relpath)
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        return f.read()

def check(condition, msg):
    if not condition:
        print(f"FAIL: {msg}")
        return False
    print(f"PASS: {msg}")
    return True

def main():
    ok = True
    dmrg = read_file("src/dmrg.cpp")
    myutil_h = read_file("src/myutil.h")
    myutil_cpp = read_file("src/myutil.cpp")
    perturb_h = read_file("src/complex_perturb_init.h")

    # --- myutil.h ---
    ok &= check("ParserComplexInitArgs" in myutil_h,
                 "myutil.h declares ParserComplexInitArgs")

    # --- myutil.cpp ---
    ok &= check("ParserComplexInitArgs" in myutil_cpp,
                 "myutil.cpp defines ParserComplexInitArgs")
    ok &= check("--complex-init-eps=" in myutil_cpp,
                 "myutil.cpp parses --complex-init-eps=")
    ok &= check("--complex-init-seed=" in myutil_cpp,
                 "myutil.cpp parses --complex-init-seed=")
    ok &= check("isfinite" in myutil_cpp,
                 "myutil.cpp validates with isfinite")
    ok &= check("ERANGE" in myutil_cpp,
                 "myutil.cpp checks ERANGE for overflow")
    ok &= check("errno" in myutil_cpp,
                 "myutil.cpp uses errno for strtod/strtoul")
    ok &= check("'-'" in myutil_cpp or "== '-'" in myutil_cpp,
                 "myutil.cpp rejects negative seed input")
    ok &= check("return true" in myutil_cpp and "return false" in myutil_cpp,
                 "myutil.cpp returns bool for parse success/failure")
    ok &= check("eps_mag = -1.0" not in myutil_cpp,
                 "myutil.cpp no longer uses -1.0 sentinel")

    # --- complex_perturb_init.h ---
    ok &= check("#ifndef USE_REAL" in perturb_h,
                 "complex_perturb_init.h guarded by #ifndef USE_REAL")
    ok &= check("SuperpositionTwoProductInitMps" in perturb_h,
                 "complex_perturb_init.h defines SuperpositionTwoProductInitMps")
    ok &= check("Centralize(0)" in perturb_h,
                 "complex_perturb_init.h calls Centralize(0)")
    ok &= check("Normalize()" in perturb_h,
                 "complex_perturb_init.h calls Normalize()")
    ok &= check("amp_b" in perturb_h,
                 "complex_perturb_init.h uses amp_b amplitude")

    # --- dmrg.cpp ---
    ok &= check("ParserComplexInitArgs" in dmrg,
                 "dmrg.cpp calls ParserComplexInitArgs")
    ok &= check("!ParserComplexInitArgs" in dmrg,
                 "dmrg.cpp checks return value of parser")
    ok &= check("const bool use_complex_perturb = (complex_init_eps > 0.0)" in dmrg,
                 "dmrg.cpp derives use_complex_perturb from eps value")
    ok &= check("#ifdef USE_REAL" in dmrg,
                 "dmrg.cpp has USE_REAL rejection block")
    ok &= check("--complex-init-eps requires complex build" in dmrg,
                 "dmrg.cpp prints USE_REAL rejection message")
    ok &= check("honeycomb_nn_links" in dmrg,
                 "dmrg.cpp extracts honeycomb NN links")
    ok &= check("nearest_neighbor_links" in dmrg,
                 "dmrg.cpp reads lattice.nearest_neighbor_links")
    ok &= check("only supports honeycomb" in dmrg,
                 "dmrg.cpp rejects non-honeycomb geometry")
    ok &= check("honeycomb_nn_links.empty()" in dmrg,
                 "dmrg.cpp checks for empty NN link table")
    ok &= check("NN link table is empty" in dmrg,
                 "dmrg.cpp has distinct error for empty NN table")
    ok &= check("No valid NN swap link" in dmrg,
                 "dmrg.cpp has distinct error for no valid swap")
    ok &= check("IsPathExist(kMpsPath)" in dmrg,
                 "dmrg.cpp has mps-exists protection")
    ok &= check("valid_swap_links" in dmrg,
                 "dmrg.cpp filters valid swap links")
    ok &= check("mt19937" in dmrg,
                 "dmrg.cpp uses mt19937 for deterministic selection")
    ok &= check("SuperpositionTwoProductInitMps" in dmrg,
                 "dmrg.cpp calls SuperpositionTwoProductInitMps")
    ok &= check("MPI_Bcast(&init_error" in dmrg,
                 "dmrg.cpp broadcasts init_error for MPI safety")
    ok &= check("complex_perturb_init.h" in dmrg,
                 "dmrg.cpp includes complex_perturb_init.h")
    ok &= check("count_before[0] != count_after[0]" in dmrg,
                 "dmrg.cpp verifies QN conservation after swap")
    ok &= check("isfinite(amp_a_val)" in dmrg or "std::isfinite(amp_a_val)" in dmrg,
                 "dmrg.cpp checks amplitude overflow")
    ok &= check("lnk.first != lnk.second" in dmrg,
                 "dmrg.cpp asserts link endpoints are distinct")
    ok &= check("orig_lab_s1" in dmrg and "orig_lab_s2" in dmrg,
                 "dmrg.cpp logs original labels of swap sites")
    ok &= check("amp_b=(" in dmrg,
                 "dmrg.cpp logs complex amplitude components")
    ok &= check("if (use_complex_perturb) {" in dmrg and
                 "MPI_Bcast" in dmrg.split("if (use_complex_perturb) {")[-1],
                 "dmrg.cpp MPI_Bcast is conditional on use_complex_perturb")
    ok &= check("catch (...)" in dmrg,
                 "dmrg.cpp has catch-all exception handler")
    ok &= check("GetNumofMps() != N" in dmrg,
                 "dmrg.cpp verifies Dump file count after construction")
    ok &= check("Manual cleanup" in dmrg,
                 "dmrg.cpp warns about manual cleanup on Dump failure")
    ok &= check("PART2" not in dmrg,
                 "dmrg.cpp has no construction comments remaining")

    # --- Normalization stability and QN checks ---
    ok &= check("std::hypot" in dmrg,
                 "dmrg.cpp uses std::hypot for stable normalization")
    ok &= check("amp_a_val" in dmrg and "amp_b_re" in dmrg,
                 "dmrg.cpp computes pre-normalized amp_a and amp_b")
    ok &= check("amp_a," in perturb_h and "amp_b" in perturb_h,
                 "complex_perturb_init.h accepts both amp_a and amp_b parameters")
    ok &= check("total_qn_a" in perturb_h and "total_qn_b" in perturb_h
                 and "total_qn_a == total_qn_b" in perturb_h,
                 "complex_perturb_init.h checks total QN equality of two product states")
    ok &= check("UINT32_MAX" in myutil_cpp,
                 "myutil.cpp checks seed against UINT32_MAX")

    # --- Structural checks beyond string presence ---

    # Preprocessor guard pairing in complex_perturb_init.h
    ifndef_count = perturb_h.count("#ifndef USE_REAL")
    endif_count = perturb_h.count("#endif // USE_REAL")
    ok &= check(ifndef_count == 1 and endif_count == 1,
                 "complex_perturb_init.h has exactly 1 #ifndef/#endif USE_REAL pair")

    # Preprocessor guard pairing in dmrg.cpp
    dmrg_ifndef = dmrg.count("#ifndef USE_REAL")
    dmrg_endif_real = dmrg.count("#endif // USE_REAL")
    ok &= check(dmrg_ifndef == dmrg_endif_real,
                 f"dmrg.cpp #ifndef USE_REAL ({dmrg_ifndef}) matches #endif USE_REAL ({dmrg_endif_real})")

    # Verify coordinate assertion exists in construction function
    ok &= check("maps to unexpected QN" in perturb_h,
                 "complex_perturb_init.h has runtime coordinate mapping assertion")

    # Verify MPI_Bcast for mps_path_exists (NFS safety)
    ok &= check("mps_path_exists" in dmrg and "MPI_Bcast(&mps_path_exists" in dmrg,
                 "dmrg.cpp broadcasts mps-exists check from master")

    # Verify eps=0 path does not touch perturb_labs or SuperpositionTwoProductInitMps
    # by checking they only appear inside #ifndef USE_REAL blocks
    first_ifndef = dmrg.index("#ifndef USE_REAL")
    ok &= check(dmrg.index("perturb_labs") > first_ifndef,
                 "dmrg.cpp perturb_labs only appears after #ifndef USE_REAL")
    ok &= check(dmrg.index("SuperpositionTwoProductInitMps") > first_ifndef,
                 "dmrg.cpp SuperpositionTwoProductInitMps only appears after #ifndef USE_REAL")

    # Verify Centralize argument is 0 (not 1 or other)
    ok &= check("Centralize(0)" in perturb_h and "Centralize(1)" not in perturb_h,
                 "complex_perturb_init.h Centralize argument is 0, not other value")

    if ok:
        print(f"\nAll checks passed.")
        return 0
    else:
        print(f"\nSome checks FAILED.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
