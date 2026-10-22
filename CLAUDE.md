# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

DMRG ground-state solver for the square-lattice and honeycomb t-J model (with t2/J2 NNN terms), built on two header-only libraries from QuantumLiquids:
- **qlten** (TensorToolkit) — QN-conserving tensors, found via `find_path(qlten)`
- **qlmps** (UltraDMRG) — MPS/MPO/DMRG algorithms, found via `find_path(qlmps)`

All binaries share one `params.json` file (see `docs/parameters.md` for the full reference). The typical workflow is: generate/run DMRG on the cluster (SLURM scripts in repo root) → measure correlations → plot with the MATLAB scripts in `plot/`.

## Build

```bash
mkdir build && cd build
cmake .. -DCMAKE_PREFIX_PATH=<path-to-qlten-and-qlmps>
make -j4
```

Key CMake options (defaults in parentheses):
- `HIGHER_SYMMETRY` (ON): U(1)×U(1) particle-number conservation (`SYM_LEVEL=1`). When OFF (`SYM_LEVEL=0`), the JSON params `mu` and `BondSingletPairPerturbation` become active.
- `RealCode` (OFF): OFF = complex tensors, ON = real. Nonzero `phi` (flux) requires a complex build — the code exits with an error otherwise.
- `QLMPS_USE_GPU` (OFF): single-GPU CUDA/cuTENSOR build, MPI not supported.
- `COMPILE_FIX_CODE` (OFF): builds the `fix_*`/`continue` debugging utilities.
- `QLTEN_TIMING_MODE` / `QLMPS_TIMING_MODE` (and `_MPI_` variants): per-part timing output.

This Windows working copy has no local `build/`; compilation happens on the SUSTech cluster (`/share/home/zhangkun/DMRGtJModel/build`, Intel oneAPI + MKL + MPI, see the repo-root `.sh` SLURM scripts for the environment setup).

## Binaries

Built by `add_executable` in `CMakeLists.txt` (the authoritative source → source-file mapping):

| Binary | Purpose | Notes |
|---|---|---|
| `dmrg` | Traditional DMRG | Builds MPO in memory (`MatReprMPO`) on every run |
| `mpogen` | Dump MPO to `./mpo/` | Must be run before `vmps`/`vmps1FM`/`noised_vmps`, which load the MPO from disk |
| `vmps` | MPO-based DMRG | Loads `./mpo` |
| `noised_vmps` | MPO-based DMRG with noise | `TwoSiteFiniteVMPSWithNoise_`, uses `noise` from JSON |
| `vmps1FM` | Same as vmps, fully-polarized initial state | |
| `measure1` | One-site + two-site correlations | Optional 2nd arg = reference site (default 0) |
| `measureSC_fix_ref_ext` | Singlet/triplet pair correlations vs. a fixed reference bond | `--refbond=i,j` or legacy `--start=i` |
| `measureS` | Total spin (CPU only) | |
| `move_center` | Move MPS center to site 0 | |
| `plot_lattice` | Dump lattice SVG | |

Run DMRG (MPI, master = rank 0):

```bash
mpirun -n 8 ./dmrg params.json --D=100,500,1000
```

`--D=` (comma-separated, no spaces, ascending) overrides `Dmin`/`Dmax`; the program runs one sweep pass per value. In `dmrg`/`vmps`, `MaxLanczIter` is ramped linearly from 3 to the JSON value across the D-list. When no `--D=` is given, Dmin/Dmax/CutOff from JSON are used.

## Runtime Data Layout

All binaries read/write plain directories in the working directory (`qlmps::kMpsPath` = `./mps`, `./mpo`, hidden `./temp`):
- `mps/` — ground-state MPS tensors. If a valid `mps/` exists, `dmrg`/`vmps` load it and continue optimizing instead of re-initializing.
- `mpo/` — MPO tensors dumped by `mpogen`.
- `.temp/` — environment tensors for the finite VMPS algorithm. Can reach TB-scale for wide systems; never delete it mid-optimization.

Measurement programs write `*.json` files (one per observable) into the working directory. Output naming convention for archived data: `docs/DATA_NAMING.md` + `plot/GetSimulationParamsPostfix.m`.

## Architecture

Type backbone (`src/tJ_hilbert_space.h`): `TenElemT` (double vs complex, `USE_REAL`/`USE_COMPLEX`), `QNT` (U1U1QN when `SYM_LEVEL=1`, U1QN otherwise), and the global `tj_site` (`qlmps::tJSite`). Everything is templated on these.

`src/params_case.h` — `CaseParams` parses `params.json` (via `qlmps::CaseParamsParserBasic`); `NormalizeGeometry` maps legacy names (`Cylinder`, `OBC`, `Rotated`, `Torus`) to `Square*` names with a warning. `tJModelParamters` is the compact t/t2/J/J2/phi holder passed to MPO construction.

Lattice classes (`src/lattice_base.h`, `src/squarelattice.h`, `src/honeycomb_lattice.h`): `LatticeBase` holds link lists — `nearest_neighbor_links`, `next_nearest_neighbor_links`, `third_nearest_neighbor_links`, plus `direct_*`/`wind_*` lists splitting links into non-winding vs. winding (for flux `phi`). Site index is column-major: `site = x*Ly + y` (`y = site % Ly`). Honeycomb lattices (`HoneycombXCCylinder`, `HoneycombYCCylinder`) have 2 sites per unit cell, so `N = 2*Lx*Ly` (`ComputeTotalSites` in `src/myutil.cpp`) and A/B sublattice = `(x+y)%2`; their NN bonds are classified into directions A/B/C.

`src/tJmodel.h` — `ConstructtJMPOGenerator<LatticeT>` adds all hopping/exchange terms to a `qlmps::MPOGenerator` from the lattice's link lists. `phi != 0` paths apply `exp(±iφ)` phases on winding links.

The geometry dispatch chain (`if Ly==2 → SquareLadder`, `SquareCylinder`, `SquareOBC`, `SquareRotated`, `SquareTorus`, `HoneycombXC`, `HoneycombYC`) is **duplicated in `dmrg.cpp` and `mpogen.cpp`** (and `measure1.cpp` has its own square-only bond measurement). Adding a geometry means updating every dispatch site.

Main-driver pattern (all `src/*.cpp` `main`s): `MPI_Init` → `CaseParams(params.json)` → build lattice → MPO (dmrg: in-memory; vmps: load `./mpo`) → initialize MPS (Néel order; holes spread evenly — honeycomb uses sublattice-resolved placement in `dmrg.cpp`) → sweep → dump to `./mps`.

`src/myutil.cpp` — `ParserBondDimension` (`--D=` parsing), `ComputeTotalSites`, `IsHoneycombGeometry`, momentum-pair generators (marked "useless for t-J model").

`src/my_measure.h` — wrappers around qlmps measurement primitives (`MeasureOneSiteOp`, `MeasureTwoSiteOpGroup`, `MeasureFourSiteOpGroup`, `MeasureBonds`); measurement mains parallelize by assigning observable groups to MPI ranks round-robin (`task_id++ % mpi_size == rank`).

`src/two_site_update_finite_vmps_mpi_impl_continue.h` — customized two-site VMPS update used by the `continue`/`fix_*` recovery utilities (compiled only with `COMPILE_FIX_CODE`).

## Known Constraints / Gotchas

- Honeycomb geometries ignore `phi` (warning printed, `phi` forced to 0).
- `phi != 0` on square geometries requires `RealCode=OFF` (complex); `dmrg`'s initial-state code also assumes complex for the winding terms.
- `vmps*` binaries hard-exit if `./mpo` is missing.
- `params.json` at repo root is the currently active config (honeycomb 12×4); `src/params.json` is an upstream example. SLURM scripts in the repo root (`submit_dmrg_zk.sh`, `measure_corr.sh`, `measure_super_corr.sh`) contain user-specific paths (`/share/home/zhangkun/...`) and scan hardcoded site/bond lists marked "只需要改这里" (only edit here).
