# Tutorial

Before following this tutorial, we recommend you first follow [`docs/quick_run.md`](quick_run.md), which demonstrates how to perform a DMRG calculation for a 2-leg $t$-$J$ model and outlines the general workflow. This tutorial covers more advanced scenarios and provides guidance for more challenging cases.

## Usage of the DMRG Sweep Program

### Argument List

The programs `dmrg` and `vmps` share a unified argument list. The first argument should be a JSON file containing all parameter information, as exemplified by [`src/params.json`](../src/params.json). For a detailed explanation of each parameter, see [`docs/parameters.md`](parameters.md).

The second argument is optional and defines a series of bond dimensions (kept states) for DMRG, in the form:
```
--D=10,30,100,300,1000
```
This argument must start with `--D=`, and the bond dimensions should be separated by commas without spaces. The order matters, as the program will perform sweeps in the order provided. If this argument is given, it overrides the `Dmin` and `Dmax` values in the JSON file. If not provided, the bond dimension settings from the JSON file are used.

### MPI Parallelization

On CPU architectures, you can use any number of MPI processes supported by your machine. If you use only one process, the non-parallel DMRG/VMPS code will run. For multi-process MPI runs, the first process acts as the master and the others as workers. We recommend using either 1 or at least 3 MPI processes for efficiency.

### Thread Number

Before launching calculations, check your machine's architecture. 
On typical CPUs, set the thread number to be less than or equal to the number of physical cores,
as using more threads than cores can not improve performance.
We also suggest the thread number should less than or equal to the number of physical cores in one socket,
as using more cores in one process usually cannot improve performance linearly.

### GPU Program

The program supports only single-GPU execution. MPI parallelization is not available for GPU runs.  
To build the GPU code, set the CUDA architectures in [`CMakeLists.txt`](../CMakeLists.txt):
```cmake
set(CMAKE_CUDA_ARCHITECTURES 70 80) 
```

### Output

As shown in `quick_run.md`, the DMRG sweep output includes lines like:
```
Site (   9,  10) E0 = -8.6334507691833000 TruncErr = 7.50e-04 D =   300 Iter =  21 LanczT =     0.14 TotT =     0.15 S =  3.2029150
```
This indicates that the two-site update between sites 9 and 10 has completed, with:
- Ground state energy $E_0 = -8.6334507691833000$
- Truncation error from SVD $\epsilon = 7.50 \times 10^{-4}$
- Bond dimension $D = 300$
- Lanczos iteration count: 21
- Lanczos time: 0.14 s
- Total update time: 0.15 s
- Entanglement entropy: 3.2029150

This information helps you judge convergence. Good convergence is indicated by stable energy, truncation error, and entanglement entropy.

### Workflow & Data Structure

For wider systems, the truncation error can be large (e.g., $7.5 \times 10^{-4}$). Increase the bond dimension to reduce the truncation error. If a sweep does not converge for the current bond dimension, run additional sweeps or increase the bond dimension.

The wave function data is saved in the `mps` directory. Sequential DMRG runs will load the existing MPS data and continue optimization. The typical workflow is to repeat steps 1–3 from `quick_run.md`, checking convergence after each sweep. If not converged, run more sweeps or increase the bond dimension.

A hidden directory `.temp` is also created, storing environment tensors for the variational MPS algorithm. 
This can require significant disk space, especially for wide systems (e.g., $L_y=8$ may use terabytes). 
Do not delete `.temp` during optimization. If deleted, the program will be broken and is not easy to relaunch.
If the `.temp` path is deleted between two running of the optimization,
tt will be regenerated automatically in the second optimization, but this takes extra time.
If disk space allows, keep `.temp` until all DMRG sweeps are complete.

## Measuring One-site and Two-site correlations.

The first argument for measurement programs should be the same JSON file used for DMRG. 
The second argument is optional and sets the reference site for two-site correlation functions.
If omitted, the reference site defaults to 0.

---

## Measuring Superconducting Pair Correlations

The program `measureSC_fix_ref` is designed to measure superconducting (SC) pair correlation functions in a memory-efficient and parallelized manner.
It computes singlet pair-pair correlation functions 
between a fixed reference bond and some other target bonds in the system.
The reference bond is $y$ bond and target bonds include 
both along the $x$ and $y$ bonds.

### Usage

```bash
mpirun -n 8 ./measureSC_fix_ref params.json [--start=START_SITE] 
```
- `params.json`: unified parameter JSON file in dmrg/vmps/measure1.
- `--start=END_SITE`: (Optional) Sets the upper site index for the reference bond. If not provided, default choose 0.

### Output files and  What Is Measured

1. **Output files**
   
The results are saved as JSON files with names as:
- `scsyya_fix_referX`, `scsyyb_fix_referX`, `scsyyc_fix_referX`, `scsyyd_fix_referX` (for $y$-direction target bonds)
- `scsyxa_fix_referX`, `scsyxb_fix_referX`, `scsyxc_fix_referX`, `scsyxd_fix_referX` (for $x$-direction target bonds)

Here, `X` is the reference site index.

2. **What Is Measured**

For a reference bond (typically a nearest-neighbor pair), the program measures the expectation value:
$$
\langle \psi | \Delta^\dagger_{i,j} \Delta_{k,l} | \psi \rangle
$$
where $\Delta_{i,j}$ is a singlet or triplet pair creation operator between sites $i$ and $j$, and $(k,l)$ runs over all other bonds in the system.

The operators measured include:
- **Singlet pairing:**  
  $\Delta_{i,j}^\dagger = \frac{1}{\sqrt{2}} \left( c^\dagger_{i,\uparrow} c^\dagger_{j,\downarrow} - c^\dagger_{i,\downarrow} c^\dagger_{j,\uparrow} \right)$

- **Triplet pairing:**  
  $\Delta_{i,j}^\dagger = \frac{1}{\sqrt{2}} \left( c^\dagger_{i,\uparrow} c^\dagger_{j,\downarrow} + c^\dagger_{i,\downarrow} c^\dagger_{j,\uparrow} \right)$

3. **How the result files are related to SC pairing correlations?**
 
Why for each direction of target bonds, it dumps four files? We denote these four files correlation data as
$A, B, C, $ and $D$. For each file, it measure the following correlaiton
 $$
     \langle \psi | c^\dagger_{i,\sigma} c^\dagger_{j,-\sigma} c_{i',\sigma'} c_{j',-\sigma }    | \psi \rangle
    $$
They are related to the singlet and triplet pair correlations by combination:

     - **Singlet:** $(A - B - C + D)/2$
     - **Triplet:** $(A + B + C + D)/2$


## Advanced: read the source files
How to find the source file for specific binary? Check `add_executable` in CMakeList.txt.




