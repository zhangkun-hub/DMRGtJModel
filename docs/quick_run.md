# Quick Run (MPI on CPU Architecture)

This guide demonstrates how to perform a simple calculation with the Square Lattice t-J Model DMRG code after building the project.  
We show how to perform a DMRG ground state calculation on a $2\times 8$ square lattice with 2 holes doped away from half-filling.  
**Before starting, ensure you have compiled the code with particle number U(1) symmetry conservation.**

## 0. Create a New Directory for the Calculation

We recommend creating a separate directory for each calculation to keep parameter files, wave function data, and results organized.

## 1. Set Up the Parameter File

Create a JSON file (e.g., `quick_params.json`) with the following content:

```json
{
  "CaseParams": {
    "Geometry": "OBC",
    "Lx": 8,
    "Ly": 2,
    "t": 1.0,
    "t2": 0,
    "J": 0.3,
    "J2": 0,
    "phi": 0,
    "mu": 0,
    "NumHole": 2,
    "Sweeps": 4,
    "Dmin": 10,
    "Dmax": 300,
    "CutOff": 1e-8,
    "LanczErr": 1e-9,
    "MaxLanczIter": 70,
    "Threads": 4,
    "noise": [],
    "Perturbation": 0,
    "BondSingletPairPerturbation": 0,
    "wavelength": 0
  }
}
```

This example sets up a $2\times 8$ lattice for the t-J model with $t=1$, $J=0.3$, and no next-nearest-neighbor terms ($t_2=0$, $J_2=0$).

For a detailed explanation of each parameter, see [`docs/parameters.md`](parameters.md).

## 2. Run DMRG Sweep

Run the following command (replace `path_to_build` with your build directory):

```bash
mpirun -n 8 path_to_build/dmrg quick_params.json --D=10,30,100,300
```

- Replace `8` with the number of MPI processes you wish to use.

You should see output similar to:

```
 ****** tJV Model Parameter List ******
NN hopping t  = 1
NNN hopping t2 = 0
NN super-exchange J  = 0.3
NNN super-exchangeJ2 = 0
Twist angle phi = 0
lattice construct
lattice construct
lattice construct
lattice construct
lattice construct
lattice construct
lattice construct
lattice construct
[	9	18	26	34	42	50	58	66	74	82	74	66	58	50	42	34	26	18	9	1]
[	9	18	26	34	42	50	58	66	74	82	74	66	58	50	42	34	26	18	9	1]
[	9	18	26	34	42	50	58	66	74	82	74	66	58	50	42	34	26	18	9	1]
[	9	18	26	34	42	50	58	66	74	82	74	66	58	50	42	34	26	18	9	1]
[	9	18	26	34	42	50	58	66	74	82	74	66	58	50	42	34	26	18	9	1]
[	9	18	26	34	42	50	58	66	74	82	74	66	58	50	42	34	26	18	9	1]
[	9	18	26	34	42	50	58	66	74	82	74	66	58	50	42	34	26	18	9	1]
[	9	18	26	34	42	50	58	66	74	82	74	66	58	50	42	34	26	18	9	1]
Initial mps as direct product state.

=====> MPI Sweep Parameters For DMRG <=====
MPS/MPO size: 	 20
Sweep times: 	 4
Bond dimension: 	 10/300
Truncation error: 	 1e-08
Lanczos max iterations 	70
MPS path: 	mps
Temp path: 	.temp
=====> Technical Parameters <=====
The number of processors(including master): 	8
The number of threads per processor: 	4
=====> Checking and Updating Boundary Tensors =====>
left boundary: 	5
right boundary: 	14
No temp directory .temp
=====> Generating the block operators =====>
right block operator number : 9
right block operator number : 18
right block operator number : 26
right block operator number : 34
right block operator number : 42
right block operator number : 50
right block operator number : 58
right block operator number : 66
right block operator number : 74
right block operator number : 82
right block operator number : 74
right block operator number : 66
right block operator number : 58

sweep 1
```

and output end with
```
sweep 4
Site (   5,   6) E0 = -7.5991299880596825 TruncErr = 9.85e-09 D =   222 Iter =  38 LanczT =     0.07 TotT =     0.07 S =  2.1186915
Site (   6,   7) E0 = -7.5991822221020540 TruncErr = 3.20e-08 D =   300 Iter =  25 LanczT =     0.06 TotT =     0.07 S =  2.3854260
Site (   7,   8) E0 = -7.5991864037555894 TruncErr = 4.75e-08 D =   300 Iter =  14 LanczT =     0.04 TotT =     0.05 S =  2.2765103
Site (   8,   9) E0 = -7.5991873036324513 TruncErr = 3.55e-08 D =   300 Iter =  11 LanczT =     0.03 TotT =     0.03 S =  2.3855835
Site (   9,  10) E0 = -7.5991876238416314 TruncErr = 9.93e-09 D =   238 Iter =   7 LanczT =     0.02 TotT =     0.02 S =  2.1201185
Site (   8,   9) E0 = -7.5991877258528016 TruncErr = 3.77e-08 D =   300 Iter =   6 LanczT =     0.02 TotT =     0.02 S =  2.3855850
Site (   7,   8) E0 = -7.5991877787750335 TruncErr = 5.04e-08 D =   300 Iter =   6 LanczT =     0.02 TotT =     0.03 S =  2.2765200
Site (   6,   7) E0 = -7.5991877298614918 TruncErr = 3.67e-08 D =   300 Iter =   6 LanczT =     0.02 TotT =     0.02 S =  2.3855737
[timing] sweep  0.31892

Moved the center of MPS to 0.
Node4 will stop.
Node5 will stop.
Node6 will stop.
Node7 will stop.
Node1 will stop.
Node2 will stop.
Node3 will stop.
E0/site: -0.47495
CPU Time : 3.87674s
E0/site: 0
E0/site: 0
CPU Time : 4.36226s
E0/site: 0
CPU Time : 4.36268s
E0/site: 0
CPU Time : 4.37261s
E0/site: 0
CPU Time : 4.36255CPU Time : 4.37238s
E0/site: 0
CPU Time : 4.37252s
E0/site: 0
CPU Time : 4.37162s
s
```

The ground state wave function will be saved in the `mps` directory, with ground state energy reported at the end.

## 3. Measure Observables

After obtaining the MPS wave function, run:

```bash
path_to_build/measure1 quick_params.json
```

This will generate several JSON files with measurement results, including:

- `sz.json`: Local spin-z expectation values $\langle \hat{S}^z_i \rangle$
- `nf.json`: Local density expectation values $\langle \hat{n}_i \rangle$
- `sz0sz.json`: Spin-z correlation functions $\langle \hat{S}^z_0 \hat{S}^z_j \rangle$
- `sp0sm.json`: Spin raising/lowering correlations $\langle \hat{S}^+_0 \hat{S}^-_j \rangle$
- `sm0sp.json`: Spin lowering/raising correlations $\langle \hat{S}^-_0 \hat{S}^+_j \rangle$
- `nf0nf.json`: Density-density correlations $\langle \hat{n}_0 \hat{n}_j \rangle$
- `single_particle_corr[1-4].json`: Single-particle correlation functions ($c^\dagger_{i,\sigma}$ and $c_{j,\sigma}$ and $c_{i,\sigma}$ and $c^\dagger_{j,\sigma}$ with $\sigma$ the spin index)

### Example: Viewing Local Spin-z

```bash
cat sz.json
```

Output:
```
[
  [[0], 2.59576799587e-05],
  [[1], -2.19478820597e-05],
  [[2], -2.26287511573e-05],
  [[3], 2.14116129123e-05],
  [[4], 2.28057212783e-05],
  [[5], -1.87932101684e-05],
  [[6], -1.30490656793e-05],
  [[7], 1.1729093732e-05],
  [[8], -9.45153504613e-06],
  [[9], 1.06538555761e-05],
  [[10], 1.96032212238e-05],
  [[11], -2.23617681307e-05],
  [[12], -2.96192012065e-05],
  [[13], 2.82371681596e-05],
  [[14], 3.35207939214e-05],
  [[15], -3.60677333142e-05]
]
```
Here, each entry gives the value of $\langle \hat{S}^z_i \rangle$ for site $i$. 
The site indices $i$ correspond to the DMRG linear ordering (column-major mapping)

In all correlation files, site 0 is used as the reference site, and $j$ runs over all other sites.
To change the reference site, addtional argument can be added, e.g.
```bash
path_to_build/measure1 quick_params.json 2
```
This will change the reference site to 2

## 4. Analyze Results

Use Python, MATLAB, gnuplot, or your preferred tool to analyze and plot the results.
