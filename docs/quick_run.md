# Quick Run (MPI on CPU Architecture)

This guide demonstrates how to perform a simple calculation with the Square Lattice t-J Model DMRG code after building the project.  
We show how to perform a DMRG ground state calculation on a $2\times 4$ square lattice with 2 holes doped away from half-filling.  
**Before starting, ensure you have compiled the code with particle number U(1) symmetry conservation.**

## 0. Create a New Directory for the Calculation

We recommend creating a separate directory for each calculation to keep parameter files, wave function data, and results organized.

## 1. Set Up the Parameter File

Create a JSON file (e.g., `quick_params.json`) with the following content:

```json
{
  "CaseParams": {
    "Geometry": "OBC",
    "Lx": 2,
    "Ly": 4,
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

This example sets up a $2\times 4$ lattice for the t-J model with $t=1$, $J=0.3$, and no next-nearest-neighbor terms ($t_2=0$, $J_2=0$).

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
Setting MaxLanczIter as : [3, 25, 47, 70]
lattice construct
Setting MaxLanczIter as : [3, 25, 47, 70]
lattice construct
Setting MaxLanczIter as : [3, 25, 47, 70]
lattice construct
Setting MaxLanczIter as : [3, 25, 47, 70]
lattice construct
Setting MaxLanczIter as : [3, 25, 47, 70]
lattice construct
Setting MaxLanczIter as : [3, 25, 47, 70]
Setting MaxLanczIter as : [3, 25, 47, 70]
lattice construct
lattice construct
Setting MaxLanczIter as : [3, 25, 47, 70]
lattice construct
[	9	18	26	34	26	18	9	1]
[	9	18	26	34	26	18	9	1]
[	9	18	26	34	26	18	9	1]
[	9	18	26	34	26	18	9	1]
[	9	18	26	34	26	18	9	1]
[	9	18	26	34	26	18	9	1]
[	9	18	26	34	26	18	9	1]
[	9	18	26	34	26	18	9	1]
Initial mps as direct product state.
D_max = 10
Initial mps as direct product state.
Initial mps as direct product state.
Initial mps as direct product state.
Initial mps as direct product state.
Initial mps as direct product state.

=====> MPI Sweep Parameters For DMRG <=====
MPS/MPO size: 	 8
Sweep times: 	 4
Bond dimension: 	 10/10
Truncation error: 	 1e-08
Lanczos max iterations 	3
MPS path: 	mps
Temp path: 	.temp
=====> Technical Parameters <=====
The number of processors(including master): 	8
The number of threads per processor: 	4
=====> Checking and Updating Boundary Tensors =====>
Initial mps as direct product state.
Initial mps as direct product state.
left boundary: 	2
right boundary: 	5
No temp directory .temp
=====> Generating the block operators =====>
right block operator number : 9
right block operator number : 18
right block operator number : 26
right block operator number : 34

sweep 1
```

and output end with
```
sweep 4
Site (   2,   3) E0 = -4.9258005218948036 TruncErr = 0.00e+00 D =    26 Iter =   2 LanczT =     0.00 TotT =     0.00 S =  2.2616350
Site (   3,   4) E0 = -4.9258005218984087 TruncErr = 0.00e+00 D =    60 Iter =   2 LanczT =     0.00 TotT =     0.01 S =  2.4634575
Site (   4,   5) E0 = -4.9258005219017429 TruncErr = 0.00e+00 D =    26 Iter =   2 LanczT =     0.00 TotT =     0.00 S =  2.2616351
Site (   3,   4) E0 = -4.9258005219049945 TruncErr = 0.00e+00 D =    60 Iter =   2 LanczT =     0.00 TotT =     0.01 S =  2.4634576
[timing] sweep  0.02455

Moved the center of MPS to 0.
Node1 will stop.
Node2 will stop.
Node3 will stop.
Node4 will stop.
Node5 will stop.
Node6 will stop.
Node7 will stop.
E0/site: -0.61573
CPU Time : 3.33102s
E0/site: 0
CPU Time : 4.06455s
E0/site: 0
CPU Time : 4.06112s
E0/site: 0
CPU Time : 4.01902s
E0/site: 0
CPU Time : 4.0338s
E0/site: 0
CPU Time : 4.05581s
E0/site: 0
CPU Time : 4.05386s
E0/site: 0
CPU Time : 4.05449s
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
  [[0], 1.29235573822e-05],
  [[1], -6.3974631298e-06],
  [[2], -6.39746312949e-06],
  [[3], 1.29235573831e-05],
  [[4], -1.29254624046e-05],
  [[5], 6.39936815239e-06],
  [[6], 6.39936815247e-06],
  [[7], -1.29254624063e-05]
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
