# Simulation Parameters (`params.json`)

The simulation parameter files for the Square Lattice t-J Model are exemplified in the `src/params.json` file.  
Below is a description of each parameter and its role in the simulation.

---

## Parameter Reference

| Parameter                     | Type    | Description                                                                                                                         |
|-------------------------------|---------|-------------------------------------------------------------------------------------------------------------------------------------|
| `Geometry`                    | string  | Lattice geometry. Options: `"SquareCylinder"`, `"SquareOBC"`, `"SquareRotated"`, `"SquareTorus"`, `"HoneycombXC"`, `"HoneycombYC"`. Legacy values (`"Cylinder"`, `"OBC"`, `"Rotated"`, `"Torus"`) remain accepted but are normalized to the corresponding `Square*` entries with a warning; prefer the new names in fresh configurations. |
| `Lx`                          | integer | Number of sites along the x-direction.                                                                                              |
| `Ly`                          | integer | Number of sites along the y-direction.                                                                                              |
| `t`                           | float   | Nearest-neighbor hopping amplitude.                                                                                                 |
| `t2`                          | float   | Next-nearest-neighbor hopping amplitude.                                                                                            |
| `J`                           | float   | Nearest-neighbor spin exchange coupling.                                                                                            |
| `J2`                          | float   | Next-nearest-neighbor spin exchange coupling.                                                                                       |
| `phi`                         | float   | Magnetic flux inserted into cylinder (if applicable). Ignored (with warning) for honeycomb geometries.                              |
| `mu`                          | float   | Chemical potential, in Hamiltonian it contribute $-\mu N$ (if no particle number conservation code).                                |
| `NumHole`                     | integer | Number of holes, counting in the way away from half-filling                                                                         |
| `Sweeps`                      | integer | Number of DMRG sweeps.                                                                                                              |
| `Dmin`                        | integer | Minimum bond dimension for DMRG.                                                                                                    |
| `Dmax`                        | integer | Maximum bond dimension for DMRG.                                                                                                    |
| `CutOff`                      | float   | Truncation cutoff for DMRG.                                                                                                         |
| `LanczErr`                    | float   | Tolerance for Lanczos eigensolver.                                                                                                  |
| `MaxLanczIter`                | integer | Maximum number of Lanczos iterations.                                                                                               |
| `Threads`                     | integer | Number of computational threads for each process to use.                                                                            |
| `noise`                       | array   | List of noise amplitudes for DMRG sweeps (e.g., `[0.1, 0.01]`) (not applicable for tranditional DMRG).                              |
| `Perturbation`                | float   | Not applicable.                                                                                                                     |
| `BondSingletPairPerturbation` | float   | Strength of a bond singlet pair perturbation on boundary of the lattice. (only applicable if no particle number conservation code). |
| `wavelength`                  | integer | Not applicable.                                                                                                                     |

---

## Example

```json
{
  "CaseParams": {
    "Geometry": "SquareOBC",
    "Lx": 3,
    "Ly": 4,
    "t": 1.0,
    "t2": -0.2,
    "J": 3.0,
    "J2": 1.0,
    "phi": 0.1,
    "mu": 0,
    "NumHole": 2,
    "Sweeps": 5,
    "Dmin": 10,
    "Dmax": 10,
    "CutOff": 1e-9,
    "LanczErr": 1e-9,
    "MaxLanczIter": 70,
    "Threads": 4,
    "noise": [
      0.1,
      0.01
    ],
    "Perturbation": 0.0,
    "BondSingletPairPerturbation": 0.0,
    "wavelength": 4
  }
}
```

---

## Usage

- Edit `src/params.json` or copy another one and edit to set up your simulation.
- Pass the file as an argument to the binaries, e.g.:
  ```bash
  ./dmrg src/params.json
  ```

---

## Command-Line Options

In addition to `params.json`, the `dmrg` binary accepts the following command-line flags:

| Option | Default | Description |
|---|---|---|
| `--D=d1,d2,...` | (use JSON `Dmin`/`Dmax`) | Bond dimension schedule. Comma-separated, ascending, no spaces. The program runs one sweep pass per value. `MaxLanczIter` is ramped linearly from 3 to the JSON value across the list. |
| `--complex-init-eps=<value>` | 0 (disabled) | Complex perturbation amplitude $\varepsilon$. Must be non-negative. When `> 0`, constructs a superposition initial MPS: $|\Psi\rangle \propto |\text{Neel}\rangle + \varepsilon_c |b\rangle$ where $\varepsilon_c = \varepsilon (1+i)/\sqrt{2}$. **Requires**: complex build (`RealCode=OFF`), honeycomb geometry, no existing `./mps` directory. |
| `--complex-init-seed=<uint>` | 12345 | Random seed for selecting the NN swap link used in the perturbation. Must be a non-negative integer in the range `[0, 4294967295]` (`UINT32_MAX`). Fixed seed ensures reproducibility across runs. |

**Notes:**

- `--complex-init-eps` is independent of the DMRG noise mechanism (`noised_vmps`). The former modifies the initial wavefunction; the latter adds noise during the optimization sweeps.
- If `--complex-init-eps > 0` and `./mps` already exists, `dmrg` exits with an error to prevent overwriting converged results. Remove or rename `./mps` first.
- The perturbation swaps a pair of opposite spins on a real honeycomb nearest-neighbor bond, preserving total particle number and total `Sz`.