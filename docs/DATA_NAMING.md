# Data Naming Convention

This document describes the standard naming convention used for simulation data files in the `data/` directory.

## File Name Structure

All data files follow a standard pattern that encodes the simulation parameters.

Format:
`[Prefix][Ly]x[Lx]t2=[t2]J=[J][FillingType][FillingValue]D[BondDim][Geometry].json`

### Components

| Component | Description | Example | Notes |
|-----------|-------------|---------|-------|
| **Prefix** | The physical quantity measured. | `nf` (Charge), `sz` (Spin), `scsyy` (Pairing) | |
| **Ly** | Lattice height (Y direction). | `6` | |
| **Lx** | Lattice width (X direction). | `16` | |
| **t2** | Next-nearest neighbor hopping amplitude. | `0`, `-0.2` | **Mandatory**. If 0, must be written as `t2=0`. |
| **J** | Nearest neighbor exchange interaction. | `0.4`, `3` | **Normalization Rule**: Integers must NOT have decimal points (e.g., use `J=3`, NOT `J=3.0`). Floats retain decimals. |
| **FillingType** | Controls particle density. | `Hole` (number of holes) or `Mu` (chemical potential) | |
| **FillingValue** | The value associated with the FillingType. | `6` (for Hole) or `1.5` (for Mu) | |
| **D** | Bond dimension (MPS accuracy parameter). | `5000` | |
| **Geometry** | Boundary conditions. | `OBC`, `Cylinder` | |

### Examples

1.  **Charge Density (nf) with Hole Doping:**
    `nf0nf6x6t2=0J=0.4Hole6D50000OBC.json`
    *   Prefix: `nf0nf` (Charge correlation/density)
    *   Size: 6x6
    *   t2: 0
    *   J: 0.4
    *   Holes: 6
    *   Bond Dim: 50000
    *   Geometry: OBC

2.  **Spin Density (sz) with Chemical Potential:**
    `sz0sz8x16t2=0J=0.4Mu1.5D15000OBC.json`
    *   Prefix: `sz0sz`
    *   Size: 8x16
    *   Filling: Mu = 1.5

### Number Formatting Rules
To ensure consistency:
- **Integers**: Must be represented without a decimal point. 
  - Correct: `J=3`
  - Incorrect: `J=3.0`
- **Floats**: Use standard decimal notation.
  - Correct: `J=0.4`

## Helper Functions

To avoid hardcoding these strings in analysis scripts, use the provided MATLAB helper:

```matlab
% In MATLAB:
postfix = GetSimulationParamsPostfix(Ly, Lx, t2, J, 'Hole', 2, 5000, 'OBC');
filename = ['nf', postfix];
```
