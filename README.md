# Square Lattice t-J Model : DMRG

## The t-J Model:

$$
\begin{aligned}
H=&-\sum_{\{i j\}, \sigma} t_{i j}\left(\hat{c}_{i, \sigma}^{\dagger} \hat{c}_{j, \sigma}+\text{H.c.}\right)\\
&+\sum_{\{i j\}} J_{i j}\left(\hat{\mathbf{S}}_{i} \cdot \hat{\mathbf{S}}_{j}-\frac{1}{4} \hat{n}_{i} \hat{n}_{j}\right),
\end{aligned},
$$
with $\{i j\}$ upto NNN pairs.

___


## Binaries and Their Usage
After compiling the programs according the **Build Instructions** below,
you can obtain the following binaries.

| Binary Name   | Functionality (a short explanation)                                       | Usage                                                  |
|---------------|---------------------------------------------------------------------------|--------------------------------------------------------|
| **dmrg**      | Traditional DMRG program for the model.                                   | mpirun -np ${num_of_processor} ./dmrg params.json      |
| **mpogen**    | Generate MPO for the model and dump into disk.                            | ./mpogen params.json                                   |
| **vmps**      | MPO-based DMRG program for the model.                                     | mpirun -np ${num_of_processor} ./vmps params.json      |
| **vmps1FM**   | MPO-based DMRG program for the model, with fully-polarized initial state. | mpirun -np ${num_of_processor} ./vmps1FM params.json   |
| **measure1**  | Measure one-site and two-site correlation functions on MPS                | ./measure1 params.json                                 |
| **measureSC** | Measure superconductivity correlation functions on MPS                    | mpirun -np ${num_of_processor} ./measureSC params.json |

For detailed usages, see docs/tutorial.md

---

## To-Do List

- [ ] 

___

## Dependence

- **C++ Compiler:** C++17 or above
- **Build System:** CMake (version 3.12 or higher)
- **Math Libraries:** Intel MKL or OpenBLAS
- **Parallelization:** MPI
- **Tensor:** [QuantumLiquids/TensorToolkit](https://github.com/QuantumLiquids/TensorToolkit)
- **DMRG:** [QuantumLiquids/UltraDMRG](https://github.com/QuantumLiquids/UltraDMRG)
- **GPU Acceleration (optional):** CUDA compiler, cuBLAS, cuSolver, cuTensor2

___

## Build Instructions

To compile the project using CMake, follow these steps:

1. Clone the repository:
   ```bash
   git clone https://github.com/QuantumLiquids/SquareLatticetJmodel.git
   cd SquareLatticetJmodel
   ```

2. Build using CMake:
   ```bash
   mkdir build && cd build
   cmake .. -DCMAKE_CXX_COMPILER=your_cxx_compiler \
            -DCMAKE_PREFIX_PATH=hint_to_find_tensortoolkit_and_ultradmrg \
            -DQLMPS_USE_GPU=ON/OFF \
            -DCUTENSOR_ROOT=path_to_cutensor_if_use_cuda
   make -j4 
   ```
   **Common CMake options:**

    - `-DHIGHER_SYMMETRY=ON` : Use particle number conserved code (default) 
      `-DHIGHER_SYMMETRY=OFF` Use particle number fluctuated code
    - `-DRealCode=ON` : Use real number tensors 
      `-DRealCode=OFF` : Use complex number tensors (default)

   **BLAS/LAPACK settings:**
    - On **x86 platforms**, CMakeList.txt will automatically use oneAPI Parallel MKL as the BLAS implementation.
    - On **ARM platforms**, it will default to OpenBLAS.
    - If you want to link against a different implementation, such as an older version of MKL, sequential MKL, or AOCL, you will need to manually set the following CMake options:
      - *BLA_VENDOR*: [CMake Manual](https://cmake.org/cmake/help/latest/module/FindBLAS.html) 
      - *BLAS_INCLUDE_DIR*: Set the header paths for the chosen BLAS/LAPACK implementation.

---

## Author

For inquiries, questions, or collaboration opportunities, please contact Hao-Xin via email:
[wanghaoxin1996@gmail.com](mailto:wanghaoxin1996@gmail.com).

## Acknowledgments
