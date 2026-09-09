#!/bin/bash
#SBATCH --job-name=Zk_tJ
#SBATCH --output=slurm-%j-%x.log
#SBATCH --mail-user=12432064@mail.sustech.edu.cn
#SBATCH --mail-type=END,FAIL
#SBATCH --partition=256G56c
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=7

# =================================================================
# 定义工作路径和参数 (Define Paths and Parameters)
# =================================================================
set -e # 让脚本在出错时立即退出

rundir=$(pwd)
rootdir="/share/home/zhangkun/DMRGtJModel"
exedir="${rootdir}/build"
paramsfile="params.json"
###### 每次要修改的参数
D_VALUES_DMRG="100,500,800,1000,1200,1500" 
# 含义：表示bond dimension,从小到大排列(中间不要加空格),例如:D_VALUES_DMRG="10,50,100"


# =================================================================
# 运行环境准备 (Environment Setup)
# =================================================================
echo "========================================================"
echo "Job Started at $(date) for Hybrid VMPS+DMRG"
echo "Job ID: ${SLURM_JOB_ID}, Job Name: ${SLURM_JOB_NAME}"
echo "Running on nodes: ${SLURM_NODELIST}"
echo "shareing Directory: $(pwd)"
echo "MPI tasks: ${SLURM_NTASKS}, CPUs per task: ${SLURM_CPUS_PER_TASK}"
echo "========================================================"
ulimit -s unlimited
ulimit -l unlimited

# **核心关键**：加载我们确定好的、正确的环境模块！
# (Critical: Load the correct environment modules!)
# 手动加载GCC 12.1.0
export PATH=/share/software/gcc/12.1.0/bin:$PATH
export LD_LIBRARY_PATH=/share/software/gcc/12.1.0/lib64:$LD_LIBRARY_PATH
# 验证环境
echo "=== GCC版本 ==="
gcc --version
echo "=== GLIBCXX支持的版本 ==="
strings /share/software/gcc/12.1.0/lib64/libstdc++.so.6 | grep GLIBCXX_3.4.29
# 设置Intel编译器环境
echo "Setting up Intel compiler environment..."
source /share/intel/2022.2.0/compiler/latest/env/vars.sh
# 设置MKL环境
echo "Setting up MKL environment..."
source /share/intel/2022.2.0/mkl/latest/env/vars.sh
# 设置MPI环境
echo "Setting up MPI environment..."
source /share/intel/2022.2.0/mpi/latest/env/vars.sh
# **核心修改**: 激活 Anaconda3 环境来使用 Python 3
# (Critical Change: Activate Anaconda3 to use Python 3)
echo "Activating Anaconda3 environment..."
source /share/software/anaconda3/bin/activate

module use /share/software/modulefiles
module load jq


echo ""
echo "========================================================"
echo "Bond Dimensions for DMRG: ${D_VALUES_DMRG}"


echo "--> Step 2.2: Running main DMRG calculation using existing MPS..."
mpirun -n ${SLURM_NTASKS} ${exedir}/dmrg "${paramsfile}" --D=${D_VALUES_DMRG}
echo "     DMRG calculation finished."


echo "========================================================"
echo "Job Finished at $(date)"
echo "========================================================"