#!/bin/bash
#SBATCH --job-name=Zk_tJ_measure
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
datadir="${rundir}/data_output"
plotsdir="${rundir}/plots"
paramsfile="params.json"

# ================================================================
# ★ 不同 bond 下多 site 扫描设置（只需要改这里！！！）
# ================================================================
bond=5000
sites=(49 50) 



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



# =================================================================
# 多个site进行各种物理量的测量
# =================================================================
echo "--> Step 2.2: Running measurement using existing MPS..."

# 创建存储整个测量数据的文件夹
storage_dir="${rundir}/${bond}_measure_phys_qty"
mkdir -p "${storage_dir}"


# 自旋、密度、单粒子关联函数、自旋关联函数和密度关联函数的计算
for site in "${sites[@]}"; do
    echo "--------------------------------------------------------"
    echo ">>> Measuring physical quantity for site = ${site}"
    echo "--------------------------------------------------------"

    # 设置计算结果所要输出的文件夹
    output_dir="${storage_dir}/measure_${site}"
    mkdir -p "${output_dir}"

    # 运行测量
    ${exedir}/measure1 "${paramsfile}" ${site}
    echo ">>> Measurement finished for site = ${site}"



    # ==============================================================================
    # 移动输出文件
    # ==============================================================================
    echo "Moving output files to measure_${site} directory..."

    # 存储需要移动的数据文件
    data_files_meas=(
        "nf.json" "nf${site}nf.json" "single_particle_corr1.json" 
        "single_particle_corr2.json" "single_particle_corr3.json" 
        "single_particle_corr4.json" "sm${site}sp.json" "sp${site}sm.json" 
        "sz.json" "sz${site}sz.json"
    )

    # 移动测量数据文件
    for file in "${data_files_meas[@]}"; do
        if [ -f "${rundir}/${file}" ]; then
            mv "${rundir}/${file}" "${output_dir}/"
            echo "Moved ${file} to measure_${site} directory"
        else
            echo "Warning: ${file} not found, skipping"
        fi
    done
done

# =================================================================
# 整理 SLURM 日志
# =================================================================
slurm_log="slurm-${SLURM_JOB_ID}-${SLURM_JOB_NAME}.log"
if [ -f "${rundir}/${slurm_log}" ]; then
    mv "${rundir}/${slurm_log}" "${storage_dir}"
fi




echo "========================================================"
echo "Job Finished at $(date)"
echo "========================================================"