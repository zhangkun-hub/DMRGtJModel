#!/bin/bash
#SBATCH --job-name=Zk_tJ_complex_measure
#SBATCH --output=slurm-%j-%x.log
#SBATCH --mail-user=12432064@mail.sustech.edu.cn
#SBATCH --mail-type=END,FAIL
#SBATCH --partition=256G56c
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=7

# =================================================================
# 复值生产复跑测量（measure1，任务 C 阶段 1）：自旋/密度/单粒子关联
# 用法：在复值 run 目录内提交 sbatch <repo>/scripts/rerun_complex_measure_corr.sh
# 与归档 measure_corr.sh 的差别：exedir 指向步骤 0 的 build_complex
# （icpc 复值二进制，测量输出 [re, im] 格式）。
# =================================================================
set -e # 让脚本在出错时立即退出

rundir=$(pwd)
rootdir="/share/home/zhangkun/DMRGtJModel"
exedir="${rootdir}/smoke_n64/run_20261104_022655/build_complex"
datadir="${rundir}/data_output"
plotsdir="${rundir}/plots"
paramsfile="params.json"

# ================================================================
# ★ 不同 bond 下多 site 扫描设置（只需要改这里！！！）
# ================================================================
bond=1500
sites=(49 50)


# =================================================================
# 运行环境准备（intel2022 口径，归档同源）
# =================================================================
echo "========================================================"
echo "Job Started at $(date) for Hybrid VMPS+DMRG (complex rerun)"
echo "Job ID: ${SLURM_JOB_ID}, Job Name: ${SLURM_JOB_NAME}"
echo "Running on nodes: ${SLURM_NODELIST}"
echo "MPI tasks: ${SLURM_NTASKS}, CPUs per task: ${SLURM_CPUS_PER_TASK}"
echo "========================================================"
ulimit -s unlimited
ulimit -l unlimited

module purge 2>/dev/null || true
export PATH=/share/software/gcc/12.1.0/bin:$PATH
export LD_LIBRARY_PATH=/share/software/gcc/12.1.0/lib64:$LD_LIBRARY_PATH
set +u
source /share/intel/2022.2.0/compiler/latest/env/vars.sh
source /share/intel/2022.2.0/mkl/latest/env/vars.sh
source /share/intel/2022.2.0/mpi/latest/env/vars.sh
source /share/software/anaconda3/bin/activate
set -u
module use /share/software/modulefiles 2>/dev/null || true
module load jq 2>/dev/null || true

[ -x "${exedir}/measure1" ] || { echo "FAIL: 复值 measure1 缺失：${exedir}/measure1"; exit 1; }

# =================================================================
# 多个site进行各种物理量的测量
# =================================================================
echo "--> Running measurement using existing complex MPS..."

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
