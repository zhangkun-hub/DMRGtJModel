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
# 复值生产复跑测量（measureSC_fix_ref_ext，任务 C 阶段 1）：超导四点关联
# 用法：在复值 run 目录内提交 sbatch <repo>/scripts/rerun_complex_measure_super.sh
# 与归档 measure_super_corr.sh 的差别：exedir 指向步骤 0 的 build_complex
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
# ★ 不同 bond 下多 refbonds 扫描设置（只需要改这里！！！）
# ================================================================
bond=1500
refbonds=("49,50" "50,51" "50,54")


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

[ -x "${exedir}/measureSC_fix_ref_ext" ] || \
  { echo "FAIL: 复值 measureSC_fix_ref_ext 缺失：${exedir}/measureSC_fix_ref_ext"; exit 1; }

# =================================================================
# 多个refbonds进行超导关联函数的计算
# =================================================================
echo "--> Running measurement using existing complex MPS..."

# 创建存储整个超导数据的文件夹
storage_dir="${rundir}/${bond}_measure_super_corr"
mkdir -p "${storage_dir}"

# 超导关联函数的计算
for ref in "${refbonds[@]}"; do
    echo "--------------------------------------------------------"
    echo ">>> Measuring superconducting correlation for refbonds = ${ref}"
    echo "--------------------------------------------------------"

    # 用于文件名标识
    tag=${ref/,/_}   #比如：0,4 -> 0_4

    # 设置计算结果所要输出的文件夹
    output_dir="${storage_dir}/super_corr_${tag}"
    mkdir -p "${output_dir}"

    # 运行测量
    mpirun -n ${SLURM_NTASKS} \
        ${exedir}/measureSC_fix_ref_ext "${paramsfile}" --refbond=${ref}

    echo ">>> Measurement finished for refbonds = ${ref}"

    # =============================================================
    # 移动输出文件
    # =============================================================
    echo "Moving output files to super_corr_${tag} directory..."

    # 存储需要移动的数据文件
    data_files_super=(
        "scsAa_fix_refbond_${tag}.json" "scsBa_fix_refbond_${tag}.json" "scsCa_fix_refbond_${tag}.json"
        "scsAb_fix_refbond_${tag}.json" "scsBb_fix_refbond_${tag}.json" "scsCb_fix_refbond_${tag}.json"
        "scsAc_fix_refbond_${tag}.json" "scsBc_fix_refbond_${tag}.json" "scsCc_fix_refbond_${tag}.json"
        "scsAd_fix_refbond_${tag}.json" "scsBd_fix_refbond_${tag}.json" "scsCd_fix_refbond_${tag}.json"
    )

    # 移动超导关联数据文件
    for file in "${data_files_super[@]}"; do
        if [ -f "${rundir}/${file}" ]; then
            mv "${rundir}/${file}" "${output_dir}/"
            echo "Moved ${file} to super_corr_${tag} directory"
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
