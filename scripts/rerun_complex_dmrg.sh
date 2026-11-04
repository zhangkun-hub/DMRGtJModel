#!/bin/bash
#SBATCH --job-name=Zk_tJ_complex
#SBATCH --output=slurm-%j-%x.log
#SBATCH --mail-user=12432064@mail.sustech.edu.cn
#SBATCH --mail-type=END,FAIL
#SBATCH --partition=256G56c
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=7

# =================================================================
# 复值生产复跑 DMRG（任务 C 阶段 1，试点 d2_25）
# 用法：在 run 目录内提交 sbatch <repo>/scripts/rerun_complex_dmrg.sh
#   （run 目录内须已有 params.json —— 以
#   scripts/rerun_complex_params_d2_25.json 为模板；脚本用 jq 生成两阶段
#   参数副本，不改动模板语义）
# 二进制：步骤 0 冒烟的 build_complex（icpc 复值 + OpenMP，已验证）
# 流程（两步绝热链，d2_25 归档口径）：
#   ① t2=J2=0.01（Néel 初态从头）→ 删 .temp → ② t2=J2=0 纯 t-J（MPS 续跑）
# 注意：含 rm -rf .temp —— 仅删可重建的临时环境张量（skill 第 3/4 节规则：
#   换物理参数前必须删；删除前校验目录存在且 basename 为 .temp）。
# =================================================================
set -euo pipefail

rundir=$(pwd)
rootdir="/share/home/zhangkun/DMRGtJModel"
COMPLEX_BUILD="${rootdir}/smoke_n64/run_20261104_022655/build_complex"
paramsfile="params.json"
###### 每次要修改的参数
D_VALUES="1500"          # d2_25 低档口径：单值 1500
SWEEPS=10                # d2_25 口径

die(){ echo "FAIL: $*" >&2; exit 1; }

# ---------- 运行环境（intel2022 口径，submit_dmrg_zk.sh 同源） ----------
ulimit -s unlimited
ulimit -l unlimited
module purge 2>/dev/null || true
export PATH=/share/software/gcc/12.1.0/bin:$PATH
export LD_LIBRARY_PATH=/share/software/gcc/12.1.0/lib64:$LD_LIBRARY_PATH
# Intel vars.sh 内部脚本会读取未定义变量，source 期间临时关闭 nounset。
set +u
source /share/intel/2022.2.0/compiler/latest/env/vars.sh
source /share/intel/2022.2.0/mkl/latest/env/vars.sh
source /share/intel/2022.2.0/mpi/latest/env/vars.sh
source /share/software/anaconda3/bin/activate
set -u
module use /share/software/modulefiles 2>/dev/null || true
module load jq 2>/dev/null || true
command -v mpirun >/dev/null 2>&1 || die "mpirun 不可用"
command -v jq >/dev/null 2>&1 || die "jq 不可用（anaconda 或 jq 模块）"
[ -x "${COMPLEX_BUILD}/dmrg" ] || die "复值 dmrg 缺失：${COMPLEX_BUILD}/dmrg"
[ -f "${paramsfile}" ] || die "params.json 缺失（请先按 rerun_complex_params_d2_25.json 准备）"

echo "========================================================"
echo "Complex rerun DMRG: $(date)"
echo "Job ID: ${SLURM_JOB_ID}, Nodelist: ${SLURM_NODELIST}"
echo "rundir: $(pwd), D: ${D_VALUES}, Sweeps: ${SWEEPS}"
echo "========================================================"

# ---------- 两阶段参数（jq 生成副本，基础模板不动） ----------
jq '.CaseParams.t2 = 0.01 | .CaseParams.J2 = 0.01 | .CaseParams.Sweeps = '"${SWEEPS}" \
  "${paramsfile}" > params_stage1.json
jq '.CaseParams.t2 = 0.0  | .CaseParams.J2 = 0.0  | .CaseParams.Sweeps = '"${SWEEPS}" \
  "${paramsfile}" > params_stage2.json

# ---------- ① 邻域（t2=J2=0.01），从头 Néel 初态 ----------
echo "--> stage 1: t2=J2=0.01, D=${D_VALUES}"
cp params_stage1.json "${paramsfile}"
mpirun -n ${SLURM_NTASKS} "${COMPLEX_BUILD}/dmrg" "${paramsfile}" --D="${D_VALUES}"
echo "     stage 1 finished."

# ---------- 换参数前删除 .temp（skill 规则：可重建，须删） ----------
TEMP_DIR="${rundir}/.temp"
if [ -d "${TEMP_DIR}" ] && [ "$(basename "${TEMP_DIR}")" = ".temp" ]; then
  echo "--> removing .temp before parameter change (rebuildable per skill rule)"
  rm -rf "${TEMP_DIR}"
fi

# ---------- ② 纯 t-J（t2=J2=0），MPS 续跑 ----------
echo "--> stage 2: pure t-J, D=${D_VALUES}"
cp params_stage2.json "${paramsfile}"
mpirun -n ${SLURM_NTASKS} "${COMPLEX_BUILD}/dmrg" "${paramsfile}" --D="${D_VALUES}"
echo "     stage 2 finished."

echo "========================================================"
echo "Complex rerun DMRG finished at $(date)"
echo "========================================================"
