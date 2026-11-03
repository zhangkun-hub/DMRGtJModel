#!/bin/bash
# =====================================================================
# N=64 冒烟 verify（C-only）：前置检查（run 完成标记/配置一致性/结果完整性）
#   → v1 复值重现实值（scs 单参考键，独立 manifest）→ v2 复值重现实值
#   （measure1 目录，generic 模式）。
# 用法：sbatch scripts/smoke_n64_verify.sh <RUN_DIR>   （RUN_DIR 由 run 脚本打印）
# 安全设计：run 未完成/配置不一致/必要目录缺失 → 明确 FAIL（不会误 PASS）；
#   每次验证使用唯一输出目录 verify_<时间>，不覆盖历史报告。
# =====================================================================
#SBATCH --job-name=smoke_n64_verify
#SBATCH --partition=256G56c
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=0:30:00
#SBATCH --output=smoke_n64_verify_%j.out

set -euo pipefail

REPO_ROOT=${REPO_ROOT:-/share/home/zhangkun/DMRGtJModel}
CONFIG_SRC=${REPO_ROOT}/scripts/smoke_n64_config.json
COMMON=${REPO_ROOT}/scripts/smoke_common.sh

# shellcheck source=/dev/null
source "${COMMON}"

# ---------- 前置检查（fail-closed；先于环境设置，未完成的 run 快速失败） ----------
RUN_DIR=${1:-}
if [ -z "${RUN_DIR}" ]; then
  echo "usage: $0 <RUN_DIR>（RUN_DIR 由 smoke_n64_run.sh 结束时打印）" >&2
  exit 2
fi
[ -d "${RUN_DIR}" ] || die "运行目录不存在：${RUN_DIR}"
[ -f "${RUN_DIR}/RUN_OK" ] || die "run 未完成（RUN_OK 标记缺失）：${RUN_DIR}"
[ -f "${RUN_DIR}/smoke_env_record.txt" ] || die "环境记录缺失：${RUN_DIR}/smoke_env_record.txt"
[ -f "${RUN_DIR}/smoke_config_used.json" ] || die "配置快照缺失：${RUN_DIR}/smoke_config_used.json"
cmp -s "${RUN_DIR}/smoke_config_used.json" "${CONFIG_SRC}" \
  || die "run 使用的配置与当前 ${CONFIG_SRC} 不一致，拒绝验证"

# ---------- 环境（与用户 ${HOME}/myenv.sh 同源） ----------
module purge 2>/dev/null || true
if [ -f "${HOME}/myenv.sh" ]; then
  source "${HOME}/myenv.sh"
else
  die "缺少 ${HOME}/myenv.sh（cmake/gcc/jq 模块 + oneAPI setvars 的集群环境初始化脚本）"
fi
command -v "${PYTHON_BIN}" >/dev/null 2>&1 \
  || die "python3 不可用（请在 ${HOME}/myenv.sh 中补充 python 模块或确认系统 python3）"

# ---------- 配置（与 run 同源：smoke_n64_config.json + params.json） ----------
PARAMS_FILE="${RUN_DIR}/run_real/params.json"
[ -f "${PARAMS_FILE}" ] || die "参数文件缺失：${PARAMS_FILE}"
GEOM=$(px "${PARAMS_FILE}" Geometry)
LY=$(px "${PARAMS_FILE}" Ly)
LX=$(px "${PARAMS_FILE}" Lx)
SINGLE_REF=$(cfg_get "${CONFIG_SRC}" single_refbond)

# ---------- 必要结果与日志完整性 ----------
for d in \
    run_real/measure1_D300 \
    run_real/measure_scs_0_1_D100 \
    run_real/measure_scs_0_1_D300 \
    run_complex/measure1_D300 \
    run_complex/measure_scs_0_1_D300; do
  [ -d "${RUN_DIR}/${d}" ] || die "结果目录缺失：${RUN_DIR}/${d}"
done
for l in run_real/dmrg_100.log run_real/dmrg_300.log \
         run_complex/dmrg_100.log run_complex/dmrg_300.log; do
  [ -f "${RUN_DIR}/${l}" ] || die "dmrg 日志缺失：${RUN_DIR}/${l}"
done

# ---------- 唯一验证输出目录（并发/重跑互不覆盖） ----------
OUT="${RUN_DIR}/verify_$(date +%H%M%S)"
mkdir "${OUT}" || die "验证输出目录创建失败（并发冲突？）：${OUT}"

VERIFY=${REPO_ROOT}/tools/verify_complex_reproduces_real.py

echo "==== v1: 复值重现实值（scs 单参考键，独立 manifest）===="
${PYTHON_BIN} "${VERIFY}" \
  --complex "${RUN_DIR}/run_complex/measure_scs_0_1_D300" \
  --real    "${RUN_DIR}/run_real/measure_scs_0_1_D300" \
  --complex-log "${RUN_DIR}/run_complex/dmrg_300.log" \
  --real-log    "${RUN_DIR}/run_real/dmrg_300.log" \
  --scs-mode --geometry "${GEOM}" --Ly "${LY}" --Lx "${LX}" \
  --refbonds "${SINGLE_REF}" \
  -o "${OUT}/verify_scs_0_1.json"

echo "==== v2: 复值重现实值（measure1 目录，generic 模式）===="
${PYTHON_BIN} "${VERIFY}" \
  --complex "${RUN_DIR}/run_complex/measure1_D300" \
  --real    "${RUN_DIR}/run_real/measure1_D300" \
  --complex-log "${RUN_DIR}/run_complex/dmrg_300.log" \
  --real-log    "${RUN_DIR}/run_real/dmrg_300.log" \
  -o "${OUT}/verify_measure1.json"

echo "SMOKE VERIFY ALL OK"
