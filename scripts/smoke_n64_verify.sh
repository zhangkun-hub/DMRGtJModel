#!/bin/bash
# =====================================================================
# N=64 冒烟 verify：前置检查（run 完成标记/配置一致性/结果完整性）→
#   v1-v3 复值重现实值（scs 独立 manifest + generic）→ v4 共轭 R3 →
#   v5 A1 schema → v6 TruncErr → v7 D 收敛敏感性。
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

module purge
module load Anaconda3
command -v "${PYTHON_BIN}" >/dev/null 2>&1 || die "python3 不可用"

# ---------- 前置检查（fail-closed） ----------
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

# ---------- 配置（与 run 同源：smoke_n64_config.json + params.json） ----------
PARAMS_FILE="${RUN_DIR}/run_real/params.json"
[ -f "${PARAMS_FILE}" ] || die "参数文件缺失：${PARAMS_FILE}"
GEOM=$(px "${PARAMS_FILE}" Geometry)
LY=$(px "${PARAMS_FILE}" Ly)
LX=$(px "${PARAMS_FILE}" Lx)
SINGLE_REF=$(cfg_get "${CONFIG_SRC}" single_refbond)
MULTI_REF=$(cfg_list "${CONFIG_SRC}" refbonds | paste -sd ';' -)
REVERSE=$(cfg_get "${CONFIG_SRC}" reverse_channels)
REF_TAG=$(echo "${SINGLE_REF}" | tr ',' '_')
MULTI_TAG=$([ "${REVERSE}" = "true" ] && echo two_rev || echo two)
REV_FLAG=()
[ "${REVERSE}" = "true" ] && REV_FLAG=(--reverse-channels)

# ---------- 必要结果与日志完整性 ----------
for d in \
    run_real/measure1_D300 \
    run_real/measure_scs_0_1_D100 \
    run_real/measure_scs_0_1_D300 \
    run_real/measure_scs_${MULTI_TAG} \
    run_complex/measure1_D300 \
    run_complex/measure_scs_0_1_D300 \
    run_complex/measure_scs_${MULTI_TAG}; do
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
CONJ=${REPO_ROOT}/tools/check_pairing_conjugacy.py
TRUNC=${REPO_ROOT}/tools/extract_trunc_err.py
BONDDIM=${REPO_ROOT}/tools/compare_bond_dim.py

echo "==== v1: 复值重现实值（scs 单参考键，独立 manifest）===="
${PYTHON_BIN} "${VERIFY}" \
  --complex "${RUN_DIR}/run_complex/measure_scs_0_1_D300" \
  --real    "${RUN_DIR}/run_real/measure_scs_0_1_D300" \
  --complex-log "${RUN_DIR}/run_complex/dmrg_300.log" \
  --real-log    "${RUN_DIR}/run_real/dmrg_300.log" \
  --scs-mode --geometry "${GEOM}" --Ly "${LY}" --Lx "${LX}" \
  --refbonds "${SINGLE_REF}" \
  -o "${OUT}/verify_scs_0_1.json"

echo "==== v2: 复值重现实值（多参考键 ± reverse，独立 manifest）===="
${PYTHON_BIN} "${VERIFY}" \
  --complex "${RUN_DIR}/run_complex/measure_scs_${MULTI_TAG}" \
  --real    "${RUN_DIR}/run_real/measure_scs_${MULTI_TAG}" \
  --complex-log "${RUN_DIR}/run_complex/dmrg_300.log" \
  --real-log    "${RUN_DIR}/run_real/dmrg_300.log" \
  --scs-mode --geometry "${GEOM}" --Ly "${LY}" --Lx "${LX}" \
  --refbonds "${MULTI_REF}" "${REV_FLAG[@]}" \
  -o "${OUT}/verify_scs_multi.json"

echo "==== v3: 复值重现实值（measure1 目录，generic 模式）===="
${PYTHON_BIN} "${VERIFY}" \
  --complex "${RUN_DIR}/run_complex/measure1_D300" \
  --real    "${RUN_DIR}/run_real/measure1_D300" \
  --complex-log "${RUN_DIR}/run_complex/dmrg_300.log" \
  --real-log    "${RUN_DIR}/run_real/dmrg_300.log" \
  -o "${OUT}/verify_measure1.json"

echo "==== v4: 共轭关系 R3（complex 与 real 的多参考键目录）===="
${PYTHON_BIN} "${CONJ}" --dir "${RUN_DIR}/run_complex/measure_scs_${MULTI_TAG}" \
  --geometry "${GEOM}" --Ly "${LY}" --Lx "${LX}" --rtol 1e-6 \
  -o "${OUT}/conjugacy_complex.json"
${PYTHON_BIN} "${CONJ}" --dir "${RUN_DIR}/run_real/measure_scs_${MULTI_TAG}" \
  --geometry "${GEOM}" --Ly "${LY}" --Lx "${LX}" --rtol 1e-6 \
  -o "${OUT}/conjugacy_real.json"

echo "==== v5: A1 spread/geometry schema 检查 ===="
${PYTHON_BIN} - "${RUN_DIR}/run_real/measure_scs_0_1_D300" \
          "${RUN_DIR}/run_complex/measure_scs_0_1_D300" \
          "${REF_TAG}" <<'EOF'
import ast
import glob
import json
import os
import sys


def load(path):
    with open(path) as f:
        return ast.literal_eval(f.read())


for measure_dir in sys.argv[1:3]:
    ref_tag = sys.argv[3]
    # 12 个主数据文件（scs{方向}{通道}_fix_refbond_...，两字符）
    main = sorted(glob.glob(os.path.join(
        measure_dir, f"scs??_fix_refbond_{ref_tag}.json")))
    assert len(main) == 12, f"{measure_dir}: 主数据文件数 {len(main)} != 12"
    # NOT_APPLICABLE 标记
    with open(os.path.join(measure_dir,
                           "repeat_measurement_scatter.json")) as f:
        mani = json.load(f)
    assert mani["status"] == "NOT_APPLICABLE", f"manifest status: {mani}"
    # spread 文件：行 [[dx,dy],[v,...]]，dx/dy 为整数类
    spread_files = sorted(glob.glob(os.path.join(
        measure_dir, f"scs?_fix_refbond_{ref_tag}_spatial_spread.json")))
    assert spread_files, f"{measure_dir}: 无 spread 文件"
    spread_keys = set()
    for sf in spread_files:
        for row in load(sf):
            key, val = row
            assert isinstance(key, list) and len(key) == 2, f"spread key: {row}"
            dx, dy = key
            assert isinstance(dx, int) and isinstance(dy, int), f"非整数类: {row}"
            assert (isinstance(val, float)
                    or (isinstance(val, list) and len(val) == 2)), \
                f"spread 值: {row}"
            spread_keys.add((dx, dy))
    # geometry 文件：行 [[r0,r1,t0,t1],[7元组]]，与 spread 类集合双向一致
    geo_files = sorted(glob.glob(os.path.join(
        measure_dir, f"scs?_fix_refbond_{ref_tag}_bond_geometry.json")))
    assert len(geo_files) == len(spread_files), \
        f"geometry 文件数 {len(geo_files)} != spread 文件数 {len(spread_files)}"
    geo_classes = set()
    for gf in geo_files:
        for row in load(gf):
            sites, val = row
            assert len(sites) == 4 and len(val) == 7, f"geometry 行: {row}"
            dcol, dx_mid, dy_mid = val[0], val[1], val[2]
            assert isinstance(dcol, int) and isinstance(dx_mid, int) \
                and isinstance(dy_mid, int), f"非整数类: {row}"
            geo_classes.add((dx_mid, dy_mid))
    assert spread_keys == geo_classes, \
        f"spread 类与 geometry 类不一致: {spread_keys ^ geo_classes}"
    print(f"A1 schema OK: {measure_dir}  "
          f"(main=12, spread_classes={len(spread_keys)})")
print("v5 A1 schema check: PASS")
EOF

echo "==== v6: TruncErr 档案（4 份 dmrg 日志）===="
${PYTHON_BIN} "${TRUNC}" \
  "${RUN_DIR}/run_real/dmrg_100.log" \
  "${RUN_DIR}/run_real/dmrg_300.log" \
  "${RUN_DIR}/run_complex/dmrg_100.log" \
  "${RUN_DIR}/run_complex/dmrg_300.log" \
  -o "${OUT}/trunc_err_record.json"
${PYTHON_BIN} - "${OUT}/trunc_err_record.json" <<'EOF'
import json
import sys
with open(sys.argv[1]) as f:
    rec = json.load(f)
for r in rec["records"]:
    n = sum(len(s["sites"]) for s in r["sweeps"])
    assert n > 0, f"{r['log_file']}: 未解析到任何 Site/TruncErr 行"
    print(f"  {r['log_file']}: sweeps={len(r['sweeps'])} "
          f"site_updates={n} E0/site={r['E0_per_site']}")
print("v6 TruncErr check: PASS")
EOF

echo "==== v7: D=100 vs D=300 收敛敏感性（A4，无自动阈值，人审）===="
${PYTHON_BIN} "${BONDDIM}" \
  --ref "${RUN_DIR}/run_real/measure_scs_0_1_D100" \
  --others "${RUN_DIR}/run_real/measure_scs_0_1_D300" \
  -o "${OUT}/bond_dim_convergence_record.json"

echo "SMOKE VERIFY ALL OK"
