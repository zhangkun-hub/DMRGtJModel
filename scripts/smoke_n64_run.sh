#!/bin/bash
# =====================================================================
# N=64 冒烟 run：环境核对 / real+complex 编译 / dmrg 链式运行 / 测量 / 旧二进制兼容 diff
#
# 安全设计（fail-closed）：
#   - 每次运行创建唯一运行目录 run_<时间戳>（原子 mkdir，已存在即 FAIL，
#     并发提交只有一个能成功）；失败重提交绝不复用旧 MPS/旧输出。
#   - 运行目录内：real/complex 各用独立工作目录；同一次运行内 D=100→D=300
#     有意共用 mps 续跑（A4 链式协议，同参数同扇区）。
#   - 测量前检查 CWD 无历史残留（params.json 除外）；测量后按
#     smoke_expected_files.py（镜像 C++ 枚举的期望文件集）核对缺失/多余，
#     再逐文件移动，无 mv -f、无宽泛通配符收集、无错误吞掉。
#   - 旧 build 只读不写（兼容性 diff 基准；旧 build 为 complex 编译，
#     其 diff 与 complex 链的 MPS/二进制配对）；diff 副本目录创建前必须不存在。
#   - 全部步骤成功后才写入 RUN_OK 标记（verify 前置检查依赖）。
# 用法：sbatch scripts/smoke_n64_run.sh
#       成功后把脚本打印的 RUN_DIR 作为参数交给 smoke_n64_verify.sh。
# 失败条件：任一步骤非零退出（set -euo pipefail）。
# =====================================================================
#SBATCH --job-name=smoke_n64
#SBATCH --partition=256G56c
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=7
#SBATCH --time=2:00:00
#SBATCH --output=smoke_n64_%j.out

set -euo pipefail

# ---------- 路径配置（REPO_ROOT 可被测试覆盖，默认集群路径） ----------
REPO_ROOT=${REPO_ROOT:-/share/home/zhangkun/DMRGtJModel}
OLD_BUILD=${REPO_ROOT}/build                      # 旧 build（含修改前二进制，只读）
SMOKE_BASE=${REPO_ROOT}/smoke_n64                 # 运行目录父目录
CONFIG_SRC=${REPO_ROOT}/scripts/smoke_n64_config.json
PARAMS_SRC=${REPO_ROOT}/scripts/smoke_n64_params.json
COMMON=${REPO_ROOT}/scripts/smoke_common.sh
EXPECTED_PY=${REPO_ROOT}/scripts/smoke_expected_files.py

# shellcheck source=/dev/null
source "${COMMON}"

# ---------- 唯一运行目录（fail-closed，原子创建；先于环境设置，保证失败也有痕迹） ----------
RUN_ID=${SMOKE_RUN_ID:-$(date +%Y%m%d_%H%M%S)}
RUN_DIR=${SMOKE_BASE}/run_${RUN_ID}
mkdir "${SMOKE_BASE}" 2>/dev/null || true   # 父目录可共享
if ! mkdir "${RUN_DIR}"; then
  die "运行目录已存在（重提交或并发冲突）：${RUN_DIR}。请使用全新运行目录，不自动复用/清理旧目录"
fi
echo "RUN_DIR=${RUN_DIR}  （verify 需以此为参数）"
exec > >(tee -a "${RUN_DIR}/step_log.txt") 2>&1

# 配置快照（verify 用其核对 run 与当前配置一致）
cp "${CONFIG_SRC}" "${RUN_DIR}/smoke_config_used.json"

NTASKS=${SLURM_NTASKS:-4}

# ---------- 环境（SMOKE_ENV=default|intel2022，两种口径） ----------
# default：与用户 ${HOME}/myenv.sh 同源（cmake/gcc/jq 模块 + oneAPI 2025
#   setvars），构建用 GNU gcc/g++ 12.1.0，张量运算串行（CMakeLists 在
#   GNU+MKL 下自动关闭 OpenMP）。
# intel2022：归档生产口径（submit_dmrg_zk.sh 同源：GCC 12.1.0 PATH +
#   Intel 2022.2.0 compiler/MKL/MPI vars.sh + anaconda + jq 模块），
#   构建用 icc/icpc，OpenMP 启用（CMakeLists 在 Intel 编译器下自动链接）。
#   用于生产复跑步骤 0 的复值链路 Intel 口径验证。
SMOKE_ENV=${SMOKE_ENV:-default}
case "${SMOKE_ENV}" in
  default)
    module purge 2>/dev/null || true
    if [ -f "${HOME}/myenv.sh" ]; then
      # setvars.sh 在环境已由登录 shell 继承时会"跳过重跑"并可能返回非零；
      # 容忍其返回码，改由下面的 command -v 硬校验兜底（环境真缺失会明确 FAIL）。
      set +e
      source "${HOME}/myenv.sh"
      MYENV_RC=$?
      set -e
      echo "note: source ${HOME}/myenv.sh rc=${MYENV_RC}（rc!=0 但下方环境校验通过时可忽略）"
    else
      die "缺少 ${HOME}/myenv.sh（cmake/gcc/jq 模块 + oneAPI setvars 的集群环境初始化脚本）"
    fi
    command -v cmake >/dev/null 2>&1 || die "cmake 不可用（检查 ${HOME}/myenv.sh）"
    command -v mpirun >/dev/null 2>&1 || die "mpirun 不可用（检查 oneAPI setvars 是否生效）"
    command -v "${PYTHON_BIN}" >/dev/null 2>&1 \
      || die "python3 不可用（请在 ${HOME}/myenv.sh 中补充 python 模块或确认系统 python3）"
    # libhptt.a 是用新版 GCC 编译的（引用 GCC 11+ libstdc++ 符号），必须用
    # myenv.sh 的 gcc 模块（12.1.0）而非系统 cc（8.3.1）构建，否则链接失败。
    GCC_VER_LINE=$(gcc --version 2>/dev/null | head -1 || true)
    GCC_MAJOR=$(echo "${GCC_VER_LINE}" | grep -oE "[0-9]+" | head -1 || true)
    [ -n "${GCC_MAJOR}" ] && [ "${GCC_MAJOR}" -ge 11 ] \
      || die "gcc 版本过旧或不可用（${GCC_VER_LINE:-无}）：需要 >= 11（libhptt.a 工具链匹配）。请确认 myenv.sh 的 gcc/12.1.0 模块已生效"
    CMAKE_CC=gcc
    CMAKE_CXX=g++
    ;;
  intel2022)
    # 与 submit_dmrg_zk.sh 环境段同源（归档生产口径）。
    # vars.sh 会前置覆盖 PATH/LD_LIBRARY_PATH；不做 module purge 之外的
    # oneAPI 2025 变量清理，若两种 setvars 残留冲突会在步骤 0 暴露并修复。
    module purge 2>/dev/null || true
    export PATH=/share/software/gcc/12.1.0/bin:$PATH
    export LD_LIBRARY_PATH=/share/software/gcc/12.1.0/lib64:$LD_LIBRARY_PATH
    # Intel vars.sh 内部脚本（如 oclfpga/fpgavars.sh）会读取未定义变量
    # （OCL_ICD_FILENAMES 等），在 set -u 下报 "unbound variable" 中断；
    # source 期间临时关闭 nounset（保留 -e：source 真失败仍会退出）。
    set +u
    source /share/intel/2022.2.0/compiler/latest/env/vars.sh
    source /share/intel/2022.2.0/mkl/latest/env/vars.sh
    source /share/intel/2022.2.0/mpi/latest/env/vars.sh
    source /share/software/anaconda3/bin/activate
    set -u
    module use /share/software/modulefiles 2>/dev/null || true
    module load jq 2>/dev/null || true
    command -v icpc >/dev/null 2>&1 || die "icpc 不可用（检查 Intel 2022.2.0 compiler vars.sh）"
    command -v mpirun >/dev/null 2>&1 || die "mpirun 不可用（检查 Intel 2022.2.0 MPI vars.sh）"
    command -v "${PYTHON_BIN}" >/dev/null 2>&1 \
      || die "python3 不可用（检查 anaconda 激活是否生效）"
    CMAKE_CC=icc
    CMAKE_CXX=icpc
    ;;
  *)
    die "未知 SMOKE_ENV=${SMOKE_ENV}（支持 default|intel2022）"
    ;;
esac

# ---------- step 0: 环境记录与依赖验证（先读后记，失败信息明确） ----------
ENV_REC="${RUN_DIR}/smoke_env_record.txt"
{
  echo "run_dir=${RUN_DIR}"
  echo "smoke_env=${SMOKE_ENV}"
  echo "repo_head=$(git -C "${REPO_ROOT}" rev-parse HEAD 2>/dev/null || echo UNKNOWN)"
  module list 2>&1
  icpc --version 2>&1 || icpx --version 2>&1 || true
  echo "cxx=$(c++ --version 2>&1 | head -1 || true)"
  echo "cmake_version=$(cmake --version 2>&1 | head -1 || true)"
} > "${ENV_REC}"

CACHE_FILE="${OLD_BUILD}/CMakeCache.txt"
[ -f "${CACHE_FILE}" ] || die "旧 build 的 CMakeCache.txt 不存在：${CACHE_FILE}"
CACHE_DEPS=$(grep -E "TENSOR_HEADER_PATH|MPS_HEADER_PATH" "${CACHE_FILE}" || true)
[ -n "${CACHE_DEPS}" ] || die "${CACHE_FILE} 缺少 TENSOR_HEADER_PATH / MPS_HEADER_PATH 条目"
echo "${CACHE_DEPS}" | tee -a "${ENV_REC}" >/dev/null

# 优先使用缓存中的 CMAKE_PREFIX_PATH（若曾被 -D 记录），否则回退到头文件路径
CACHED_PREFIX=$(grep "^CMAKE_PREFIX_PATH" "${CACHE_FILE}" | head -1 | cut -d= -f2 || true)
HEADER_PATH=$(grep "TENSOR_HEADER_PATH:PATH" "${CACHE_FILE}" | head -1 | cut -d= -f2 || true)
PREFIX="${CACHED_PREFIX:-${HEADER_PATH}}"
[ -n "${PREFIX}" ] || die "无法从 ${CACHE_FILE} 确定依赖前缀（CMAKE_PREFIX_PATH / TENSOR_HEADER_PATH 均缺失）"
if [ ! -d "${PREFIX}/qlten" ] || [ ! -d "${PREFIX}/qlmps" ]; then
  die "依赖头目录缺失：需要 ${PREFIX}/qlten 与 ${PREFIX}/qlmps（请核对集群实际路径后修正缓存提取）"
fi
echo "CMAKE_PREFIX_PATH=${PREFIX}" | tee -a "${ENV_REC}" >/dev/null

# hptt 路径同样从旧缓存提取（新构建目录无缓存，需显式传入 cmake）
# （过滤 -ADVANCED 等内部缓存变量，只取真实路径条目）
HPTT_FLAGS=()
while IFS= read -r line; do
  key=${line%%=*}
  key=${key%%:*}
  val=${line#*=}
  HPTT_FLAGS+=("-D${key}=${val}")
done < <(grep -E "^hptt_" "${CACHE_FILE}" | grep -v -- "-ADVANCED" || true)
echo "hptt_cache_entries=${HPTT_FLAGS[*]:-（无）}" | tee -a "${ENV_REC}" >/dev/null

# ---------- 配置（唯一事实来源：smoke_n64_config.json + params.json） ----------
GEOM=$(px "${PARAMS_SRC}" Geometry)
LY=$(px "${PARAMS_SRC}" Ly)
LX=$(px "${PARAMS_SRC}" Lx)
SINGLE_REF=$(cfg_get "${CONFIG_SRC}" single_refbond)
N_SITES=$(( 2 * LY * LX ))

# ---------- 构建（独立目录，fail-closed 创建） ----------
# 显式用 gcc/g++（myenv 模块 12.1.0，与 libhptt.a 工具链匹配）；
# 只构建冒烟所需目标（plot_lattice 的 std::filesystem 链接问题与本任务
# 无关，跳过全目标构建）；不启用 OpenMP（qlten 头文件 pragma 写法在
# GNU 下会报错，pragma 被忽略即张量运算串行，只影响性能不影响正确性）。
BUILD(){ # BUILD <rel_dir> [cmake flags...]
  local dir=$1; shift
  ensure_new_dir "${RUN_DIR}/${dir}"
  cd "${RUN_DIR}/${dir}"
  cmake "${REPO_ROOT}" -DCMAKE_PREFIX_PATH="${PREFIX}" \
    -DCMAKE_C_COMPILER="${CMAKE_CC}" -DCMAKE_CXX_COMPILER="${CMAKE_CXX}" \
    "${HPTT_FLAGS[@]}" "$@" \
    | tee cmake.log
  make -j8 dmrg measure1 measureSC_fix_ref_ext | tee make.log
}

echo "==== step 1: real build 编译 ===="
BUILD build_real -DRealCode=ON
echo "real build OK: ${RUN_DIR}/build_real"

echo "==== step 2: complex build 编译（默认 RealCode OFF = USE_COMPLEX）===="
BUILD build_complex
echo "complex build OK: ${RUN_DIR}/build_complex"

# ---------- 计算与测量 ----------
RUN_DMRG(){ # RUN_DMRG <workdir> <dmrg_binary> <--D=...>
  local wd=$1; local bin=$2; local dflag=$3
  # 工作目录可已存在：同一次运行内 D=100→D=300 有意共用 mps 续跑；
  # 跨运行复用由唯一运行目录机制阻止。
  mkdir -p "${wd}" || die "工作目录创建失败：${wd}"
  cd "${wd}"
  cp "${PARAMS_SRC}" params.json
  mpirun -n ${NTASKS} "${bin}" params.json "$dflag" \
    | tee "dmrg_${dflag#--D=}.log"
}
MEAS1(){ # MEAS1 <workdir> <measure1_binary> <tag>
  local wd=$1; local bin=$2; local tag=$3
  cd "${wd}"
  [ -d mps ] || die "mps 目录缺失：${wd}"
  check_no_residue
  mpirun -n ${NTASKS} "${bin}" params.json 0
  require_files "${wd}" nf.json sz.json sz0sz.json sp0sm.json sm0sp.json \
    nf0nf.json single_particle_corr1.json single_particle_corr2.json \
    single_particle_corr3.json single_particle_corr4.json
  ensure_new_dir "measure1_${tag}"
  mv nf.json sz.json sz0sz.json sp0sm.json sm0sp.json nf0nf.json \
     single_particle_corr1.json single_particle_corr2.json \
     single_particle_corr3.json single_particle_corr4.json \
     "measure1_${tag}/"
}
MEAS_SCS(){ # MEAS_SCS <workdir> <measureSC_binary> <tag> [args...]
  local wd=$1; local bin=$2; local tag=$3; shift 3
  cd "${wd}"
  [ -d mps ] || die "mps 目录缺失：${wd}"
  check_no_residue
  mpirun -n ${NTASKS} "${bin}" params.json "$@"
  # 期望文件集：镜像 C++ 枚举（几何+过滤+方向+通道），与实测数据无关
  mapfile -t EXPECTED < <(${PYTHON_BIN} "${EXPECTED_PY}" list "${GEOM}" "${LY}" "${LX}" "$@")
  for f in "${EXPECTED[@]}"; do
    [ -f "$f" ] || die "测量输出缺失：${wd}/${f}"
  done
  local unexpected
  unexpected=$(comm -23 \
      <(ls scs*.json 2>/dev/null | sort) \
      <(printf '%s\n' "${EXPECTED[@]}" | sort) || true)
  [ -z "${unexpected}" ] || die "存在不在期望集合中的输出：${unexpected}"
  ensure_new_dir "measure_scs_${tag}"
  mv "${EXPECTED[@]}" "measure_scs_${tag}/"
}

echo "==== step 3: real dmrg D=100 -> 测量 -> D=300（同目录有意续跑）-> 测量 ===="
RUN_DMRG "${RUN_DIR}/run_real" "${RUN_DIR}/build_real/dmrg" --D=100
MEAS1 "${RUN_DIR}/run_real" "${RUN_DIR}/build_real/measure1" D100
MEAS_SCS "${RUN_DIR}/run_real" "${RUN_DIR}/build_real/measureSC_fix_ref_ext" \
        0_1_D100 --refbond="${SINGLE_REF}"
RUN_DMRG "${RUN_DIR}/run_real" "${RUN_DIR}/build_real/dmrg" --D=300
MEAS1 "${RUN_DIR}/run_real" "${RUN_DIR}/build_real/measure1" D300
MEAS_SCS "${RUN_DIR}/run_real" "${RUN_DIR}/build_real/measureSC_fix_ref_ext" \
        0_1_D300 --refbond="${SINGLE_REF}"

echo "==== step 4: complex dmrg D=100,300 + 测量（实初态，虚部应≈0）===="
RUN_DMRG "${RUN_DIR}/run_complex" "${RUN_DIR}/build_complex/dmrg" --D=100
RUN_DMRG "${RUN_DIR}/run_complex" "${RUN_DIR}/build_complex/dmrg" --D=300
MEAS1 "${RUN_DIR}/run_complex" "${RUN_DIR}/build_complex/measure1" D300
MEAS_SCS "${RUN_DIR}/run_complex" "${RUN_DIR}/build_complex/measureSC_fix_ref_ext" \
        0_1_D300 --refbond="${SINGLE_REF}"

echo "==== step 5: 旧二进制 vs 新二进制 单参考键字节级 diff（complex 对）===="
# 旧 build 为 complex 编译（CMakeCache RealCode=OFF，2026-06-18），只能读
# complex MPS；故 diff 在 complex 链完成后执行，旧/新均为 complex 二进制，
# 12 文件输出同为 [re,im] 格式，字节级比较有效。
for d in run_old_diff run_new_diff; do
  ensure_new_dir "${RUN_DIR}/${d}"
done
for d in run_old_diff run_new_diff; do
  cp -r "${RUN_DIR}/run_complex/mps" "${RUN_DIR}/${d}/mps"
  SRC_CNT=$(ls "${RUN_DIR}/run_complex/mps"/mps_ten*.qlten 2>/dev/null | wc -l)
  DST_CNT=$(ls "${RUN_DIR}/${d}/mps"/mps_ten*.qlten 2>/dev/null | wc -l)
  { [ "${SRC_CNT}" -ge "${N_SITES}" ] && [ "${DST_CNT}" -eq "${SRC_CNT}" ]; } \
    || die "mps 复制不完整：${d}（src=${SRC_CNT} dst=${DST_CNT} N=${N_SITES}）"
done
(cd "${RUN_DIR}/run_old_diff" && cp "${PARAMS_SRC}" params.json \
  && mpirun -n ${NTASKS} "${OLD_BUILD}/measureSC_fix_ref_ext" params.json --refbond="${SINGLE_REF}")
(cd "${RUN_DIR}/run_new_diff" && cp "${PARAMS_SRC}" params.json \
  && mpirun -n ${NTASKS} "${RUN_DIR}/build_complex/measureSC_fix_ref_ext" params.json --refbond="${SINGLE_REF}")
DIFF_CNT=0
for f in scsA{a,b,c,d} scsB{a,b,c,d} scsC{a,b,c,d}; do
  fn="${f}_fix_refbond_$(echo "${SINGLE_REF}" | tr ',' '_').json"
  if ! cmp -s "${RUN_DIR}/run_old_diff/${fn}" "${RUN_DIR}/run_new_diff/${fn}"; then
    echo "FAIL: 12 文件兼容 diff 失败: ${fn}"
    DIFF_CNT=$((DIFF_CNT+1))
  fi
done
if [ ${DIFF_CNT} -ne 0 ]; then die "${DIFF_CNT} 个文件与旧版不一致"; fi
echo "old-vs-new 12 文件 diff: 全部一致"

echo "==== step 6: 完成标记（最后写入，失败路径不会产生）===="
touch "${RUN_DIR}/RUN_OK"
echo "smoke_n64_run: ALL STEPS OK"
echo "RUN_DIR=${RUN_DIR}"
echo "下一步：sbatch scripts/smoke_n64_verify.sh ${RUN_DIR}"
