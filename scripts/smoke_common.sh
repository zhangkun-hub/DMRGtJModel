# =====================================================================
# 冒烟脚本共享库（纯函数定义，无副作用，可被测试 source）。
# 提供：JSON 配置读取（cfg_get / cfg_list / px）与 fail-closed 报错（die）。
# 注意：本文件不设置 set -e/-u（由调用脚本自行设置），不执行任何命令。
# PYTHON_BIN 可覆盖（集群默认 python3；本地测试可注入 python 真身）。
# =====================================================================

PYTHON_BIN=${PYTHON_BIN:-python3}

# cfg_get <json_file> <key>       -> 输出标量值（布尔转为 true/false）
cfg_get() {
  ${PYTHON_BIN} - "$1" "$2" <<'PY'
import json
import sys
sys.stdout.reconfigure(newline="\n")  # 防止 Windows 文本模式把 \n 写成 \r\n
cfg = json.load(open(sys.argv[1], encoding="utf-8"))
v = cfg[sys.argv[2]]
if isinstance(v, bool):
    print("true" if v else "false")
else:
    print(v)
PY
}

# cfg_list <json_file> <key>      -> 每行输出一个列表元素
cfg_list() {
  ${PYTHON_BIN} - "$1" "$2" <<'PY'
import json
import sys
sys.stdout.reconfigure(newline="\n")
cfg = json.load(open(sys.argv[1], encoding="utf-8"))
print("\n".join(str(x) for x in cfg[sys.argv[2]]))
PY
}

# px <params_json> <key>          -> params.json 标量读取（等价 cfg_get）
px() {
  cfg_get "$1" "$2"
}

# die <message...>                -> stderr 报错并退出 1
die() {
  echo "FAIL: $*" >&2
  exit 1
}

# ensure_new_dir <path>           -> 目录已存在即 FAIL（拒绝覆盖/合并），否则创建
ensure_new_dir() {
  [ ! -e "$1" ] || die "目录已存在，拒绝覆盖/合并：$1"
  mkdir "$1" || die "目录创建失败：$1"
}

# require_files <dir> <file...>    -> 任一文件缺失即 FAIL
require_files() {
  local d=$1
  shift
  local f
  for f in "$@"; do
    [ -f "${d}/${f}" ] || die "文件缺失：${d}/${f}"
  done
}

# check_no_residue                 -> 当前目录除 params.json 外不得有历史 JSON 残留
check_no_residue() {
  local residue
  residue=$(ls *.json 2>/dev/null | grep -v '^params.json$' || true)
  [ -z "${residue}" ] || die "工作目录存在历史输出残留，拒绝继续：${residue}"
}
