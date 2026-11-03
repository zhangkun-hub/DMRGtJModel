"""Local static + synthetic tests for the N=64 smoke scripts.

Covers the 17-item checklist from the safety-review round where feasible
locally (cluster-only items are NOT_VERIFIED by design):

 1  bash -n on run/verify/common
 2  params/config JSON parse
 3  refbond single source of truth (no hardcoded refbond in scripts)
 4  run/verify read the same config via the shared readers
 5  pre-existing run dir -> run FAILs (fail-closed)
 6  real/complex work dirs do not overlap
 7  D=100 -> D=300 chained MPS reuse order preserved (real and complex)
 8  existing diff dir -> FAIL (ensure_new_dir)
 9  residue scs files before measurement -> FAIL (check_no_residue)
 10 missing expected measurement file -> FAIL (require_files /
    smoke_expected_files.py check)
 11 extra measurement file -> FAIL (smoke_expected_files.py check)
 12 no overwrite: no `mv -f`, dirs created via ensure_new_dir
 13 CMakeLists static dependency integrity (test_cmake_sources.py)
 14 verify on incomplete run -> FAIL (RUN_OK missing)
 15 verify with config mismatch -> FAIL
 16 config changed -> readers return the new value (no stale hardcode)
 17 same run dir twice -> second creation fails (exclusivity)

Run: python tools/tests/test_smoke_scripts.py
"""
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, os.pardir, os.pardir))
SCRIPTS = os.path.join(ROOT, "scripts")
RUN_SH = os.path.join(SCRIPTS, "smoke_n64_run.sh")
VERIFY_SH = os.path.join(SCRIPTS, "smoke_n64_verify.sh")
COMMON_SH = os.path.join(SCRIPTS, "smoke_common.sh")
CONFIG = os.path.join(SCRIPTS, "smoke_n64_config.json")
PARAMS = os.path.join(SCRIPTS, "smoke_n64_params.json")
EXPECTED_PY = os.path.join(SCRIPTS, "smoke_expected_files.py")

# Windows 下 python 的 subprocess 会把 "bash" 解析成 WSL 启动器；显式使用
# Git Bash 的 bash.exe，并对传入路径做 posix 化。
BASH_EXE = os.environ.get("GIT_BASH_EXE", r"D:\Git\usr\bin\bash.exe")


def posix(p):
    return p.replace("\\", "/")


def bash_env(extra=None):
    """环境：注入真实 python（Git Bash 的 python3 是商店占位程序），
    并禁用 MSYS 参数路径转换（否则命令串中的 ';' 会被改写）。"""
    env = dict(os.environ)
    env["PYTHON_BIN"] = posix(sys.executable)
    env["MSYS2_ARG_CONV_EXCL"] = "*"
    if extra:
        env.update(extra)
    return env


def read(path):
    with open(path, encoding="utf-8") as f:
        return f.read()


def bash_n(path):
    subprocess.run([BASH_EXE, "-n", posix(path)], check=True)


def run_bash(script, extra_env, script_args=None):
    """Run a script with stubbed `module` (absent locally) and env overrides."""
    env = bash_env(extra_env)
    args = [BASH_EXE, "-c",
            'module(){ :; }; export -f module; exec bash "$@"',
            "bash", posix(script)]
    if script_args:
        args += list(script_args)
    return subprocess.run(args, env=env, capture_output=True, text=True,
                          encoding="utf-8", errors="replace")


def make_fake_tree():
    """Temp tree: scripts/{config,params,common,expected_py} + smoke_n64 base."""
    tree = tempfile.mkdtemp(prefix="smoke_tree_")
    sdir = os.path.join(tree, "scripts")
    os.makedirs(sdir)
    for src in (CONFIG, PARAMS, COMMON_SH):
        shutil.copy(src, os.path.join(sdir, os.path.basename(src)))
    shutil.copy(EXPECTED_PY, os.path.join(sdir, "smoke_expected_files.py"))
    os.makedirs(os.path.join(tree, "smoke_n64"))
    return tree


def test_bash_syntax():
    bash_n(RUN_SH)
    bash_n(VERIFY_SH)
    bash_n(COMMON_SH)
    print("  [PASS] 1. bash -n (run/verify/common)")


def test_json_parse():
    for p in (CONFIG, PARAMS):
        with open(p, encoding="utf-8") as f:
            json.load(f)
    print("  [PASS] 2. params/config JSON parse")


def test_refbond_single_source():
    for script in (RUN_SH, VERIFY_SH):
        text = read(script)
        assert "--refbond=0,1" not in text, f"{script}: hardcoded refbond"
        assert "${SINGLE_REF}" in text, f"{script}: refbond not config-driven"
        assert "${MULTI_REF}" not in text and "reverse_channels" not in text, \
            f"{script}: 残留 B 特性引用"
        assert re.search(r"cfg_get\s+\"\$\{CONFIG_SRC\}\"", text), \
            f"{script}: not reading CONFIG_SRC"
    with open(CONFIG, encoding="utf-8") as f:
        cfg = json.load(f)
    assert cfg["single_refbond"] == "0,1"
    assert "refbonds" not in cfg and "reverse_channels" not in cfg
    print("  [PASS] 3. refbond 唯一事实来源（脚本无硬编码，无 B 残留）")


def test_shared_readers():
    out = subprocess.run(
        [BASH_EXE, "-c",
         f'source "{posix(COMMON_SH)}"; cfg_get "{posix(CONFIG)}" '
         f'single_refbond; '
         f'px "{posix(PARAMS)}" Geometry'],
        capture_output=True, text=True, env=bash_env(), encoding="utf-8", errors="replace")
    lines = out.stdout.split()
    assert lines == ["0,1", "HoneycombXC"], lines
    print("  [PASS] 4. run/verify 共用读取器输出正确")


def test_run_dir_exists_fails():
    tree = make_fake_tree()
    run_id = "fixed01"
    os.makedirs(os.path.join(tree, "smoke_n64", f"run_{run_id}"))
    proc = run_bash(RUN_SH, {"REPO_ROOT": tree, "SMOKE_RUN_ID": run_id})
    assert proc.returncode != 0, proc.stdout
    assert "已存在" in proc.stderr or "已存在" in proc.stdout, proc.stderr
    print("  [PASS] 5. 运行目录已存在 -> run FAIL")


def test_workdirs_disjoint():
    text = read(RUN_SH)
    assert "${RUN_DIR}/run_real" in text and "${RUN_DIR}/run_complex" in text
    assert "build_real" in text and "build_complex" in text
    print("  [PASS] 6. real/complex 工作与构建目录互不重叠")


def test_chain_order():
    text = read(RUN_SH)
    real100 = text.index('RUN_DMRG "${RUN_DIR}/run_real" "${RUN_DIR}/build_real/dmrg" --D=100')
    real300 = text.index('RUN_DMRG "${RUN_DIR}/run_real" "${RUN_DIR}/build_real/dmrg" --D=300')
    cplx100 = text.index('RUN_DMRG "${RUN_DIR}/run_complex" "${RUN_DIR}/build_complex/dmrg" --D=100')
    cplx300 = text.index('RUN_DMRG "${RUN_DIR}/run_complex" "${RUN_DIR}/build_complex/dmrg" --D=300')
    # diff（旧/新二进制）必须在 complex 链之后：旧 build 为 complex 编译
    # （RealCode=OFF），只能读 complex MPS，diff 须与 complex 对配对。
    diff_dir = text.index("run_old_diff")
    assert real100 < real300 < cplx100 < cplx300 < diff_dir
    assert 'cp -r "${RUN_DIR}/run_complex/mps"' in text
    diff_region = text[diff_dir:]
    assert '${RUN_DIR}/build_complex/measureSC_fix_ref_ext' in diff_region
    assert '${RUN_DIR}/build_real/measureSC_fix_ref_ext' not in diff_region
    print("  [PASS] 7. 链顺序与 diff 配对（real -> complex -> diff，complex 对）")


def test_helpers():
    # 8. ensure_new_dir：已存在 FAIL / 新建成功
    tmp = tempfile.mkdtemp(prefix="helpers_")
    ex = os.path.join(tmp, "ex")
    os.makedirs(ex)
    proc = subprocess.run(
        [BASH_EXE, "-c", f'source "{posix(COMMON_SH)}"; ensure_new_dir "{posix(ex)}"'],
        capture_output=True, text=True, env=bash_env(), encoding="utf-8", errors="replace")
    assert proc.returncode != 0 and "已存在" in proc.stderr
    new = os.path.join(tmp, "new")
    subprocess.run(
        [BASH_EXE, "-c", f'source "{posix(COMMON_SH)}"; ensure_new_dir "{posix(new)}"'],
        check=True, capture_output=True, env=bash_env())
    assert os.path.isdir(new)
    # 9. check_no_residue：残留 FAIL / 干净通过
    dirty = os.path.join(tmp, "dirty")
    os.makedirs(dirty)
    with open(os.path.join(dirty, "scsAa_x.json"), "w") as f:
        f.write("[]")
    proc = subprocess.run(
        [BASH_EXE, "-c", f'source "{posix(COMMON_SH)}"; cd "{posix(dirty)}"; check_no_residue'],
        capture_output=True, text=True, env=bash_env(), encoding="utf-8", errors="replace")
    assert proc.returncode != 0 and "残留" in proc.stderr
    clean = os.path.join(tmp, "clean")
    os.makedirs(clean)
    with open(os.path.join(clean, "params.json"), "w") as f:
        f.write("{}")
    subprocess.run(
        [BASH_EXE, "-c", f'source "{posix(COMMON_SH)}"; cd "{posix(clean)}"; check_no_residue'],
        check=True, capture_output=True, env=bash_env())
    # 10. require_files：缺失 FAIL
    proc = subprocess.run(
        [BASH_EXE, "-c", f'source "{posix(COMMON_SH)}"; require_files "{posix(clean)}" nope.json'],
        capture_output=True, text=True, env=bash_env(), encoding="utf-8", errors="replace")
    assert proc.returncode != 0 and "缺失" in proc.stderr
    print("  [PASS] 8/9/10. ensure_new_dir / check_no_residue / require_files")


def test_expected_files():
    # list 模式：单参考键 12 文件（3 方向 × 4 通道）
    out = subprocess.run(
        [sys.executable, EXPECTED_PY, "list", "HoneycombXC", "4", "8",
         "--refbond=0,1"], capture_output=True, text=True, check=True)
    assert len(out.stdout.split()) == 12, out.stdout
    # 多参考键必须被拒绝（C-only）
    proc = subprocess.run(
        [sys.executable, EXPECTED_PY, "list", "HoneycombXC", "4", "8",
         "--refbond=0,1;0,4"],
        capture_output=True, text=True, encoding="utf-8", errors="replace")
    assert proc.returncode != 0 and "单个参考键" in proc.stderr
    # check 模式：完整 -> 0；缺失 -> 非 0；多余 -> 非 0
    tmp = tempfile.mkdtemp(prefix="expected_")
    for f in out.stdout.split():
        with open(os.path.join(tmp, f), "w") as fh:
            fh.write("[]")
    subprocess.run(
        [sys.executable, EXPECTED_PY, "check", tmp, "HoneycombXC", "4", "8",
         "--refbond=0,1"], check=True)
    os.remove(os.path.join(tmp, "scsAa_fix_refbond_0_1.json"))
    proc = subprocess.run(
        [sys.executable, EXPECTED_PY, "check", tmp, "HoneycombXC", "4", "8",
         "--refbond=0,1"],
        capture_output=True, text=True, env=bash_env(), encoding="utf-8", errors="replace")
    assert proc.returncode != 0 and "缺失" in proc.stderr
    # 恢复缺失文件后再测"多余"（checker 先报缺失会掩盖多余，需独立场景）
    with open(os.path.join(tmp, "scsAa_fix_refbond_0_1.json"), "w") as fh:
        fh.write("[]")
    with open(os.path.join(tmp, "scsZz_fix_refbond_9_9.json"), "w") as fh:
        fh.write("[]")
    proc = subprocess.run(
        [sys.executable, EXPECTED_PY, "check", tmp, "HoneycombXC", "4", "8",
         "--refbond=0,1"],
        capture_output=True, text=True, env=bash_env(), encoding="utf-8", errors="replace")
    assert proc.returncode != 0 and "多余" in proc.stderr
    print("  [PASS] 11/12. 期望文件集：缺失/多余均 FAIL，无覆盖")


def test_no_mv_f():
    for script in (RUN_SH, VERIFY_SH):
        # 只检查命令形态（注释中的说明不算）
        assert not re.search(r"^\s*mv\s+-f", read(script), re.M), \
            f"{script}: 存在 mv -f 命令"
    print("  [PASS] 12b. 无 mv -f 命令（逐文件移动）")


def test_verify_preflight():
    tree = make_fake_tree()
    # 14. 无 RUN_OK -> FAIL
    rund = os.path.join(tree, "smoke_n64", "run_x")
    os.makedirs(rund)
    proc = run_bash(VERIFY_SH, {"REPO_ROOT": tree}, [posix(rund)])
    assert proc.returncode != 0 and "RUN_OK" in proc.stderr, proc.stderr
    # 15. RUN_OK 存在但配置快照不一致 -> FAIL
    open(os.path.join(rund, "RUN_OK"), "w").close()
    with open(os.path.join(rund, "smoke_env_record.txt"), "w") as f:
        f.write("x")
    with open(os.path.join(rund, "smoke_config_used.json"), "w") as f:
        json.dump({"single_refbond": "9,9", "refbonds": ["9,9"],
                   "reverse_channels": False}, f)
    proc = run_bash(VERIFY_SH, {"REPO_ROOT": tree}, [posix(rund)])
    assert proc.returncode != 0 and "不一致" in proc.stderr, proc.stderr
    print("  [PASS] 14/15. verify 拒绝未完成/配置不一致的 run")


def test_config_change_propagates():
    tree = make_fake_tree()
    cfg = os.path.join(tree, "scripts", "smoke_n64_config.json")
    with open(cfg, "w", encoding="utf-8") as f:
        json.dump({"single_refbond": "2,3"}, f)
    out = subprocess.run(
        [BASH_EXE, "-c", f'source "{posix(COMMON_SH)}"; cfg_get "{posix(cfg)}" single_refbond'],
        capture_output=True, text=True, env=bash_env(), encoding="utf-8", errors="replace")
    assert out.stdout.split() == ["2,3"], out.stdout
    # 两个脚本都从 CONFIG_SRC 经 cfg_get/cfg_list 取值（无旧值残留）
    for script in (RUN_SH, VERIFY_SH):
        text = read(script)
        assert '"2,3' not in text and '"0,1' not in text.replace(
            'cfg["single_refbond"]', "")
    print("  [PASS] 16. 配置修改后读取新值（无旧硬编码）")


def test_env_switch():
    text = read(RUN_SH)
    assert "SMOKE_ENV=${SMOKE_ENV:-default}" in text
    # intel2022 分支：归档生产口径（Intel 2022.2.0 + icc/icpc）
    assert "/share/intel/2022.2.0" in text
    assert "CMAKE_CC=icc" in text and "CMAKE_CXX=icpc" in text
    # BUILD 编译器已参数化（防 gcc/g++ 硬编码回归）
    assert '-DCMAKE_C_COMPILER="${CMAKE_CC}"' in text
    assert '-DCMAKE_CXX_COMPILER="${CMAKE_CXX}"' in text
    assert "-DCMAKE_C_COMPILER=gcc" not in text
    # default 分支仍走 myenv.sh + gcc 版本校验
    assert "CMAKE_CC=gcc" in text and "CMAKE_CXX=g++" in text
    print("  [PASS] 18. SMOKE_ENV 口径开关（default GNU / intel2022 icpc）")


def test_run_dir_exclusive():
    # 17. 同一 RUN_ID 二次创建失败（fail-closed 互斥）
    tree = make_fake_tree()
    proc1 = run_bash(RUN_SH, {"REPO_ROOT": tree, "SMOKE_RUN_ID": "dup01"})
    # proc1 会因依赖验证失败而退出（无 CMakeCache），但运行目录应已成功创建
    run_dir = os.path.join(tree, "smoke_n64", "run_dup01")
    assert os.path.isdir(run_dir), "首次运行应成功创建运行目录"
    proc2 = run_bash(RUN_SH, {"REPO_ROOT": tree, "SMOKE_RUN_ID": "dup01"})
    assert proc2.returncode != 0 and "已存在" in proc2.stderr, proc2.stderr
    print("  [PASS] 17. 同一运行目录第二次创建 FAIL（互斥）")


def main():
    test_bash_syntax()
    test_json_parse()
    test_refbond_single_source()
    test_shared_readers()
    test_run_dir_exists_fails()
    test_workdirs_disjoint()
    test_chain_order()
    test_helpers()
    test_expected_files()
    test_no_mv_f()
    test_verify_preflight()
    test_config_change_propagates()
    test_env_switch()
    test_run_dir_exclusive()
    subprocess.run([sys.executable,
                    os.path.join(HERE, "test_cmake_sources.py")], check=True)
    print("test_smoke_scripts: ALL PASS")


if __name__ == "__main__":
    main()
