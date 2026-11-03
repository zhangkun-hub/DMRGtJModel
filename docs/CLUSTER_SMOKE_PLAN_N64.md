# 集群 N=64 冒烟方案（C-only，待批准，不自动提交）

> 对象：任务 C（复值基态与复数配对关联）的集群验收冒烟。
> 脚本：`scripts/smoke_n64_run.sh`（运行+编译+测量）、
> `scripts/smoke_n64_verify.sh`（验证）。
> 几何：HoneycombXC Ly=4 Lx=8（N=64），参数 `scripts/smoke_n64_params.json`
> （t=0.3/t2=0.01/J=0.1/J2=0.01、NumHole=2、Sweeps=4、noise 数组同生产）。
> 参考键唯一事实来源：`scripts/smoke_n64_config.json`（single_refbond="0,1"）。
> 总资源：1 节点 256G56c，run 作业 ≤2h（实际预计 <30 min），verify ≤30 min。

## 安全设计

- **唯一运行目录**：每次 run 创建 `smoke_n64/run_<时间戳>`（原子 mkdir，
  已存在即 FAIL；并发提交只有一个成功）。失败重提交**绝不复用**旧 MPS/旧
  输出；不自动清理旧目录（无 `rm -rf`）。
- **MPS 隔离**：real/complex 各用独立工作目录；同一次运行内 D=100→D=300
  有意共用 mps 续跑；diff 副本目录创建前必须不存在，复制后核对 mps 张量
  文件数量。
- **测量完整性**：测量前检查工作目录无历史 JSON 残留；测量后按
  `scripts/smoke_expected_files.py`（镜像基线 C++ 枚举：几何+过滤 j0>i1+
  方向+通道 a-d，单参考键）核对缺失/多余，再逐文件移动（无 `mv -f`、
  无宽泛通配符、无错误吞掉）。
- **旧二进制保护**：旧 build 目录只读（兼容性 diff 基准），新构建目录与其
  完全分离。
- **完成标记**：run 全部成功后才写 `RUN_OK`；verify 前置检查 RUN_OK、
  环境记录、配置快照一致性（`cmp`）、必要结果目录与日志，任一缺失即 FAIL。
- **配置同源**：参考键只存在于 `smoke_n64_config.json`，run 与 verify 均经
  `scripts/smoke_common.sh` 读取。

## 提交方式（不自动执行）

```bash
cd /share/home/zhangkun/DMRGtJModel
sbatch scripts/smoke_n64_run.sh
# run 作业日志末尾打印 RUN_DIR=.../smoke_n64/run_<时间戳>；随后：
sbatch scripts/smoke_n64_verify.sh /share/home/zhangkun/DMRGtJModel/smoke_n64/run_<时间戳>
```

verify 无 `--dependency` 前置（手动顺序提交）；若 run 未完成/失败，verify
会因 RUN_OK 或结果缺失而**显式 FAIL**，不会误 PASS。

## 步骤与验收

| # | 步骤 | 输入 | 输出 | 失败条件 | 资源 |
|---|---|---|---|---|---|
| 1 | 环境核对 | `source ${HOME}/myenv.sh`（**与用户日常环境同源**：cmake/4.0.0-rc1、gcc/12.1.0、jq/1.7.1 模块 + oneAPI setvars）、`OLD_BUILD/CMakeCache.txt` | `smoke_env_record.txt`（run_dir、repo HEAD、module 列表、icpc/icpx 与 cmake 版本、依赖路径、hptt 缓存条目）；依赖前缀必须同时含 `qlten/` 与 `qlmps/` 目录 | myenv.sh 缺失、cmake/mpirun/python3 不可用、缓存条目缺失或依赖目录缺失 → 明确 FAIL（fail-fast） | 无 |
| 2 | real build（C1 已修复） | 仓库源码 + 前缀 + `-DRealCode=ON` | `build_real/`、cmake/make 日志 | 编译错误 → 退出 1（CMakeLists 失效 target 已移除，全部 target 源文件存在性经静态测试核对） | ~10 min |
| 3 | complex build | 同上，默认选项（`USE_COMPLEX`） | `build_complex/` | 同上 | ~10 min |
| 4 | 旧/新二进制单参考键字节级 diff | 同一 `run_real/mps`（D=300 后）复制到 run_old_diff/run_new_diff（复制前必须不存在，复制后数量核对）；旧/新二进制同参测 `--refbond=0,1` | 两份 12 文件 | 任一 `scs{X}{a..d}_fix_refbond_0_1.json` `cmp` 不一致 → 退出 1 | ~1 min |
| 5 | real dmrg 链 D=100→测量→D=300→测量 | `run_real/` | `dmrg_100/300.log`；`measure1_D{100,300}/`、`measure_scs_0_1_D{100,300}/`（各 12 主文件） | 任一步骤非零退出；测量后期望文件集核对不过 → 退出 | 4×7 核 ~5 min/段 |
| 6 | complex dmrg 链 + 测量（实初态，虚部应≈0） | `run_complex/` | 同第 5 步复值版（值 `[re, im]`） | 同上 | 同第 5 步 |
| 7 | 复值重现实值（验收 1/2） | v1：scs 单参考键目录对 + 两日志（独立 manifest）；v2：measure1 目录对（generic 模式） | `verify_<时间>/verify_{scs_0_1,measure1}.json` | 完整性任一失败或数值超容差（atol=rtol=1e-8 组合判据，实际容差入报告）→ FAIL；若仅数值超差，检查 scaled 值再定论 | 单核 <1 min |

## 全部通过后的判定

- run：`smoke_n64_run: ALL STEPS OK` + `RUN_OK` 标记存在；
- verify：`SMOKE VERIFY ALL OK` → 任务 C 的复值计算链路（实值态重现 +
  虚部≈0）验收通过，具备进入复值生产复跑方案拟定的条件；
- 任一 FAIL → 保留 `smoke_n64/run_<时间戳>/` 全部产物与 `step_log.txt`，
  回报后决定修复与重跑范围，不擅自放宽阈值或跳过检查。

## 结果归档

- `smoke_env_record.txt` 回填 `docs/PROTOCOL_MEASURE_SC_LRO_V1.md` 集群信息栏；
- 各 `verify_<时间>/` 报告与两份作业日志作为冒烟验收记录随协议保存；
- 冒烟数据不进入归档数据目录（归档目录一律不动）。

## 已知限制（如实声明）

- 集群依赖真实环境项（module/编译器/qlten/qlmps 路径与版本、编译通过性、
  旧二进制 diff、DMRG 数值、端到端验证）在本地均为 **NOT_VERIFIED**，以
  冒烟作业为准；
- 冒烟只构建所需目标（dmrg/measure1/measureSC_fix_ref_ext）：`plot_lattice`
  在 GCC 8 下因 `std::filesystem` 缺少 `-lstdc++fs` 链接失败（与本任务
  无关的既有问题，待单独修复）；
- 集群实际编译器为系统 GNU 8.3.1；**不启用 OpenMP**（GNU 8 拒绝 qlten
  头文件的 OpenMP shared 子句写法），张量运算串行——只影响性能不影响
  正确性，N=64 冒烟可接受；生产协议若要复现归档数据的 Intel 编译器口径
  （含 OpenMP），需单独评估；
- `dmrg.cpp` 按 mps 文件数复用 MPS（不校验扇区）——本方案通过唯一运行目录
  与 fail-closed 目录创建规避跨运行复用；同一次运行内 D100→D300 复用为
  有意设计；
- 本冒烟只验证**复值计算能力在实值态上的数值一致性**；不验证手性态选择、
  不构成 d+id 证据（见 `docs/COMPLEX_RUN_PROTOCOL.md` 的手性态选择边界）；
- verify 的 generic 模式（measure1 目录）无独立 manifest（两侧同缺同漏的
  理论可能），scs 数据使用独立 manifest 的 scs-mode；
- 任务 A/B 已还原：本冒烟不再包含误差记录、共轭检查、多参考键等步骤。
