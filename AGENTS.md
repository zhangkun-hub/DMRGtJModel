# AGENTS.md

本文件为 AI 编码助手提供项目指令，作用域覆盖仓库全部目录。

## 项目概述

DMRG 基态求解器，适用于正方晶格与蜂窝晶格 t-J 模型（含 NNN 项 t2/J2），基于两个 header-only 库：
- **qlten** (TensorToolkit) — QN 守恒张量，`find_path(qlten)` 定位
- **qlmps** (UltraDMRG) — MPS/MPO/DMRG 算法，`find_path(qlmps)` 定位

典型工作流：生成参数 → 集群跑 DMRG → 测量关联 → MATLAB 绘图（`plot/`）。

## 构建

本地 Windows 为工作副本，编译与计算在 SUSTech 集群（`/share/home/zhangkun/DMRGtJModel/build`，Intel oneAPI + MKL + MPI）。

```bash
mkdir build && cd build
cmake .. -DCMAKE_PREFIX_PATH=<path-to-qlten-and-qlmps>
make -j4
```

### CMake 选项（默认值）

| 选项 | 默认 | 说明 |
|---|---|---|
| `HIGHER_SYMMETRY` | ON | U(1)×U(1) 粒子数守恒（`SYM_LEVEL=1`）；OFF 时 `mu`/`BondSingletPairPerturbation` 生效 |
| `RealCode` | OFF | OFF=复数张量，ON=实数；非零 `phi` 必须复数构建 |
| `QLMPS_USE_GPU` | OFF | 单 GPU CUDA/cuTENSOR，不支持 MPI |
| `COMPILE_FIX_CODE` | OFF | 编译 `fix_*`/`continue` 调试工具 |
| `QLTEN_TIMING_MODE` / `QLMPS_TIMING_MODE` | OFF | 分部计时输出 |

### 集群构建要点

- 登录 shell 需先 `source ${HOME}/myenv.sh`（cmake/4.0.0-rc1 + gcc/12.1.0 + jq/1.7.1 + oneAPI setvars；setvars 返回 rc=3 可忽略）
- 必须 `-DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++`（系统 cc 8.3.1 与 libhptt.a ABI 不匹配）
- 不启用 OpenMP（qlten `#pragma omp` 在 GNU 下编译失败，串行仅影响性能）
- 集群 build 为复数构建（`RealCode:BOOL=OFF`）

## 二进制程序

| Binary | 用途 | 备注 |
|---|---|---|
| `dmrg` | 传统 DMRG | 内存构建 MPO，每次运行重建 |
| `mpogen` | 导出 MPO 到 `./mpo/` | `vmps`/`vmps1FM`/`noised_vmps` 的前置 |
| `vmps` | MPO-based DMRG | 从 `./mpo` 加载 |
| `noised_vmps` | 带噪声的 MPO-based DMRG | 使用 JSON 的 `noise` 字段 |
| `vmps1FM` | 全极化初态 DMRG | |
| `measure1` | 单点量+两点关联 | 可选第 2 参数指定源点（默认 0） |
| `measureSC_fix_ref_ext` | 单态/三重态配对关联 | `--refbond=i,j` 或 `--start=i` |
| `measureS` | 总自旋（仅 CPU） | |
| `move_center` | MPS 中心移至 site 0 | |
| `plot_lattice` | 导出晶格 SVG | GNU 下不链接（`std::filesystem` 需 `-lstdc++fs`） |

运行示例：`mpirun -n 8 ./dmrg params.json --D=100,500,1000`

`--D=` 逗号分隔、无空格、升序，程序逐值跑 sweep。`dmrg`/`vmps` 中 `MaxLanczIter` 沿 D-list 线性递增。未指定时使用 JSON 的 Dmin/Dmax/CutOff。

## 运行时数据布局

所有程序在工作目录下读写：
- `mps/` — 基态 MPS 张量。已有合法 `mps/` 时自动加载续算
- `mpo/` — `mpogen` 导出的 MPO 张量
- `.temp/` — 有限 VMPS 环境张量，可达 TB 级，优化中途不可删除

测量程序输出 `*.json` 文件到工作目录。归档命名见 `docs/DATA_NAMING.md`。

## 架构

- 类型骨架 `src/tJ_hilbert_space.h`：`TenElemT`（double/complex）、`QNT`（U1U1QN/U1QN）、全局 `tj_site`
- 参数 `src/params_case.h`：`CaseParams` 解析 `params.json`；`NormalizeGeometry` 映射旧名到 `Square*`；`tJModelParamters` 紧凑参数包
- 晶格 `src/lattice_base.h`、`src/squarelattice.h`、`src/honeycomb_lattice.h`：`LatticeBase` 持有 NN/NNN/3NN 链接表及 direct/wind 分类；格点 column-major `site = x*Ly + y`；蜂窝 `N = 2*Lx*Ly`，A/B 子晶格 = `(x+y)%2`
- MPO 构建 `src/tJmodel.h`：`ConstructtJMPOGenerator<LatticeT>` 从链接表添加项，`phi≠0` 时加缠绕相位
- 几何分派链（`if Ly==2 → SquareLadder` ...）在 `dmrg.cpp`/`mpogen.cpp` 中重复；添加几何须更新所有分派位置
- 驱动模式：`MPI_Init` → `CaseParams` → 建晶格 → MPO → 初始化 MPS（Néel + 均匀空穴） → sweep → dump
- 工具 `src/myutil.cpp`：`ParserBondDimension`、`ComputeTotalSites`、`IsHoneycombGeometry`
- 测量 `src/my_measure.h`：qlmps 测量原语封装，MPI round-robin 分配

## 计算工作流

### 项目分工

| 任务 | 项目 |
|---|---|
| DMRG 计算、绝热演化、测量 | 本仓库 + 集群 |
| 数据提取、作图、拟合、相态判断 | `DMRG-data-processing`（有 `dmrg-density-plot`/`dmrg-spin-corr-plot`/`dmrg-supercorr-plot` 三个 skill） |
| UltraDMRG 库算法/接口开发 | `UltraDMRG`（有 `ultradmrg-algorithm-dev` skill） |

### 环境与同步

- 本地改代码后 `git commit` + `git push`，集群 `git pull` 后 `sbatch`
- 每个计算独立 run 目录：含 `params.json`、`mps/`、`mpo/`、`.temp/` 与 SLURM 日志
- 参数格式见 `docs/parameters.md`；根目录 `params.json` 是当前活跃示例

### MPS 续跑与 .temp 规则

- run 目录已有 `mps/` 且张量数匹配时，程序自动加载旧 MPS 续算；否则从 Néel 态 + 均匀空穴初始化
- `.temp/` 是可重建的临时环境张量，非续跑/测量必需
- **同一物理参数**（加大 D、重提交、测量）：`.temp` 可删
- **换物理参数绝热演化**：切换参数前**必须删除** `.temp`
- 删前核对路径，勿误删 `mps/` 或测量结果

### 绝热演化

从"容易算对"的邻近参数出发，逐步演化到目标参数，每步确认收敛后再走下一步。直接在新参数从 Néel 初态出发容易陷入亚稳态。

**典型路径（honeycomb t-J，目标 t=0.3 J=0.1 纯 t-J）**：
1. `t=0.3, J=0.1, t2=0.01, J2=0.01` → `sbatch` 跑收敛
2. 删除 `.temp/`
3. 改 `params.json` 为 `t2=0, J2=0`，同目录再次 `sbatch`（自动加载上一步 MPS）
4. 检查收敛：能量/截断误差/纠缠熵随 sweep 稳定；不收敛则插入中间参数点

### 提交 DMRG（submit_dmrg_zk.sh）

- 改脚本中 `D_VALUES_DMRG`（升序、逗号分隔、无空格，如 `"100,500,800,1000,1200,1500"`）
- `sbatch submit_dmrg_zk.sh`（默认 8 MPI × 7 threads）
- 收敛判据：日志中 `E0`/`TruncErr`/`S` 随 sweep 稳定

### 测量

**measure_corr.sh（单点+两点关联）**
- 改脚本中 `sites` 列表（两点关联源点）；单点量 `nf.json`/`sz.json` 与源点无关
- 源点选在离边界约 1/4 处（N=200 → 49/50，N=400 → 97/98）
- 输出移入 `${bond}_measure_phys_qty/measure_${site}/`

**measure_super_corr.sh（四点关联 → 超导关联）**
- 改 `refbonds` 列表（格式 `"i,j"`），选离边界 1/4 处
- 输出 12 个文件 `scs{X}{ch}_fix_refbond_{i}_{j}.json`，移入 `${bond}_measure_super_corr/super_corr_{i}_{j}/`
- 超导关联合成 `<Δ†Δ> = 0.5(a-b-c+d)` 由后处理项目完成

### 数据交接

后处理项目按 `Data_N{200|400}/{t_J}/{doping}/{bond_dim}_measure_*/` 读取，集群 run 目录命名应与此对应。

## 复值协议（d+id 探索）

- 默认构建即复值（`USE_COMPLEX`）；测量输出 `[re, im]`
- 复值二进制读不了实值 MPS，反之亦然；real/complex 各用独立 build 与 run 目录
- 验收：N=64 冒烟 `scripts/smoke_n64_run.sh` + `scripts/smoke_n64_verify.sh`（2026-11-03 已通过）
- 工具：`tools/verify_complex_reproduces_real.py`（scs 硬判据 / generic re-diff 报告）、`tools/diag_measure1_diff.py`（逐行排查）
- 复值能跑 ≠ 手性态已被选择 ≠ d+id 证据；完整协议见 `docs/COMPLEX_RUN_PROTOCOL.md`

## 已知约束

- 蜂窝晶格忽略 `phi`（强制置 0）；非零 `phi` 需复数构建
- `vmps*` 系列必须先跑 `mpogen` 生成 `./mpo`
- 换系统大小时旧 MPS 被覆盖（张量数不匹配自动重新初始化），需保留先拷贝
- 根目录与 run 目录的 `params.json` 不是一回事，SLURM 脚本用 run 目录里那份
- 蜂窝晶格 `N = 2*Lx*Ly`，`Lx`/`Ly` 语义与正方晶格不同（2 格点/原胞）
- `plot_lattice` 在 GNU 下不链接（`std::filesystem` 需 `-lstdc++fs`，已知问题）
- 几何分派链在 `dmrg.cpp`/`mpogen.cpp` 中重复，新增几何须更新所有分派位置

## 用户角色

项目维护者是凝聚态计算物理方向的硕士研究生（研三），项目目标是完成论文和毕业。讨论科研问题、给出技术方案时，应站在硕士层级沟通，表述清晰易懂，不做超出该层级的过度深入或抽象展开。

## 关联项目查阅规则

普通参数配置、提交计算、测量、数据交接等流程问题，优先查阅本仓库文档与脚本。只有确实涉及 UltraDMRG 内部接口或算法实现、且现有文档不足以判断时，才扩大到 UltraDMRG 源码。
