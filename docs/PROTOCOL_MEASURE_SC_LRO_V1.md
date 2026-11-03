# Measure SC-LRO 协议锁定清单 v1

> 本文件是 `measure_sc_lro_protocol_v1` 协议的本地可追溯记录。协议 tag 已创建于：
> - DMRGtJModel：`measure_sc_lro_protocol_v1` @ `431aa08`
> - UltraDMRG（本地审计参照）：`measure_sc_lro_protocol_v1` @ `2ac051f`
> 随数据保存时，本清单与 tag 哈希一并归档。

## 仓库与版本

| 仓库 | commit | 状态 |
|---|---|---|
| DMRGtJModel | `431aa08`（"Add project CLAUDE.md, dmrgtj-calculation skill, ..."） | 工作区干净 |
| UltraDMRG（本地） | `2ac051f`（fork of QuantumLiquids/UltraDMRG） | 工作区干净 |
| DMRG-data-processing（消费侧） | `418be2c` | 3 个未跟踪文件 |

注意：实际编译使用的 qlmps 安装在集群，与本机 UltraDMRG **版本不一致**
（`tJSite`/`tJOperators`/`MeasureElectronPhonon4PointFunction` 仅存在于集群
qlmps）。集群侧信息需回填（见下）。

## 编译配置

- `RealCode`：归档数据 `ON`（`-DUSE_REAL=1`）；复值协议用默认 `OFF`
  （`-DUSE_COMPLEX=1`）。
- `HIGHER_SYMMETRY`：`ON`（SYM_LEVEL=1，U1U1QN 守恒 N、Sz）。
- `QLMPS_USE_GPU`：`OFF`；`COMPILE_FIX_CODE`：`OFF`；timing 选项全 `OFF`。
- 依赖：BLAS/LAPACK（MKL）、MPI、OpenMP、hptt。

## 集群环境（已回填：N=64 冒烟，2026-11-03，smoke_env_record.txt）

- 环境入口：`source ${HOME}/myenv.sh`（cmake/4.0.0-rc1、gcc/12.1.0、jq/1.7.1
  模块 + Intel oneAPI `setvars.sh`；冒烟脚本与其同源）。
- qlten 安装路径与版本：`/share/home/zhangkun/.local/include/qlten`
  （header-only；qlten.h 时间戳 2025-10-17，无版本号）。
- qlmps 安装路径与版本：`/share/home/zhangkun/.local/include/qlmps`
  （header-only，集群安装；版本号未记录，与本地 UltraDMRG 不一致，见上）。
- module 列表（`module list`，冒烟时记录）：
  `1) cmake/4.0.0-rc1  2) gcc/12.1.0  3) jq/1.7.1`
- 编译器版本：`icpc` 不存在（oneAPI 2025 已移除经典 icpc）；`icpx` =
  Intel oneAPI DPC++/C++ Compiler 2025.2.1（2025.2.0.20250806）。
  冒烟构建显式指定 `gcc/g++ 12.1.0`（与 libhptt.a 工具链匹配）。
- python3 来源：系统 `/usr/bin/python3`（Python 3.6.8）。
- 其他依赖：hptt `~/.local/include` + `~/.local/lib/libhptt.a`；
  CMAKE_PREFIX_PATH=`/share/home/zhangkun/.local/include`。

## 复值链路冒烟验收记录（2026-11-03，N=64 HoneycombXC 8×4）

- 运行目录：`smoke_n64/run_20261103_224728/`（run 作业 `ALL STEPS OK`，
  repo_head `f065209`；verify 重跑时 repo_head `5090aef`）。
- 编译：real（RealCode=ON）与 complex（默认）均通过（GNU 12.1.0，无 OpenMP）。
- DMRG：real 链 D=100→300、complex 链 D=100→300，E0/site 均为 −0.11183
  （5 位小数一致；complex 测量虚部恒 0）。
- 兼容性：旧 complex 二进制（2026-06-18 构建）vs 新 complex 二进制
  单参考键 12 文件字节级一致。
- verify（`verify_232524/`）：v1 scs 单参考键 PASS（360 行，scaled re-diff
  8.8e-3，虚部 0）；v2 measure1 PASS（632 行，scaled re-diff 1.23e4 按
  generic 口径仅报告；虚部 0）；`SMOKE VERIFY ALL OK`。
- 口径说明：generic 模式（v2）比较两条独立收敛链，re-diff 受 D=300 收敛
  精度限制（sz 涨落差异 ~1.2e-4），仅报告不作判据；scs 模式（v1）保留
  完整数值硬判据（见 `docs/COMPLEX_RUN_PROTOCOL.md`）。
- 结论：复值计算链路在实值态上重现实值结果（验收 1/2 通过）；不验证手性
  态选择、不构成 d+id 证据。

## D 调度与收敛标准

- `--D=100,500,800,1000,1200,1500`（逗号分隔、无空格、从小到大；每段固定
  Dmin=Dmax=D 跑满 Sweeps；`MaxLanczIter` 从 3 线性爬升到 JSON 值）。
- 收敛标准（params.json）：`Sweeps`、`CutOff`、`LanczErr`、`MaxLanczIter`；
  noise 数组直接来自 JSON。

## 测量版本与口径

- 测量二进制：`measureSC_fix_ref_ext`（**单参考键基线语义**，
  `--refbond=i,j`，与归档数据逐字节兼容——B1.1 已核实的算符定义与弦因子
  约定未变）；`measure1` 及其余测量二进制未改动。
- 复值协议（任务 C）仅引入构建配置变更（`RealCode OFF` = 默认复值路径）
  与验收工具 `tools/verify_complex_reproduces_real.py`（单参考键 manifest
  门控 + generic 双向完备性 + E0 比对），不改变测量输出口径。
- 任务 A/B 实现已按用户决定还原（提交 `f268fac` 历史中保留，可恢复）。

## 源点 / refbond 选择

- 归档口径：N200 源点 49/50，refbonds `49,50`、`50,51`、`50,54`；
  N400 源点 97/98（消费侧审计 P1-G.23 记录几何/参考键/量子态三个层次）。
- 配对矩阵协议：refbond 子集与 `--reverse-channels` 组合在集群执行方案中
  按对称性逐案明确。

## 数据生成参数

- 几何：`HoneycombXC`（当前 params.json：Lx=12, Ly=4 → N=96，NumHole=2）；
  t=0.3, t2=0.01, J=0.1, J2=0.01, phi=0（honeycomb 忽略 phi）。
- 消费侧归档层级：`Data_N{200,400}/{t_J}/{doping}/{bond_dim}_measure_{...}`；
  Data_Energy 保存 slurm 优化日志（逐格点 TruncErr 在 DMRG stdout 中，
  随日志归档）。

## 纪律

- 协议 tag 创建后，参与复跑的数据生成代码若发生任何改动，必须升级协议
  版本号并记录 diff；
- 归档数据目录中的文件一律不动；
- 集群计算任务先出执行方案（D 调度、N 序列、doping 清单、refbond 子集）
  经批准后提交。
