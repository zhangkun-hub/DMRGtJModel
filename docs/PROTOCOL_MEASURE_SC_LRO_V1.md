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

## 集群环境（待回填）

- [ ] qlten 安装路径与版本：`__________`
- [ ] qlmps 安装路径与版本：`__________`
- [ ] module 列表（提交脚本中的 `module load ...` 序列）：`__________`
- [ ] 编译器版本（`icpc --version`）：`__________`

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
