---
name: dmrgtj-calculation
description: DMRGtJModel 项目计算工作流指南：在 honeycomb/平方晶格 t-J 模型上准备参数（params.json）、绝热演化续跑、提交集群 DMRG 任务（submit_dmrg_zk.sh）与运行关联测量（measure_corr.sh / measure_super_corr.sh）。当用户要跑新参数的 DMRG 基态、换参数重算、继续/重启已有计算、选择测量源点或参考键、修改 SLURM 脚本或 params.json 时使用——只要涉及本仓库的计算或测量流程就应先读本 skill，即使用户没提"DMRG"或"skill"。数据提取、作图、拟合与相态判断不在本 skill 范围（由 DMRG-data-processing 项目的 skills 负责）。
---

# DMRGtJModel 计算工作流

本 skill 覆盖 DMRGtJModel 仓库的**计算侧**完整流程：准备参数 → 绝热演化得到正确基态 → 提交集群 DMRG 任务 → 运行关联测量 → 把数据交给后处理项目。代码架构与构建细节见仓库根目录 CLAUDE.md，本 skill 只讲"怎么算对、怎么跑"。

## 1. 项目分工（先想清楚任务属于哪边）

三个相关项目各司其职，不要在错误的项目里做重复工作：

| 任务 | 项目/位置 |
|---|---|
| 跑 DMRG、绝热演化、测量（本 skill） | 本仓库 + 集群 `/share/home/zhangkun/DMRGtJModel` |
| 数据提取、作图、拟合、相态判断 | `C:\Users\ZhangKun\Desktop\DMRG-data-processing`（有自己的 CLAUDE.md 和 3 个 skill：`dmrg-density-plot`、`dmrg-spin-corr-plot`、`dmrg-supercorr-plot`） |
| UltraDMRG 库算法/接口开发 | `C:\Users\ZhangKun\Desktop\UltraDMRG`（有 `ultradmrg-algorithm-dev` skill） |

## 2. 环境与同步

- 本机 Windows 只有源码，**编译与计算都在 SUSTech 集群**（Intel oneAPI + MKL + MPI，build 目录 `/share/home/zhangkun/DMRGtJModel/build`）。
- 同步方式：本地改代码/脚本后 `git commit` + `git push`，集群上 `git pull`，然后 `sbatch`。
- 每个计算放在独立的 run 目录（`rundir=$(pwd)`）：内含自己的 `params.json`、`mps/`、`mpo/`、`.temp/` 与 slurm 日志。SLURM 脚本按当前目录工作，所以必须在 run 目录内提交。
- 参数文件格式见 `docs/parameters.md`；`params.json` 在仓库根是当前活跃示例，各 run 目录用自己的副本。

## 3. 关键机制：MPS 续跑与 .temp

- 若 run 目录已有 `mps/` 且张量数量与系统大小一致，`dmrg` 会**自动加载旧 MPS 作为初态继续优化**，不会重新初始化；否则从 Néel 态 + 等间距空穴开始。
- `.temp/` 存有限 VMPS 算法的环境张量，每次跑 DMRG 都会生成，占磁盘较多（宽系统可达 TB 级）。它是**可重建的临时状态，不是续跑或测量所必需的状态文件**：同一物理参数下删除 `.temp`，不影响加大 bond dimension 的续跑，也不影响物理量测量（程序会自动重建，代价是重建时间）。
- **删除规则**：
  - 同一物理参数（仅加大 D、重新提交、测量）：`.temp` **可以删**，不需要为续跑或测量而保留；
  - 换物理参数绝热演化：每次切换到下一组参数前，**必须删除**上一组参数计算产生的 `.temp`；
  - 删除前确认目标文件确实是 `.temp`（核对路径与生成方式），避免误删 `mps/`、测量结果或其他必要数据；不要把规则泛化为删除整个临时目录或其他计算数据。

## 4. 绝热演化（换参数算基态的正确方式）

**原理**：直接在新参数点上从 Néel 初态出发，DMRG 容易陷入局部极小/亚稳态，得不到真正的基态。正确做法是从一个"容易算对"的相邻参数点出发，把已收敛的基态 MPS 作为新参数点的初态，逐步演化到目标参数。每步参数变化要小，且每步都确认收敛后再走下一步；参数变化较大时不要假设一步就能到位，应插入更细的中间参数点。演化路径不唯一，最终是否得到真正的基态以能量、截断误差、纠缠熵稳定为准，并可与其他初态（如从头 Néel 态）或已有计算的结果交叉验证——取能量更低、判据更一致者。**不要断言一步演化必然足够或必然得到基态。**

**案例（honeycomb t-J，目标 t=0.3, J=0.1 纯 t-J 基态）**：

1. 第一步：`params.json` 设为 `t=0.3, J=0.1, t2=0.01, J2=0.01`（带 NNN 微扰的"好算"邻域），`sbatch submit_dmrg_zk.sh` 跑出收敛的 MPS。
2. 第二步：**删除 run 目录的 `.temp/`**。
3. 第三步：同一 run 目录改 `params.json` 为 `t=0.3, J=0.1, t2=0, J2=0`，再次 `sbatch submit_dmrg_zk.sh`。程序加载第一步的 MPS 作为初态继续优化。
4. 检查收敛并交叉验证：能量、截断误差、纠缠熵随 sweep 稳定才算收敛。本例 t2/J2 变化小，一步通常足够；若未收敛或疑似亚稳态，插入中间参数点（如 t2=J2=0.005）逐步演化，并可另跑一次从头 Néel 初态的计算对比能量（取更低者）。

## 5. 提交 DMRG 计算（submit_dmrg_zk.sh）

- 脚本里标"每次要修改的参数"处改 `D_VALUES_DMRG`（bond dimension 序列，**升序、逗号分隔、无空格**，如 `"100,500,800,1000,1200,1500"`）；程序逐值跑 sweep 收敛。
- 提交：`sbatch submit_dmrg_zk.sh`（8 MPI × 7 threads，SLURM 配置按需改）。
- 判断收敛：看 slurm 日志的 `Site (i, j) E0 = ... TruncErr = ... S = ...` 行，能量/截断误差/熵随 sweep 稳定即收敛；不收敛则加大 D 或增加 sweep。
- 注意事项：honeycomb 几何忽略 `phi`；非零 `phi` 需复数构建（`RealCode=OFF`）；`NumHole` 影响初态空穴放置（子晶格等间距）。

## 6. 测量

### 6.1 measure_corr.sh（单点量 + 两点关联）

- `sbatch measure_corr.sh`，改脚本中"sites 扫描设置"处的 `sites` 列表。
- `sites` 是**两点关联函数的源点**（measure1 的第 2 个参数）；单点量 `nf.json`/`sz.json` 与 sites 无关，任何取值都不影响。
- 源点选择：处理的是 OBC 边界条件的 t-J 模型，边界效应大，源点应选在**离边界约 1/4 系统长度处**。已有数据约定：N=200 → sites 49/50，N=400 → sites 97/98（与后处理项目解析的 left_bound 规则匹配）。
- 输出移入 `${bond}_measure_phys_qty/measure_${site}/`：`nf.json`、`sz.json`、`nf{site}nf`、`sz{site}sz`、`sp{site}sm`、`sm{site}sp`、`single_particle_corr1..4`。

### 6.2 measure_super_corr.sh（四点关联 → 超导关联）

- `sbatch measure_super_corr.sh`，改脚本中 `refbonds` 列表（格式 `"i,j"`）。
- refbond 同样**选在离边界 1/4 处**（N=200 → `49,50`/`50,51`/`50,54`；N=400 → `97,98`/`98,99`/`98,102`）。
- 输出 12 个文件 `scs{X}{ch}_fix_refbond_{i}_{j}.json`（X=A/B/C 键方向 × ch=a/b/c/d 通道），移入 `${bond}_measure_super_corr/super_corr_{i}_{j}/`。
- 超导关联 `<Δ†Δ> = 0.5(a - b - c + d)` 的合成、SU(2) 对称性检验与作图在**后处理项目**完成（`dmrg-supercorr-plot` skill），本仓库只负责产出这 12 个文件。

## 7. 数据交接给后处理项目

- 后处理项目按目录结构读取：`Data_N{200|400}/{t_J}/{doping}/{bond_dim}_measure_phys_qty/measure_{version}/` 与 `.../{bond_dim}_measure_super_corr/super_corr_{i}_{j}/`。集群 run 目录的命名应与此对应（doping 目录层按空穴率组织）。
- 完整目录约定、哪些组合实际有数据、数据质量事实见后处理项目 CLAUDE.md——**交给它之前先读一遍**。

## 8. 涉及 UltraDMRG 接口时

需要新增/修改算法接口、排查库行为（如 `TwoSiteFiniteVMPS`、`FiniteDMRG`、MPOGenerator）时，先读 `C:\Users\ZhangKun\Desktop\UltraDMRG` 的 CLAUDE.md（分层架构）与 `ultradmrg-algorithm-dev` skill（算法开发模式：单进程 API、MPI master/slave 三文件模式、测试注册规范）。**查阅范围控制**：普通参数配置、提交计算、测量、数据交接等流程问题，优先查阅本仓库 CLAUDE.md、docs/ 与 SLURM 脚本，不查库源码；只有确实涉及内部接口或算法实现、且现有文档与代码证据不足以判断时，才扩大到 UltraDMRG 源码检查。

## 9. 常见坑

- **`.temp` 两条规则**：同一物理参数下可以删（不影响加大 D 续跑和测量，程序自动重建）；换物理参数绝热演化前必须删。删前确认目标是 `.temp`，勿误删 `mps/` 或测量结果。
- **换系统大小时**：`mps/` 张量数与新 N 不一致会自动重新初始化——旧 MPS 会被覆盖（`Dump` 覆盖写），需要保留旧结果就先拷贝 `mps/`。
- **`vmps*` 系列**必须先跑 `mpogen` 生成 `./mpo`；`dmrg` 不需要（内存中构建 MPO）。
- **根目录与 run 目录的 params.json 不是一回事**：SLURM 脚本用的是 run 目录里那份。
- 蜂窝晶格 `N = 2*Lx*Ly`；`Lx`/`Ly` 语义与平方晶格不同（2 格点/原胞）。
