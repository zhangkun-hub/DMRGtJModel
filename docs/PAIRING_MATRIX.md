# 跨参考键配对矩阵（任务 B）

## 定义

P_{αβ}(i,j) = ⟨Δ†_{i,α} Δ_{j,β}⟩，其中 i、j 为（最近邻）键，α、β ∈ {A, B, C}
为键方向。Δ†_{i,α} 为键 i 上的配对产生算符（自旋组合由通道给出）。

## 键方向 A/B/C 的映射（已核实）

方向分类规则即 `HoneycombXC/YC::ClassifyNearestNeighbor`
（`src/honeycomb_lattice.h:179-241`），输入键按站点序 (i<j)：

- **HoneycombXC**（y 周期、x 开放）：
  - 同列（x0==x1，竖直键）：`y_min = ((y0+Ly−y1)%Ly==1 ? y0 : y1)`；
    `(y_min+x0)%2==0 → A`，否则 B。
  - 同行（y0==y1，水平键）：C。
- **HoneycombYC**（y 周期、x 开放，zigzag 方向与 XC 互换）：
  - 同行（y0==y1，水平键）：`(min(x0,x1)+y0)%2==0 → A`，否则 B。
  - 同列（x0==x1，竖直键）：C。

注意：A/B/C 的几何含义在 XC 与 YC 之间不同（C 在 XC 是水平键、在 YC 是竖直键）。
Python 镜像参考实现：`tools/honeycomb_geometry.py`（已与 C++ 逐条对比验证）。

## 通道 a/b/c/d（及 e/f/g/h）与算符串

参考对 (r0,r1)、目标对 (t0,t1)，局域算符约定：0=↑、1=↓、2=空穴；
`bupc=b†↑`、`bupa=b↑`、`bdnc=b†↓`、`bdna=b↓`；弦算符 f = diag(−1,−1,+1) 只插
参考对内与目标对内（`finite_mps_measu_disk.h:325-344`，**B1.1 已核实的语义，
本扩展未改变**）。

通道 a-d（`src/measureSC_fix_ref_extension.cpp` main，产生@参考对、湮灭@目标对）：

| 通道 | phys_ops = {@r0, @r1, @t0, @t1} |
|---|---|
| a | {bupc, bdnc, bupa, bdna} |
| b | {bdnc, bupc, bupa, bdna} |
| c | {bupc, bdnc, bdna, bupa} |
| d | {bdnc, bupc, bdna, bupa} |

通道 e-h（`--reverse-channels` 时输出，湮灭@参考对、产生@目标对，为 a-d 的
角色互换共轭）：

| 通道 | phys_ops |
|---|---|
| e | {bupa, bdna, bupc, bdnc} |
| f | {bdna, bupa, bupc, bdnc} |
| g | {bupa, bdna, bdnc, bupc} |
| h | {bdna, bupa, bdnc, bupc} |

单重态/三重态组合仍在二进制外由消费侧进行（口径不变）。

## 方向语义

- α = **参考键**方向，β = **目标键**方向；
- 文件名标签 X（`scs{X}...`）= **目标键**方向（目标键枚举按方向分组）；
- 参考键方向 α **不在文件名中**，须由消费侧按参考键索引用同一分类规则
  （`tools/honeycomb_geometry.py` 的 `classify_nn_bond`）获得；
- 通道 a-d/e-h 的自旋组合映射见上表，共轭配对为 a↔e、b↔f、c↔g、d↔h。

## 测量覆盖范围与 API 约束

`MeasureFourSiteOpGroup` 要求目标键严格位于参考键之后：**j0 > i1**（现有按
格点序阈值过滤 `t0 ≤ r1` 剔除，恰好就是该 API 约束）。据此对全部 NN 键对
（真实几何枚举，`tools/honeycomb_geometry.py` 镜像 C++ 键表）分类，得到
矩阵元素的**覆盖状态语义**：

| 状态 | 含义 |
|---|---|
| `MEASURED_DIRECT` | j0 > i1，直接测得 P(i→j) |
| `RECONSTRUCTED_HERMITIAN` | 反向元素 P(j,i)，由 R1（伴随恒等式）从直接测量值复共轭得到 |
| `NOT_MEASURED_SHARED_SITE` | 两键共享格点：两个方向均不满足 j0 > i1；R2/R3 失效，**禁止用共轭补齐** |
| `NOT_MEASURED_MPS_INTERLEAVE` | 两键**物理上不重叠（R1/R2/R3 均成立）**，但 MPS 索引区间重叠（嵌套 `i0<j0<j1<i1` 或交叉 `i0<j0<i1<j1`）：排序约束使两个方向均不满足 j0 > i1，且 e-h 通道不能绕过该约束 → **没有可用于重建的已测元素**，禁止用共轭凭空补齐 |
| `NOT_MEASURED_OTHER` | 对角元素 P(i,i)（键上配对自关联，目标=参考被过滤排除） |

**真实几何枚举统计**（矩阵元素占比，对角线 = 全部键对数 M²）：

| 几何 | N | M | MEASURED_DIRECT | RECONSTRUCTED_HERMITIAN | SHARED_SITE | INTERLEAVE | OTHER(对角) |
|---|---|---|---|---|---|---|---|
| XC Ly=4 Lx=8 | 64 | 94 | 45.87% | 45.87% | 4.16% | 3.03% | 1.06% |
| YC Ly=4 Lx=8 | 64 | 92 | 44.34% | 44.34% | 4.16% | 6.07% | 1.09% |
| XC Ly=4 Lx=12 | 96 | 142 | 47.24% | 47.24% | 2.78% | 2.04% | 0.70% |
| YC Ly=4 Lx=12 | 96 | 140 | 46.23% | 46.23% | 2.78% | 4.05% | 0.71% |

区间交错实例：XC 的 ((0,3),(1,2))（嵌套）、YC 的 ((0,4),(1,5))（交叉）。
**MPS 区间交错类在两种几何中都存在（约 2–6%），消费侧组装矩阵时必须把该
类元素显式标记为 NOT_MEASURED_MPS_INTERLEAVE，不得用共轭关系凭空补齐。**

## 共轭关系（数学推导）

定义（对任意键对 i、j）：

P_{αβ}(i,j) = ⟨Δ†_{i,α} Δ_{j,β}⟩；Q_{αβ}(i,j) = ⟨Δ_{i,α} Δ†_{j,β}⟩

**R1（伴随恒等式，恒成立，与键是否重叠无关）**：取 O = Δ†_{i,α} Δ_{j,β}，
O† = Δ†_{j,β} Δ_{i,α}，由 ⟨O⟩\* = ⟨O†⟩：

P_{βα}(j,i) = P\*_{αβ}(i,j)

**R1 是配对矩阵 Hermiticity 关系，也是反向矩阵元素的补全依据**：
下三角 = 上三角复共轭。R1 不需要 e-h 通道，也不需要键不重叠。

**R2（物理不重叠键对易）**：当两个配对算符作用在**不同的物理格点集合**上
（i ∩ j = ∅，**不要求 MPS 索引区间完全分离**）时，四格点全不同，
{c_{iσ}, c†_{jτ}} = 0，且 Δ、Δ† 都是偶费米子算符，故
[Δ_{i,α}, Δ†_{j,β}] = 0：

P_{βα}(j,i) = ⟨Δ†_{j,β} Δ_{i,α}⟩ = ⟨Δ_{i,α} Δ†_{j,β}⟩ = Q_{αβ}(i,j)

**R3（物理不重叠键的 Q-P 共轭，e-h 通道的检验对象）**：由 R1、R2 立得：

P_{αβ}(i,j) = Q\*_{αβ}(i,j)

**R3 对物理上不重叠的配对键成立，不要求 MPS 索引区间完全分离。**
对测量域内的同一有序键对 (i→j)（j0 > i1 保证 i∩j = ∅）：
**e-h 通道测得 Q，其复共轭必须等于 a-d 通道测得的 P**。算子级核对过排列
符号：⟨c_{t1↓}c_{t0↑}c_{r1↓}c_{r0↑}⟩ 到 ⟨c†_{r0↑}c†_{r1↓}c_{t0↑}c_{t1↓}⟩
需 6 次相邻交换（符号 +1），四组通道 (a↔e, b↔f, c↔g, d↔h) 置换图样相同，
无额外费米子负号、无弦因子修正。

适用边界：

1. **R1 是矩阵补全依据**（对所有键对成立，与键是否重叠、MPS 索引是否
   交错无关）；
2. **e-h 通道用于 R3 的数值一致性检验，不是 R1 的来源**；
3. **R2/R3 的适用条件只是物理格点集合不重叠**。仅当两键**共享格点**时
   R2/R3 失效（{c_{sσ}, c†_{sτ}} = δ_{στ}(1−n_{s,−σ}) ≠ 0，存在交换子
   修正）；
4. **MPS 索引交错的物理不重叠键对：R1/R2/R3 均成立**。但当前
   `MeasureFourSiteOpGroup` 的排序约束（j0 > i1）使两个方向都无法直接
   测量，e-h 通道也不能绕过该排序约束——因此数据中**没有可用于 R1 重建
   的已测元素**。此类元素标记 `NOT_MEASURED_MPS_INTERLEAVE`，不得写成
   "R2/R3 不适用"或"共轭关系不成立"，也禁止凭空补齐；
5. 数值验证：`tools/check_pairing_conjugacy.py` 对每个被测 (i→j) 检查
   `P == conj(Q)`（R3）残差，并校验文件名标签与目标键实际方向一致：

```bash
python tools/check_pairing_conjugacy.py --dir <measure_dir> \
    --geometry HoneycombXC --Ly 4 --Lx 12 --rtol 1e-6 \
    -o pairing_conjugacy_check.json
```

## 多参考键运行（配对矩阵生成）

扩展后的 `measureSC_fix_ref_ext` 支持：

```bash
mpirun -n 8 ./measureSC_fix_ref_ext params.json --refbond=49,50;50,51;50,54
mpirun -n 8 ./measureSC_fix_ref_ext params.json --all-refbonds --reverse-channels
```

- 每个参考键沿用完全相同 12 文件命名（`scs{X}{ch}_fix_refbond_{i}_{j}.json`），
  单参考键运行输出与旧版逐字节一致；
- `--reverse-channels` 追加 e-h 12 文件；
- `--all-refbonds` 将全部 NN 键作为参考键（N=96 时约 140+ 个参考键，适合小
  系统/冒烟；生产建议对称不等价参考键子集，用时在集群执行方案中明确）。
- 效率：每个参考键一次全目标键扫描，R 个参考键 ≈ R 倍单参考键耗时。

## 矩阵组装（消费侧）

- **上三角 P(i→j)**（j0 > i1）：a-d 通道文件直接读得；
- **下三角 P(j,i)**：R1 补全——上三角对应元素的复共轭（不需要额外测量，
  也**不是**由 e-h 文件补全）；
- **e-h 文件**：R3 一致性检验数据（`check_pairing_conjugacy.py` 消费），
  组装矩阵时不参与填充；
- **不可测元素**：按覆盖状态语义标注（共享格点 / MPS 索引交错 / 对角），
  消费侧矩阵中显式写 NOT_MEASURED，禁止用共轭凭空补齐；
- 每个元素的方向对 (α, β)：β 来自文件名标签，α 由参考键索引按
  `classify_nn_bond` 计算。

**样本独立性（必须遵守）**：`RECONSTRUCTED_HERMITIAN` 元素由 R1 从直接
测量值复共轭得到，**不是独立测量样本**——做统计、拟合或有效样本数计数
时，重建元素与 `MEASURED_DIRECT` 元素不得等同计数（一对 (i,j) 只贡献
一个独立样本）。

**缺失元素零填充风险**：消费侧当前管线以稀疏键对列表处理测量行
（`SC_Pair_Structure_Factor` 的 C 列表，无稠密矩阵组装），未见零填充。
若未来任何步骤把键对组装成稠密矩阵（如 dict.get(key, 0) 或 numpy 数组
默认值），**缺失元素会被静默置零**——零关联与未测量在物理上完全不同，
会制造虚假的零信号。任何矩阵化步骤必须使用带掩码的数据结构
（如 masked array / NaN 占位 + NOT_MEASURED 掩码），并在统计中排除
缺失元素。

## 与消费侧距离口径的兼容性（已数值核实）

消费侧（DMRG-data-processing `HoneycombXCGeometry.sup_corr_distance`，
`Physical_Quantity/data_records/Honeycomb_Lattice.py:95`）的 SC 距离 =
参考键与目标键**键中点**的欧氏距离，其坐标模型 = **√3 × 本仓库 XC
`SiteCoordinate`**（对全部 64 个基元格点数值验证一致，最大偏差 4e-15）。
注意其实现**不取 y 方向最小镜像**：跨 y 缝的键对距离按未缠绕坐标计算而
被高估，与本仓库 geometry 文件的 dy（最小镜像）不同——两者差异仅限于
跨缝键对。消费侧现有距离可完全由其自身坐标模型从行索引恢复，本仓库
`bond_geometry.json` 的 `dcol/dx_mid/dy_mid/dx/dy` 提供独立的最小镜像
口径，供对照使用。

## 本征分析可行性评估

数据侧可行：把配对矩阵按通道组合成 Hermitian 矩阵（P 与 conj(Q) 拼装）后
对角化，可提取主要配对关联模式及其空间结构。**结论边界**：

1. 本征矢的空间对称性**不能单独**作为晶格配对对称性或 SC-LRO 的证明，必须
   结合键方向基底、相对相位与尺寸标度分析；
2. 配对关联长程且参考键一致，**不代表已证明超导长程序**；SC-LRO 结论仍需
   质量门（`docs/MEASURE_ERROR_RECORDS.md`）与尺寸标度证据支持。
