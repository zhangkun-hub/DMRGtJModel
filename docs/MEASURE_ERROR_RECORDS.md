# 测量端误差记录（任务 A）

本文档定义 DMRGtJModel 测量端新增的四类记录：`spatial_group_spread`、
`repeat_measurement_scatter`、`truncation_error_record`、
`bond_dimension_convergence_record`。四类记录**分别存储、分别命名、不合并**。

## 语义边界（必须遵守）

四类记录都是**数据质量与收敛证据档案**，**不是独立统计噪声**：

1. `spatial_group_spread` 是同距离组内目标键的**空间位置离散度**（纯几何量），
   不是统计误差；
2. `truncation_error_record` 是优化过程的**截断误差档案**（MPS 质量历史），
   不是关联函数误差的直接估计；
3. `bond_dimension_convergence_record` 是不同 bond dimension 下同一观测量的变化，
   用于评估**数值收敛敏感性**，不是独立统计噪声；
4. `repeat_measurement_scatter` 需要真实存在的重复采样机制。当前测量程序对每个
   目标键做**单次确定性缩并**（`MeasuResElem{sites, avg}` 只存聚合值），不存在
   重复采样机制，因此**标记为 NOT_APPLICABLE，不编造任何采样方差**。

在缺少独立噪声证据的情况下，不得将关联函数尾部的微小数值直接解释为有效物理
信号；尾部信噪比及其统计显著性应标记为 **UNKNOWN / NOT VERIFIED**。

## 记录清单与格式

以下新文件由 `measureSC_fix_ref_ext`（扩展版）在测量运行时输出到当前工作目录，
与主数据文件同目录。**现有 12 个主数据文件（scs{X}{a..d}_fix_refbond_{i}_{j}.json）
的格式与口径零改动**。

### 1. `scs{X}_fix_refbond_{i}_{j}_spatial_spread.json`（spatial_group_spread）

每个键方向 X（A/B/C）一份，每个**键中点晶格位移类**一行，结构（复值构建）：

```
[
  [[1, 0], [0.501, 0]],        ← 位移类 (dx_mid=1, dy_mid=0)
  [[1, -2], [0.487, 0]],
  ...
]
```

- 第一元素 `[dx_mid, dy_mid]`：目标键键中点相对参考键键中点的**精确晶格位移类**（整数，纯整数算术，无浮点取整）：
  - `dx_mid = (x_t0 + x_t1) − (x_r0 + x_r1)`，x = site / Ly（x 方向开放，无缠绕）；
  - `dy_mid = wrap((y_t0 + y_t1) − (y_r0 + y_r1), 2·Ly)`，y = site % Ly，
    缠绕到 [−Ly, Ly)（y 方向周期，规范最小镜像类；±Ly 边界取同一类 −Ly）。
- 第二元素 `[spread, 0]`：该位移类内全部目标键**键中点位置**的 RMS 偏离
  （SiteCoordinate 物理坐标）。单键类 spread = 0；沿周期方向平移的等价键
  构成多键类，spread 反映其边界/缠绕不对称性。
- 实值构建下值与主数据同构（标量）；复值构建下写 `[spread, 0]`。
- 该分组按位移类区分 Δy 类别（同一 `dcol = t0/Ly − r0/Ly` 列位移可含多个
  dy_mid 类，本文件不混并）。**这是纯几何离散度，不涉及 MPS，不含统计含义**。

### 2. `scs{X}_fix_refbond_{i}_{j}_bond_geometry.json`（逐键几何描述）

每键一行，行序与主数据一致（按目标键首格点排序）：

```
[
  [[49, 50, 97, 98], [12, 25, 0, 21.650, 0.000, 41.570, 2.000]],
  ...
]
```

第二元素为 7 元实数组 `[dcol, dx_mid, dy_mid, dx, dy, mx, my]`：

- `dcol`：格点级列位移 `t0/Ly − r0/Ly`（与消费侧 1D 距离
  `(t0−r0)/Ly` 的整数部分一致）；
- `dx_mid, dy_mid`：同上定义的精确晶格位移类（分组键，可无损恢复）；
- `dx, dy`：物理位移 `dx = dx_mid · (px(2Ly)−px(0))/4`、
  `dy = dy_mid · (py(2)−py(0))/4`（同子格格点对推导的半列/半行步长；
  由整数类推导，**不做浮点最小镜像**，边界处无舍入歧义）；
- `mx, my`：目标键键中点的物理坐标（`SiteCoordinate`）。

约束：该几何记录要求 Lx ≥ 2 且 Ly ≥ 4（需同子格格点对推导步长），
退化几何跳过并告警。消费侧可据此用任意自定义口径重新分组。

### 3. `repeat_measurement_scatter.json`（NOT_APPLICABLE 标记）

每次运行由 rank 0 输出一次：

```json
{
  "status": "NOT_APPLICABLE",
  "reason": "no repeat sampling mechanism exists: each target bond is measured by a single deterministic contraction (MeasureFourSiteOpGroup); no stochastic sampling or independent repeated measurements are performed",
  "refbonds": [[49, 50], [50, 51]]
}
```

若未来实现真实重复采样机制（如独立随机初始态重复优化、随机测量基底等），
此文件将替换为真实波动记录；在此之前**不得编造**。

### 4. `trunc_err_record.json`（truncation_error_record，工具生成）

DMRG 优化日志（stdout/slurm）中每个两格点更新行已含 `TruncErr`（只打印到
master stdout，dmrg 二进制不落盘）。用 `tools/extract_trunc_err.py` 抽取：

```bash
python tools/extract_trunc_err.py run1.log run2.log -o trunc_err_record.json
```

输出为逐日志文件、逐 sweep、逐格点的档案：

```json
{"records": [
  {"log_file": "run1.log", "E0_per_site": -0.45,
   "D_max_sequence": [100, 500],
   "sweeps": [
     {"sweep": 1, "sites": [
        {"site_l": 0, "site_r": 1, "E0": -0.9, "TruncErr": 1.2e-05,
         "D": 100, "Iter": 3, "LanczT": 0.12, "TotT": 0.34, "S": 1.234}]}]}]}
```

**协议要求**：每次 dmrg 运行必须保留 stdout 日志，并执行抽取。
该档案是 MPS 质量历史，**不得**并入噪声底。

### 5. `bond_dim_convergence_record.json`（工具生成）

协议：同一 (几何, doping, 协议) 下，用单值 `--D=<d>` 链式跑
（`dmrg --D=d1` → 测量 → `dmrg --D=d2`（复用 mps 续跑）→ 测量 → …），
各 D 数据落不同目录，然后：

```bash
python tools/compare_bond_dim.py --ref dir_D100 --others dir_D500 dir_D1500 \
    -o bond_dim_convergence_record.json
```

输出按文件、按行（按格点索引对齐）的值差档案与每文件 `max_abs_diff`。
这是**数值收敛敏感性**记录，不是独立统计噪声。

## 兼容性

- 主数据 12 文件由原 `MeasureFourSiteOpGroup` + `DumpMeasuRes` 路径输出，
  枚举与命名未变；单参考键、默认参数的运行输出与旧版逐字节一致。
- 所有新文件均为新命名，不覆盖任何现有输出。

## 物理结论边界

- 已有收敛性与尺寸标度证据可独立使用；
- 但尾部信噪比及其统计显著性在获得独立噪声证据前标记 UNKNOWN / NOT VERIFIED。
