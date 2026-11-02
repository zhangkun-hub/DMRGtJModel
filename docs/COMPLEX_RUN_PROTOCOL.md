# 复值基态与复数配对关联运行协议（任务 C）

## 构建配置

`RealCode` 默认 **OFF**，即默认构建就是复值路径（`-DUSE_COMPLEX=1`，
`CMakeLists.txt:9-14`）。归档数据是 `-DRealCode=ON`（`-DUSE_REAL=1`）实值版。

复值构建（与现有集群构建命令一致，**不加** `-DRealCode=ON`）：

```bash
cmake .. -DCMAKE_PREFIX_PATH=<path-to-qlten-and-qlmps>
make -j4
```

`HIGHER_SYMMETRY` 保持默认 ON（SYM_LEVEL=1，U(1)×U(1) 守恒 N 与 Sz）。
测量与优化代码对 `TenElemT` 全模板化；复值下测量 JSON 值写 `[re, im]`
（实值下写标量），消费侧解析器已兼容两者。

## 验收标准（严格分离三层能力）

1. **复值计算能力**：复值程序在相同物理态、相同参数及一致测量约定下重现实值
   计算结果；
2. **实值态保持**：复值 MPS 表示的实值态，其测量虚部在数值误差范围内为零；
3. **手性态选择能力** 与 **d+id 物理验证能力** 是两个独立问题：复值 DMRG
   能运行**不代表**体系已形成 d+id 手性态。

验收 1/2 用 `tools/verify_complex_reproduces_real.py` 自动比对（**完整性门控
先行，数值比对在后；任一完整性错误均 FAIL**）：

```bash
# 配对测量目录（scs 文件）：独立 expected manifest 门控
python tools/verify_complex_reproduces_real.py \
    --complex <complex_scs_dir> --real <real_scs_dir> \
    --complex-log <complex_dmrg_log> --real-log <real_dmrg_log> \
    --scs-mode --geometry HoneycombXC --Ly 4 --Lx 12 \
    --refbonds 49,50;50,51;50,54 [--reverse-channels] \
    -o complex_reproduces_real_report.json

# 其他测量目录（如 measure1 输出）：双向完备性门控（无独立 manifest）
python tools/verify_complex_reproduces_real.py \
    --complex <complex_measure1_dir> --real <real_measure1_dir> \
    --complex-log <complex_dmrg_log> --real-log <real_dmrg_log> \
    -o complex_reproduces_real_report.json
```

完整性门控（缺/多文件、缺/多行、重复行、JSON 解析错误、非有限数值、
scs 模式的方向标签不一致，任一不得 PASS）；dmrg 日志**必填**（`E0/site`
比对是验收一部分）。数值判据：`|Re(c)−r|`、`|Im(c)|` 与 E0 差均须 ≤
`atol + rtol·|r|`（默认 atol=rtol=1e-8，实际容差写入报告）。
scs 模式的 expected manifest 由参考键列表 + honeycomb 几何枚举 + 目标键
过滤规则（j0 > i1）+ 方向分类 + 通道组合独立构造（`--refbonds ALL` 可取
全部 NN 键）。

## 运行协议

- 用与归档数据完全相同的物理参数（几何、Lx/Ly、t/t2/J/J2、NumHole、D 调度、
  Sweeps/CutOff/LanczErr/MaxLanczIter、noise 数组）跑复值 dmrg + 测量；
- 初始态沿用实值 Néel+等间距空穴（`dmrg.cpp:134-182`，未改动）——这是
  **实值可表示态**，收敛后测量虚部应 ≈ 0（验收 2）；
- 数据落新目录，**不覆盖归档实值数据**；
- d+id 判定基于**方向分辨的配对关联**：用 A/B/C 方向 × a-d（及 e-h）通道
  构造 d_{x²−y²} 与 d_{xy} 两个配对基底，分析其复数关联与相对相位；
  **不得**将单个空间关联函数的相位平台直接等同于 d+id 判据。

## 手性态选择（开放项，待物理决定，不阻塞代码交付）

复值 DMRG 从实初态出发通常收敛到实值态（差一个整体相位）。要寻找手性态需
显式打破时间反演对称的初态/扰动，候选方案（需单独立项或协议决定）：

1. **noise 注入**（已支持）：`noised_vmps` / JSON `noise` 数组在优化中注入
   随机扰动，复值算术下可能自发进入手性扇区；
2. **相位扭曲初态**：dmrg.cpp 增加初态模式（环流/相位涡旋直积态）——需小改
   `dmrg.cpp`，本任务未实施；
3. **phi 插值**：φ≠0 通量构型优化后取 φ→0（复值专属路径，`tJmodel.h:63-67`
   在实值下强制报错）。

无论采用哪种方案，均需按验收 1/2 先确认复值计算能力，再讨论手性态证据。

## 数据格式与兼容性

- 复值测量输出 `[[索引], [re, im]]`；实值输出标量。消费侧
  `LoadSCData.m/LoadLocalSCData.m` 已按 `numel==1` 判定兼容；
- 新数据用新目录命名（集群执行方案中明确），归档数据目录一律不动。

## 已知风险

- 集群 qlmps 与本地 UltraDMRG 版本不一致；复值路径在集群的实际编译状态需
  冒烟验证（见集群执行方案）。若暴露 real-only 编译错误，按实情最小修复并
  记录（预期 ≤ 1 文件）。
