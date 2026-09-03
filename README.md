# real_data — mixFunMap 手稿五个分析的简化复现代码

本目录收录 mixFunMap 手稿中五个实证/模拟分析的**大幅简化复现脚本**。原始分析代码
（`manuscript/code/`、`real_data_analysis/`、`ril/scripts/`、`Fig/` 等）保持原样不动；
这里的五个脚本是独立、紧凑的复现版本，去掉了原流水线中的检查点、日志、审计与
集群调度层，只保留必要的处理、分析与可视化代码，模型设定与数值实现与原稿完全一致。

每个脚本分两部分：

1. **分析段（`full` 模式）**：从 `inputs/` 中的分析就绪输入重新运行拟合与全基因组
   扫描，复现 `results/` 中的结果。
2. **绘图段（默认模式）**：不重新拟合任何模型，直接读取 `results/` 中冻结的结果，
   即时重绘手稿全部图版（Figure 1–5 及补充图 S1–S5 的 PDF + TIFF）。

## 目录结构

```
real_data/
├── 01_simulation.R            # V15 形式化 logistic 模拟（Figure 1, S1–S5）
├── 02_wheat_height.R          # FIP1 小麦株高 GWAS（Figure 2）
├── 03_staph_mic0.R            # 金黄色葡萄球菌 MIC=0 生长曲线 GWAS（Figure 3）
├── 04_rice_ril.R              # 水稻 RIL 双环境（对照/低水）分析（Figure 4）
├── 05_wheat_canopy_cover.R    # FIP1 小麦冠层覆盖度 LOP 分析（Figure 5）
├── inputs/                    # 分析就绪输入（rds / csv）
│   ├── wheat_height/          #   FPWW012 株高表型+基因型（含 FM 用 0/1/2 硬判定）
│   ├── staph_mic0/            #   growdata_mic0.csv / snpdata.csv / FASTSTRUCTURE.csv
│   ├── rice_ril/              #   RiceCGM/RDP1 双环境联合输入
│   └── wheat_canopy_cover/    #   FPWW012 冠层覆盖度 LOP 输入
├── results/                   # 冻结的分析结果（绘图段读取；full 模式可覆盖重建）
│   ├── simulation/summary/    #   模拟汇总（17,200 数据集的指标汇总）
│   ├── wheat_height/          #   三方法扫描结果、阈值、Top4 可解释性 rds
│   ├── staph_mic0/            #   三方法扫描结果、Top4 可解释性 rds
│   ├── rice_ril/              #   control / lowwater / gxe 及各基准方法子目录
│   ├── wheat_cc/              #   LOP 选阶、混合模型扫描、比较方法、区域注释
│   └── annotation/            #   Table S 显著位点注释（v2 布局，供绘图标注基因名）
└── figures/                   # 运行绘图段后输出的全部图版（PDF + TIFF）
```

## 用法

所有脚本默认只绘图（秒级完成）：

```bash
Rscript 01_simulation.R
Rscript 02_wheat_height.R
Rscript 03_staph_mic0.R
Rscript 04_rice_ril.R
Rscript 05_wheat_canopy_cover.R
```

重新运行完整分析（耗时较长，需安装 mixFunMap 包与 GMMAT）：

```bash
Rscript 01_simulation.R full [cores] [n_rep_null] [n_rep_alt]   # 小样本可做试点
Rscript 02_wheat_height.R full [cores] [output_dir]
Rscript 03_staph_mic0.R  full [cores] [output_dir]
Rscript 04_rice_ril.R    full [cores] [output_dir]
Rscript 05_wheat_canopy_cover.R full [cores] [output_dir]
```

注意：01 的正式设计为 28 个零假设单元 × 300 重复 + 44 个备择单元 × 200 重复
= 17,200 个数据集，full 模式请传小重复数做试点（如 `full 4 2 1`）。

## 依赖

- **绘图模式**：R ≥ 4.4，`ggplot2`、`ggrepel`、`readxl`（01 另需 `gridExtra`）。
- **full 模式**：额外需要 mixFunMap 包（仓库根目录 `./mixFunMap`，安装或置于
  `.libPaths`）与 `GMMAT` 包（minP 基准）。

## 模型与阈值概要

| 分析 | 均值曲线 | mixFunMap 设定 | 阈值 |
|---|---|---|---|
| 模拟 (V15) | 三参数 logistic (1, 0.65, 7) | 3 df P3D Wald，固定 τ² 零模型 | 0.05/m |
| 小麦株高 | 三参数 logistic | Q=PC1–5 + VanRaden K，3 df P3D Wald | Bonf. 0.05/M_eff=1.53e-5（M_eff=3,270）；建议性 1/M=5.38e-5（M=18,583） |
| 金黄色葡萄球菌 | 标准 logistic | Q=PC1–3 + K，3 df P3D Wald | Bonf. 0.05/M；建议性 1/M |
| 水稻 RIL | 三参数 logistic（单环境）/ 联合 G×E | 单环境 3 df Wald；联合 `fit_mixfunmap_joint` 3 df Wald | Bonf. 0.05/M_eff（M_eff=6,904）；建议性 1/M（M=33,697） |
| 小麦冠层覆盖度 | LOP（BIC 选阶 → 4 阶） | Q=PC1–5 + K，5 df P3D Wald | 同小麦株高 |

每个脚本头部注释给出了与该分析完全一致的模型规格、QC 步骤与数据来源细节。

## 数据来源

- **小麦株高 / 冠层覆盖度**：FIP1（GABI-WHEAT）FPWW012 群体，90K SNP 芯片。
- **金黄色葡萄球菌**：99 个菌株，MIC=0 条件下 14 个时间点生长曲线
  （growdata/snpdata/FASTSTRUCTURE 三个原始 csv）。
- **水稻**：RiceCGM/RDP1 RIL 群体（349 个家系，21 天，33,697 个标记，
  对照与低水两个环境）。
- **模拟**：手稿 V15 配置（P3D-Wald、相对 QTL 方向 (1,1,−1)、K 中性因果标记、
  n=100、目标家系大小 5、Q/噪声比 0.10、K/噪声比 0.50、Henderson REML）。
