# HyBsNet

`HyBsNet` 是 HyBs Figure 3 网络分析的独立 R 包。包函数只处理数据和
R 对象，不依赖当前工作目录，也不会在未明确指定时覆盖已有结果。

当前保留的科学参数为：

- Figure 3 DEG 集合：BH 校正后 `P < 0.05`；
- 网络边：Jaccard similarity `> 0.10`；
- 社群：Louvain；
- 核心/外围：每个疾病网络复合中心性中位数；
- 默认布局随机种子：`42`。

## 安装

从 `HyBs_Fig3` 项目根目录运行：

```bash
R CMD INSTALL packages/HyBsNet
```

也可以安装 `R CMD build` 生成的 `HyBsNet_0.1.0.9000.tar.gz`。

项目中的 `.r-lib` 保存的是 Windows 构建，不应在 macOS 上放到
`.libPaths()` 的首位。安装和验证请使用当前 macOS R library。

### 从 GitHub 安装

首次使用时安装一次，之后直接 `library()`：

```r
install.packages("pak")
pak::pak("maxbond233/HyBsNet")

library(HyBsNet)
```

更新到 GitHub 最新版本时重新运行：

```r
pak::pak("maxbond233/HyBsNet")
```

## 分析与绘图分离

```r
library(HyBsNet)

deg_raw <- read_hybs_csv("MAST_deg_summary.csv")
deg <- prepare_deg_sets(deg_raw, condition = "Obesity", padj_cutoff = 0.05)

edge_raw <- read_hybs_csv("deg_similarity_all_clusters_jaccard_obesity.csv")
network <- build_similarity_network(
  edge_raw,
  deg_data = deg,
  min_similarity = 0.10,
  condition = "Obesity"
)

layout <- calculate_network_layout(network, method = "fr", seed = 42)
plot <- plot_similarity_network(
  network,
  layout = layout,
  labels = select_network_labels(network, mode = "top_degree", n = 10),
  style = theme_fig3_current()
)
```

只改颜色时，复用同一个 `layout`：

```r
new_style <- hybs_network_style(
  region_colors = c(BS = "#0072B2", SH = "#D55E00", TH = "#009E73")
)

plot_similarity_network(network, layout = layout, style = new_style)
```

## 冻结新版网络坐标

手动调整或 community-first 布局应先转换成带拓扑签名的正式布局。之后
绘制任意基因或通路时都复用这个对象：

```r
frozen_layout <- as_network_layout(
  redesigned_nodes,
  network,
  x_col = "x_redesign",
  y_col = "y_redesign",
  layout_id = "log2fc_1p5_community_first_v1",
  layout_method = "community_first_frozen"
)

save_network_layout(frozen_layout, "diabetes_frozen_layout.csv")
```

`load_network_layout()` 默认严格要求节点和拓扑一致。如果阈值变化导致
节点替换，必须显式调用 `reconcile_network_layout()`；共同节点坐标保持不变，
新增节点按已定位邻居加权锚定。

## 叠加基因或通路结果

```r
gene_overlay <- prepare_gene_overlay(
  mast_table,
  gene = "HLA-E",
  condition = "Diabetes"
)

plot_network_overlay(
  network,
  gene_overlay,
  layout = frozen_layout,
  size_by = "overlay"
)

pathway_overlay <- prepare_pathway_overlay(
  gsea_table,
  pathway = "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  condition = "Diabetes"
)

plot_network_overlay(
  network,
  pathway_overlay,
  layout = frozen_layout,
  size_by = "deg_count"
)
```

默认 gene overlay 的填色是 `avg_log2FC`，不是绝对表达量；大小默认是
`pct.1`。GSEA 的 `NES` 表示疾病差异排序中的富集方向，也不应直接命名为
单细胞 pathway activity。UCell、AUCell 或 GSVA 等分数可以通过
`prepare_node_overlay()` 输入，并用准确的 `score_type` 标记。

如果输入的是平均表达量或 pseudobulk 表达量，使用语义更明确的
`prepare_expression_overlay()`：

```r
expression_overlay <- prepare_expression_overlay(
  expression_table,
  gene = "SST",
  condition = "Diabetes",
  value_col = "avg_expression",
  size_col = "pct_expressing"
)
```

## 提取并绘制子网络

```r
local <- extract_subnetwork(
  network,
  mode = "ego",
  center_node = "SH_IN_10_ESR1",
  order = 1,
  edge_mode = "induced"
)

local_frozen <- calculate_subnetwork_layout(
  local,
  mode = "frozen",
  parent_layout = frozen_layout
)

plot_subnetwork(local, local_frozen, overlay = gene_overlay)
plot_subnetwork_context(network, frozen_layout, local, overlay = gene_overlay)
```

`extract_subnetwork()` 还支持 `nodes`、`pattern`、`community` 和 `group`
查询。局部图可以复用全局冻结坐标，也可以通过 `mode = "compact"` 生成
独立的紧凑布局；紧凑布局不会写回或改变全局坐标。

完整、可在 RStudio 中 Source 的入口位于：

```text
inst/scripts/run_fixed_network_views.R
```

也可以直接修改安装包中的 YAML 模板副本，然后读取：

```r
style_cfg <- read_network_config("figure3_style.yml")
style <- style_from_config(style_cfg)
```

## 当前 Figure 3 一键重建

安装包后运行：

```bash
Rscript inst/scripts/rebuild_figure3_networks.R \
  /path/to/HyBs_Fig3 \
  /path/to/new_output_directory
```

脚本从冻结的 `results/current/` 表开始，生成 A、B、G 的 PDF、PNG、SVG
及 source data，同时生成 F 的中心性/DEG burden 数据。不会重跑 MAST。

## 验证

```r
testthat::test_package("HyBsNet")
```

真实数据回归：

```bash
Rscript tools/validate_current_figure3.R /path/to/HyBs_Fig3
Rscript tools/validate_network_views_v1.R /path/to/HyBs_Fig3
```
