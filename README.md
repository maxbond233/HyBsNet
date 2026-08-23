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
```
