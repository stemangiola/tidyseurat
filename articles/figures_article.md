# Code for producing the figures in the article

[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)

\
`# Article workflow`\
\
[`library`](https://rdrr.io/r/base/library.html)`(`[`tidyverse`](https://tidyverse.tidyverse.org)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`Seurat`](https://satijalab.org/seurat)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(``SingleR``)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`plotly`](https://plotly-r.com)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`tidyHeatmap`](https://www.r-project.org)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`ggalluvial`](http://corybrunson.github.io/ggalluvial/)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`ggplot2`](https://ggplot2.tidyverse.org)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`tidyseurat`](https://github.com/stemangiola/tidyseurat)`)`\
[`options`](https://rdrr.io/r/base/options.html)`(``future.globals.maxSize ``=`` ``50068`` ``*`` ``1024``^``2``)`\
\
`# Use colourblind-friendly colours`\
`friendly_cols`` ``<-`` ``dittoSeq``::``dittoColors``(``)`\
\
`# Set theme`\
`custom_theme`` ``<-`\
`  `[`list`](https://rdrr.io/r/base/list.html)`(`\
`    `[`scale_fill_manual`](https://ggplot2.tidyverse.org/reference/scale_manual.html)`(``values ``=`` ``friendly_cols``)``,`\
`    `[`scale_color_manual`](https://ggplot2.tidyverse.org/reference/scale_manual.html)`(``values ``=`` ``friendly_cols``)``,`\
`    `[`theme_bw`](https://ggplot2.tidyverse.org/reference/ggtheme.html)`(``)`` ``+`\
`      `[`theme`](https://ggplot2.tidyverse.org/reference/theme.html)`(`\
`        panel.border ``=`` `[`element_blank`](https://ggplot2.tidyverse.org/reference/element.html)`(``)``,`\
`        axis.line ``=`` `[`element_line`](https://ggplot2.tidyverse.org/reference/element.html)`(``)``,`\
`        panel.grid.major ``=`` `[`element_line`](https://ggplot2.tidyverse.org/reference/element.html)`(``size ``=`` ``0.2``)``,`\
`        panel.grid.minor ``=`` `[`element_line`](https://ggplot2.tidyverse.org/reference/element.html)`(``size ``=`` ``0.1``)``,`\
`        text ``=`` `[`element_text`](https://ggplot2.tidyverse.org/reference/element.html)`(``size ``=`` ``9``)``,`\
`        legend.position ``=`` ``"bottom"``,`\
`        strip.background ``=`` `[`element_blank`](https://ggplot2.tidyverse.org/reference/element.html)`(``)``,`\
`        axis.title.x ``=`` `[`element_text`](https://ggplot2.tidyverse.org/reference/element.html)`(``margin ``=`` `[`margin`](https://ggplot2.tidyverse.org/reference/element.html)`(``t ``=`` ``10``, r ``=`` ``10``, b ``=`` ``10``, l ``=`` ``10``)``)``,`\
`        axis.title.y ``=`` `[`element_text`](https://ggplot2.tidyverse.org/reference/element.html)`(``margin ``=`` `[`margin`](https://ggplot2.tidyverse.org/reference/element.html)`(``t ``=`` ``10``, r ``=`` ``10``, b ``=`` ``10``, l ``=`` ``10``)``)``,`\
`        axis.text.x ``=`` `[`element_text`](https://ggplot2.tidyverse.org/reference/element.html)`(``angle ``=`` ``30``, hjust ``=`` ``1``, vjust ``=`` ``1``)`\
`      ``)`\
`  ``)`\
\
`PBMC_clean_scaled_UMAP_cluster_cell_type`` ``<-`` `[`readRDS`](https://rdrr.io/r/base/readRDS.html)`(``"dev/PBMC_clean_scaled_UMAP_cluster_cell_type.rds"``)`

\
`p1`` ``=`` `\
`  ``PBMC_clean_scaled_UMAP_cluster_cell_type`` `[`%>%`](../reference/pipe.md)\
`  `[`pivot_longer`](../reference/pivot_longer.md)`(`\
`    `[`c`](https://rdrr.io/r/base/c.html)`(``mito.fraction``, ``S.Score``, ``G2M.Score``)``, `\
`    names_to``=``"property"``, `\
`    values_to``=``"Value"`\
`  ``)`` `[`%>%`](../reference/pipe.md)\
`  `[`mutate`](../reference/mutate.md)`(``property ``=``  `[`factor`](https://rdrr.io/r/base/factor.html)`(``property``, levels ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"mito.fraction"``, ``"G2M.Score"``, ``"S.Score"``)``)``)`` `[`%>%`](../reference/pipe.md)\
`  `[`ggplot`](../reference/ggplot.md)`(`[`aes`](https://ggplot2.tidyverse.org/reference/aes.html)`(``sample``, ``Value``)``)`` ``+`` `\
`  `[`geom_boxplot`](https://ggplot2.tidyverse.org/reference/geom_boxplot.html)`(``outlier.size ``=`` ``0.5`` ``)`` ``+`` `\
`  `[`facet_wrap`](https://ggplot2.tidyverse.org/reference/facet_wrap.html)`(``~``property``, scales ``=`` ``"free_y"`` ``)`` ``+`\
`  ``custom_theme`` ``+`\
`  `[`theme`](https://ggplot2.tidyverse.org/reference/theme.html)`(``aspect.ratio``=``1``)`

\
`p2`` ``=`` `\
`  ``PBMC_clean_scaled_UMAP_cluster_cell_type`` `[`%>%`](../reference/pipe.md)\
`  `[`sample_n`](../reference/sample_n.md)`(``20000``)`` `[`%>%`](../reference/pipe.md)\
`  `[`ggplot`](../reference/ggplot.md)`(`[`aes`](https://ggplot2.tidyverse.org/reference/aes.html)`(``UMAP_1``, ``UMAP_2``, color``=``seurat_clusters``)``)`` ``+`\
`  `[`geom_point`](https://ggplot2.tidyverse.org/reference/geom_point.html)`(``size``=``0.05``, alpha``=``0.2``)`` ``+`\
`  ``custom_theme`` ``+`\
`  `[`theme`](https://ggplot2.tidyverse.org/reference/theme.html)`(``aspect.ratio``=``1``)`\
\
`PBMC_clean_scaled_UMAP_cluster_cell_type`` `[`%>%`](../reference/pipe.md)\
`  `[`sample_n`](../reference/sample_n.md)`(``20000``)`` `[`%>%`](../reference/pipe.md)\
`  `[`plot_ly`](../reference/plotly.md)`(`\
`    x ``=`` ``~``` `UMAP_1` ```,`\
`    y ``=`` ``~``` `UMAP_2` ```,`\
`    z ``=`` ``~``` `UMAP_3` ```,`\
`    color ``=`` ``~``seurat_clusters``,`\
`    colors ``=`` ``friendly_cols``[``1``:``24``]``,sizes ``=`` ``50``, size ``=`` ``1`\
`  ``)`\
\
`markers`` ``=`` `[`readRDS`](https://rdrr.io/r/base/readRDS.html)`(``"dev/PBMC_marker_df.rds"``)`

\
`p3`` ``=`` `\
`  ``PBMC_clean_scaled_UMAP_cluster_cell_type`` `[`%>%`](../reference/pipe.md)\
`  `[`arrange`](../reference/arrange.md)`(``first.labels``)`` `[`%>%`](../reference/pipe.md)\
`  `[`mutate`](../reference/mutate.md)`(``seurat_clusters ``=`` ``fct_inorder``(``seurat_clusters``)``)`` `[`%>%`](../reference/pipe.md)\
`  `[`join_features`](../reference/join_features.md)`(``features``=`[`c`](https://rdrr.io/r/base/c.html)`(``"CD3D"``, ``"HLA-DRB1"``)``)`` `[`%>%`](../reference/pipe.md)\
`  `[`ggplot`](../reference/ggplot.md)`(`[`aes`](https://ggplot2.tidyverse.org/reference/aes.html)`(``y``=``seurat_clusters`` , x``=``.abundance_SCT``, fill``=``first.labels``)``)`` ``+`\
`  ``geom_density_ridges``(``bandwidth ``=`` ``0.2``)`` ``+`\
`  `[`facet_wrap`](https://ggplot2.tidyverse.org/reference/facet_wrap.html)`(``~`` ``.feature``, nrow ``=`` ``2``)`` ``+`\
`  `[`coord_flip`](https://ggplot2.tidyverse.org/reference/coord_flip.html)`(``)`` ``+`\
`  ``custom_theme`

\
`# Plot heatmap`\
`p4`` ``=`` `\
`  ``PBMC_clean_scaled_UMAP_cluster_cell_type`` `[`%>%`](../reference/pipe.md)\
`  `[`sample_n`](../reference/sample_n.md)`(``2000``)`` `[`%>%`](../reference/pipe.md)\
`  `[`DoHeatmap`](https://satijalab.org/seurat/reference/DoHeatmap.html)`(`\
`    features ``=`` ``markers``$``gene``,`\
`    group.colors ``=`` ``friendly_cols`\
`  ``)`

\
`p5`` ``=`` `\
`  ``PBMC_clean_scaled_UMAP_cluster_cell_type`` `[`%>%`](../reference/pipe.md)\
`  `[`sample_n`](../reference/sample_n.md)`(``1000``)`` `[`%>%`](../reference/pipe.md)\
`  `[`join_features`](../reference/join_features.md)`(``features``=``markers``$``gene``)`` `[`%>%`](../reference/pipe.md)\
`  `[`mutate`](../reference/mutate.md)`(``seurat_clusters ``=`` `[`as.integer`](https://rdrr.io/r/base/integer.html)`(``seurat_clusters``)``)`` `[`%>%`](../reference/pipe.md)\
`  `[`filter`](../reference/filter.md)`(``seurat_clusters``<``10``)`` `[`%>%`](../reference/pipe.md)\
`  `[`group_by`](../reference/group_by.md)`(``seurat_clusters``)`` `[`%>%`](../reference/pipe.md)\
`  `\
`  ``# Plot heatmap`\
`  `[`heatmap`](https://rdrr.io/r/stats/heatmap.html)`(`\
`    .row ``=`` ``.feature``,`\
`    .column ``=`` ``.cell``, `\
`    .value ``=`` ``.abundance_SCT``, `\
`    palette_grouping ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(`[`rep`](https://rdrr.io/r/base/rep.html)`(``"black"``,``9``)``)``, `\
`    palette_value ``=`` ``circlize``::``colorRamp2``(`[`c`](https://rdrr.io/r/base/c.html)`(``-``1.5``, ``0``, ``1.5``)``, `[`c`](https://rdrr.io/r/base/c.html)`(``"purple"``, ``"black"``, ``"yellow"``)``)``,`\
`    `\
`    ``# ComplexHeatmap parameters`\
`    row_gap ``=`` `[`unit`](https://rdrr.io/r/grid/unit.html)`(``0.1``, ``"mm"``)``, column_gap ``=`` `[`unit`](https://rdrr.io/r/grid/unit.html)`(``0.1``, ``"mm"``)`\
`  ``)`` `[`%>%`](../reference/pipe.md)\
`    `\
`  ``# Add annotation`\
`  ``add_tile``(``sample``, palette ``=`` ``friendly_cols``[``1``:``7``]``)`` `[`%>%`](../reference/pipe.md)\
`  ``add_point``(``PC_1``)`` `

\
`p6`` ``=`` `\
`  ``PBMC_clean_scaled_UMAP_cluster_cell_type`` `[`%>%`](../reference/pipe.md)\
`  `[`unite`](../reference/unite.md)`(``"cluster_cell_type"``, `[`c`](https://rdrr.io/r/base/c.html)`(``first.labels``, ``seurat_clusters``)``, remove``=``FALSE``)`` `[`%>%`](../reference/pipe.md)\
`  `[`pivot_longer`](../reference/pivot_longer.md)`(`\
`    `[`c`](https://rdrr.io/r/base/c.html)`(``seurat_clusters``, ``first.labels_single``)``,`\
`    names_to ``=`` ``"classification"``, values_to ``=`` ``"value"`\
`  ``)`` `[`%>%`](../reference/pipe.md)\
`  `\
`  `[`ggplot`](../reference/ggplot.md)`(`[`aes`](https://ggplot2.tidyverse.org/reference/aes.html)`(``x ``=`` ``classification``, stratum ``=`` ``value``, alluvium ``=`` ``cell``,`\
`           fill ``=`` ``first.labels``, label ``=`` ``value``)``)`` ``+`\
`  `[`scale_x_discrete`](https://ggplot2.tidyverse.org/reference/scale_discrete.html)`(``expand ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``1``, ``1``)``)`` ``+`\
`  ``geom_flow``(``)`` ``+`\
`  ``geom_stratum``(``alpha ``=`` ``.5``)`` ``+`\
`  ``# geom_text(stat = "stratum", size = 3) +`\
`  ``geom_text_repel``(``stat ``=`` ``"stratum"``, size ``=`` ``3``,`\
`             nudge_x      ``=`` ``0.05``,`\
`             direction    ``=`` ``"y"``,`\
`             angle        ``=`` ``0``,`\
`             vjust        ``=`` ``0``,`\
`             segment.size ``=`` ``0.2`\
`         ``)`` ``+`\
`  `[`scale_fill_manual`](https://ggplot2.tidyverse.org/reference/scale_manual.html)`(``values ``=`` ``friendly_cols``)`` ``+`\
`  ``#guides(fill = FALSE) +`\
`  `[`coord_flip`](https://ggplot2.tidyverse.org/reference/coord_flip.html)`(``)`` ``+`\
`  `[`theme_bw`](https://ggplot2.tidyverse.org/reference/ggtheme.html)`(``)`` ``+`\
`  `[`theme`](https://ggplot2.tidyverse.org/reference/theme.html)`(`\
`    panel.border ``=`` `[`element_blank`](https://ggplot2.tidyverse.org/reference/element.html)`(``)``,`\
`    axis.line ``=`` `[`element_line`](https://ggplot2.tidyverse.org/reference/element.html)`(``)``,`\
`    panel.grid.major ``=`` `[`element_line`](https://ggplot2.tidyverse.org/reference/element.html)`(``size ``=`` ``0.2``)``,`\
`    panel.grid.minor ``=`` `[`element_line`](https://ggplot2.tidyverse.org/reference/element.html)`(``size ``=`` ``0.1``)``,`\
`    text ``=`` `[`element_text`](https://ggplot2.tidyverse.org/reference/element.html)`(``size ``=`` ``9``)``,`\
`    legend.position ``=`` ``"bottom"``,`\
`    strip.background ``=`` `[`element_blank`](https://ggplot2.tidyverse.org/reference/element.html)`(``)``,`\
`    axis.title.x ``=`` `[`element_text`](https://ggplot2.tidyverse.org/reference/element.html)`(``margin ``=`` `[`margin`](https://ggplot2.tidyverse.org/reference/element.html)`(``t ``=`` ``10``, r ``=`` ``10``, b ``=`` ``10``, l ``=`` ``10``)``)``,`\
`    axis.title.y ``=`` `[`element_text`](https://ggplot2.tidyverse.org/reference/element.html)`(``margin ``=`` `[`margin`](https://ggplot2.tidyverse.org/reference/element.html)`(``t ``=`` ``10``, r ``=`` ``10``, b ``=`` ``10``, l ``=`` ``10``)``)``,`\
`    axis.text.x ``=`` `[`element_text`](https://ggplot2.tidyverse.org/reference/element.html)`(``angle ``=`` ``30``, hjust ``=`` ``1``, vjust ``=`` ``1``)`\
`  ``)`
