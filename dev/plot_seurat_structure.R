
library( DataExplorer )

plot_str(pbmc_small, type = "r" )
plot_str(pbmc_small , type="d")

plot_str(pbmc_small %>%  join_transcripts("CD3G"), type = "r" )
plot_str(pbmc_small %>%  join_transcripts("CD3G"), type = "d" )
