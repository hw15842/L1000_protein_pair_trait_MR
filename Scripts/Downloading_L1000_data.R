
devtools::install_github("cmap/cmapR")

library("cmapR")

df_col <- read_gctx_meta("/Users/hw15842/Downloads/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx", dim="col")


df_row <- read_gctx_meta("/Users/hw15842/Downloads/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx", dim="row")
