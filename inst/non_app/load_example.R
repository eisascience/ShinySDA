library(ggplot2)

envv = list()
envv$SDAres = readRDS("/Volumes/Maggie/Work/OHSU/Bimber/Expts/RIRA_manuscript/data/sda/primeseq/sda.506842.rds")
list.files("/Volumes/Maggie/Work/OHSU/Bimber/Expts/RIRA_manuscript/data/sda/primeseq/sda.506842")

envv$UMAP_CS_batch = readRDS("/Volumes/Maggie/Work/OHSU/Bimber/Expts/RIRA_manuscript/data/sda/primeseq/sda.506842/sda_results_UMAP_CellScore_QCfil_300_1_3_4_6_7_8_16_17_18_19_24_25_27_28_29_30.rds")

envv$tsne_CS_batch = readRDS("/Volumes/Maggie/Work/OHSU/Bimber/Expts/RIRA_manuscript/data/sda/primeseq/sda.506842/sda_results_tSNE_CellScore_QCfil_1000_pp1001_3_4_6_7_8_16_17_18_19_24_25_27_28_29_30.rds" )


tempDFX <- as.data.frame(envv$tsne_CS_batch$Y)
rownames(tempDFX)  <- rownames(envv$SDAres$scores)
colnames(tempDFX) <- c("tSNE1_batch", "tSNE2_batch")

tempDFX$GeneExpr <- rep(0, nrow(tempDFX))


SDAres <- envv$SDAres

GeneSet = "CD19"

GeneSet <- GeneSet[GeneSet %in% colnames(SDAres$loadings[[1]])]

GeneExpr <- SDAres$scores %*% SDAres$loadings[[1]][,as.character(GeneSet)]
colnames(SDAres$loadings[[1]])
nrow(SDAres$loadings[[1]])

if(is.null(rownames(SDAres$loadings[[1]]))){
  rownames(SDAres$loadings[[1]]) = paste0("C", 1:nrow(SDAres$loadings[[1]]))
}

TitleX = paste0("Expr of :", GeneSet )
LoadOrdVal <- round(SDAres$loadings[[1]][,as.character(GeneSet)][order(abs(SDAres$loadings[[1]][,as.character(GeneSet)]), decreasing = T)], 3)

tempDFX[rownames(GeneExpr), ]$GeneExpr <- GeneExpr[,1]


ggplot(tempDFX, aes(tSNE1_batch, tSNE2_batch,  color=cut(asinh(GeneExpr^3),
                                                         breaks = c(-Inf, -1, -.5, 0, .5, 1, Inf)))) +
  geom_point(size = 1) + theme_bw() +
  scale_color_manual("Expr", values = rev(c("red", "orange", "yellow", "lightblue", "dodgerblue", "blue")) ) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha=1))) +
  theme(legend.position = "bottom", aspect.ratio=1) + 
  coord_cartesian(xlim = NULL, ylim = NULL, expand = FALSE)  +
  labs(title = paste0("SDA-Batch-removed DGE\n", TitleX), 
       subtitle = paste("Found in comps: \n",
                        paste(names(LoadOrdVal)[1:5], collapse = ", "), 
                        "\n",
                        paste(LoadOrdVal[1:5], collapse = ", "), 
                        "\n",
                        paste(names(LoadOrdVal)[6:10], collapse = ", "), 
                        "\n",
                        paste(LoadOrdVal[6:10], collapse = ", "), 
                        "\n"), 
       caption = "Caption here") + 
  ylab("asinh(GeneExpr^3)")
