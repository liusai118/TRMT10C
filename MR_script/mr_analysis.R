library(data.table)
Sys.getenv("R_ENVIRON_USER")
ieugwasr::user()
library(TwoSampleMR)
library(dplyr) 
library(ieugwasr)
library(Matrix)
library(data.table)
library(dplyr)
library(TwoSampleMR)
library(VariantAnnotation)

library(MendelianRandomization) 
library(TwoSampleMR) 

library(data.table)
library(dplyr)
library(VariantAnnotation)

library(ieugwasr)
###############
exp <- fread("mvAge.summary.EUR.txt")
library(CMplot)
a_all <- as.data.frame(exp)
a_all <- a_all[,c(1,2,3,9)]
a_all$CHR <- paste0("chr",a_all$CHR)
library(ggsci)
CMplot(a_all,
       plot.type="m",
       LOG10=TRUE,
       cex = 0.4,
       col= pal_nejm()(8),
       #highlight = SNPs,
       highlight.col = NULL,
       highlight.cex = 0.1,
       #highlight.pch = c(15:17),
       #highlight.text = genes,
       #highlight.text.col = "black",
       ##threshold=c(1e-5,1e-3),
       #threshold.lty = 2,
       #threshold.col = "black",
       #threshold.lwd = 3,
       amplify = FALSE,
       dpi = 720,
       file.output = T,
       verbose = F,height = 6,width = 12)
library(dplyr)
library(ChIPseeker)
#library(annotatr)
library(dplyr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
a_all2 <- a_all %>%
  dplyr::filter(Pvalue < 1e-6)
regions <- GRanges(seqnames = a_all2$CHR, ranges = IRanges(start = a_all2$BP, end =a_all2$BP))
annotated_regions <- annotatePeak(regions, TxDb = txdb)
annotated_df <- as.data.frame(annotated_regions)
plotAnnoPie(annotated_regions)
plotAnnoBar(annotated_regions)


library(ggplot2)
library(clusterProfiler)
library(yulab.utils)
ids <- bitr(annotated_df$geneId,'ENTREZID','SYMBOL','org.Hs.eg.db')

ego_ALL <- enrichGO(gene = annotated_df$geneId,
                    OrgDb=org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    minGSSize = 1,
                    pvalueCutoff = 1,
                    qvalueCutoff = 1,
                    readable = TRUE) 

annotated_df <- annotated_df %>%
  mutate(GeneSymbol = mapIds(
    org.Hs.eg.db,
    keys = geneId,
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"  
  ))
table(annotated_df$GeneSymbol)
gene <- as.data.frame(table(annotated_df$GeneSymbol))
write.csv(ego_ALL,"go1.csv")
target <- c("MHC protein complex assembly",
            "regulation of T cell activation",
            "positive regulation of lymphocyte activation",
            "dopamine uptake",
            "energy homeostasis",
            "learning or memory",
            "regulation of synapse organization",
            "excitatory postsynaptic potential")


p5 <- ggplot(ego_ALL@result[ego_ALL@result$Description %in% target,], 
             aes(x = Count, y = Description, color = -log(pvalue), size = Count)) +
  geom_point() +
  scale_color_gradientn(
    name = "-log10(pvalue)",
    values = seq(0, 1, 0.5),
    colors = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25")
  ) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
  scale_x_continuous(limits = c(1, 30)) +
  geom_point(shape = 21, color = "black", stroke = 1, alpha = 1) +  # 添加黑色边框
  scale_size_continuous(name = "Count", range = c(4, 8), breaks = c(2, 10, 20, 40)) +  # 调整点的大小范围和刻度，添加图例标题
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # 添加边框
    axis.line = element_line(colour = "black"),  # 添加坐标轴线
    panel.background = element_rect(fill = "white"),  # 添加面板背景
    axis.text.x = element_text(angle = 45, hjust = 1),  # 旋转 x 轴标签
    axis.title.x = element_blank(),  # 移除 x 轴标题
    axis.title.y = element_blank(),  # 移除 y 轴标题
    axis.text = element_text(family = "Arial", size = 12),  # 设置坐标轴文本字体为 Arial 12 号
    panel.grid.major = element_blank(),  # 去除主要网格线
    panel.grid.minor = element_blank()  # 去除次要网格线
  )

print(p5)
ggsave(plot = p5,filename = "go1.svg",width = 5.7,height =3.2 )
##################################################
