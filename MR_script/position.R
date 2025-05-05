library(ChIPseeker)
#library(annotatr)
library(dplyr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
BiocManager::install("ChIPseeker")
exp2 <- a_all2 %>%
  dplyr::filter(Pvalue < 1e-6)
#exp2$CHR <- paste0("chr",exp2$CHR)
regions <- GRanges(seqnames = exp2$CHR, ranges = IRanges(start = exp2$BP, end =exp2$BP))
annotated_regions <- annotatePeak(regions, TxDb = txdb)
annotated_df <- as.data.frame(annotated_regions)
plotAnnoPie(annotated_regions)
plotAnnoBar(annotated_regions)
#################
annotated_df <- annotated_df %>%
  mutate(GeneSymbol = mapIds(
    org.Hs.eg.db,
    keys = geneId,
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"  # 如果一个 geneId 对应多个 GeneSymbol，则取第一个
  ))
a112 <- as.data.frame(table(annotated_df$GeneSymbol))

################
ids <- bitr(annotated_df$geneId,'ENTREZID','SYMBOL','org.Hs.eg.db')

library(clusterProfiler)
library(yulab.utils)
ego_ALL <- enrichGO(gene = annotated_df$geneId,#我们上面定义???
                    OrgDb=org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "ALL",#富集的GO类型
                    pAdjustMethod = "BH",#这个不用管，???般都用的BH
                    minGSSize = 1,
                    pvalueCutoff = 1,#P值可以取0.05
                    qvalueCutoff = 1,
                    readable = TRUE) 

annotated_df <- annotated_df %>%
  mutate(GeneSymbol = mapIds(
    org.Hs.eg.db,
    keys = geneId,
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"  # 如果一个 geneId 对应多个 GeneSymbol，则取第一个
  ))
gene <- as.data.frame(table(annotated_df$GeneSymbol))
write.csv(ego_ALL,"go1.csv")
target <- c("regionalization","pattern specification process","gland development","neuron projection guidance","cell fate specification","glial cell differentiation","axonogenesis","negative regulation of neuron differentiation","gliogenesis","axon development","positive regulation of neurogenesis","neuron migration","eye development","glial cell fate commitment","signal release")
p5 <- ggplot(ego_ALL@result[ego_ALL@result$Description %in% target,], 
             aes(x = Count, y = Description, color = -log(pvalue), size = Count)) +
  geom_point() +
  scale_color_gradientn(
    name = "-log10(pvalue)",
    values = seq(0, 1, 0.5),
    colors = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25")
  ) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
  scale_x_continuous(limits = c(10, 60)) +
  geom_point(shape = 21, color = "black", stroke = 1, alpha = 1) +  # 添加黑色边框
  scale_size_continuous(name = "Count", range = c(4, 8), breaks = c(10, 40, 80, 120, 160)) +  # 调整点的大小范围和刻度，添加图例标题
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