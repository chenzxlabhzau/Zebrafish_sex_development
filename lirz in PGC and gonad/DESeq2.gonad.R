###########gonad 
gonad=read.csv("/home/lirz/shui/gonad.count",sep = "", header = T)
row.names(gonad) <- gonad[,1]
gonad = gonad[,-c(1,20:27)]
gonad = gonad[-c(32521:32525),]
gonad = unique(gonad)
gonad = gonad[rowSums(gonad)!=0,]
gonad <- as.matrix(gonad)


library(DESeq2)
condition<-factor(c("big","big","big","small","small","small","big","big","big","small","small","small","big","big","big","small","small","small"))
name<-c(colnames(gonad))
time=factor(c("20","20","20","20","20","20","25","25","25","25","25","25","30","30","30","30","30","30"))

dds <- DESeqDataSetFromMatrix(gonad,DataFrame(condition,time), design= ~time + condition + time:condition)
dds <- DESeq(dds)

resdata <- as.data.frame(counts(dds, normalized=TRUE))
res <- results(dds)
res = res[order(res$padj),]
res = as.data.frame(res)
summary(res)

#Normalized 
nor_gonad <- as.data.frame(counts(dds,normalized = T))
write.csv(nor_gonad,file = "/home/lirz/shui/gonad.nor_total.csv")

#diff genes
write.csv(res,file = "/home/lirz/shui/gonad.diffshui.csv")
table(res$padj<0.05)
diff_gene_deseq2 <-subset(res, padj < 0.05 & abs(log2FoldChange) > 1)

#heatmap
select<-order(rowMeans(counts(dds,normalized=TRUE)),decreasing = TRUE)
nt<-normTransform(dds)
log2.norm.counts<-assay(nt)[select,]
df<-as.data.frame(colData(dds))


library(pheatmap)
pheatmap(log2.norm.counts,cluster_rows = TRUE,show_rownames = FALSE,cluster_cols = TRUE,annotation_col = df)
a = pheatmap(cor(resdata,method = "spearman"),display_numbers = TRUE,number_color = "black")

#PCA
library(pheatmap)
library(Rtsne)
library(ggplot2)
library(ggfortify)
library(mvtnorm)
df = rbind(as.data.frame(nor_gonad),group=c(colnames(nor_gonad))) %>%
  t()
pca = prcomp(t(nor_gonad)) 
autoplot(pca,data = df,colour = 'group',size=5 )+theme_bw()
