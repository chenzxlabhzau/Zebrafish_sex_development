#Read file
Buc_1 <- read.csv("/home/lirz/shui/Buc_1.htseqcount",sep = "", header = F)
Buc_2 <- read.csv("/home/lirz/shui/Buc_2.htseqcount",sep = "", header = F)
dndMo_1 <- read.csv("/home/lirz/shui/dndMo_1.htseqcount",sep = "", header = F)
dndMo_2 <- read.csv("/home/lirz/shui/dndMo_2.htseqcount",sep = "", header = F)
Wt_1 <- read.csv("/home/lirz/shui/Wt_1.htseqcount",sep = "", header = F)
Wt_2 <- read.csv("/home/lirz/shui/Wt_2.htseqcount",sep = "", header = F)

Buc_1=Buc_1[-c(32521:32525),]
Buc_2=Buc_2[-c(32521:32525),]
dndMo_1=dndMo_1[-c(32521:32525),] 
dndMo_2=dndMo_2[-c(32521:32525),]
Wt_1=Wt_1[-c(32521:32525),]
Wt_2=Wt_2[-c(32521:32525),]

Buc = merge(Buc_1,Buc_2,by="V1")
dndMo = merge(dndMo_1,dndMo_2,by="V1")
Wt = merge(Wt_1,Wt_2,by="V1")

Buc_dnd=merge(Buc,dndMo,by="V1") 
total=merge(Buc_dnd,Wt,by="V1")

names(total)[1:7] <- c("gene id",c(rep("Buc",2)),c(rep("dndMo",2)),c(rep("Wt",2)))
row.names(total) <- total[,1]
total = total[,-1]

total = unique(total)
total = total[rowSums(total)!=0,]
total <- as.matrix(total)

library(DESeq2)
condition<-factor(c("Buc","Buc","dndMo","dndMo","Wt","Wt"))
dds <- DESeqDataSetFromMatrix(total, DataFrame(condition), design= ~ condition)
dds <- DESeq(dds)

resdata <- as.data.frame(counts(dds, normalized=TRUE))
res <- results(dds)
res = res[order(res$padj),]
res = as.data.frame(res)
summary(res)

#Normalized 
nor_total <- as.data.frame(counts(dds,normalized = T))
write.csv(nor_total,file = "/home/lirz/shui/nor_total.csv")

#diff genes
write.csv(res,file = "/home/lirz/shui/diffshui.csv")
table(res$padj<0.05)
diff_gene_deseq2 <-subset(res, padj < 0.05 & abs(log2FoldChange) > 1)


select<-order(rowMeans(counts(dds,normalized=TRUE)),decreasing = TRUE)
nt<-normTransform(dds)
log2.norm.counts<-assay(nt)[select,]
df<-as.data.frame(colData(dds))

#heatmap
library(pheatmap)
pheatmap(log2.norm.counts,cluster_rows = TRUE,show_rownames = FALSE,cluster_cols = TRUE,annotation_col = df)
bk <- c(seq(-9,-0.1,by=0.01),seq(0,9,by=0.01))
a = pheatmap(cor(resdata,method = "spearman"),display_numbers = TRUE,number_color = "black",
             color = c(colorRampPalette(colors = c("yellow","pink"))(length(bk)/2),colorRampPalette(colors = c("pink","red"))(length(bk)/2)))

#PCA
library(pheatmap)
library(Rtsne)
library(ggplot2)
library(ggfortify)
library(mvtnorm)
pca = prcomp(df[,1:(ncol(df)-1)])
df = rbind(as.data.frame(nor_total),group=c(rep('Buc',2),rep('dndMo',2),rep('Wt',2))) %>% 
  t()
autoplot(pca,data = df,colour = 'group',size = 7)+theme_bw()


##Volcano map(Buc vs dnd)
library(ggplot2)
data=res_Buc_dnd
head(data)
data$color <- ifelse(data$padj<0.05 & abs(data$log2FoldChange)>= 1,ifelse(data$log2FoldChange > 1,'red','blue'),'gray')
color <- c(red = "red",gray = "gray",blue = "blue")

p <- ggplot(data, aes(log2FoldChange, -log10(padj), col = color)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = color) +
  labs(x="log2 (fold change)",y="-log10 (p-value)") +
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
  geom_vline(xintercept = c(-1, 1), lty=4,col="grey",lwd=0.6) +
  theme(legend.position = "none",
        panel.grid=element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
p

##Volcano map(Buc vs wt)
library(ggplot2)
data=res_Buc_wt
head(data)
data$color <- ifelse(data$padj<0.05 & abs(data$log2FoldChange)>= 1,ifelse(data$log2FoldChange > 1,'red','blue'),'gray')
color <- c(red = "red",gray = "gray",blue = "blue")

p <- ggplot(data, aes(log2FoldChange, -log10(padj), col = color)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = color) +
  labs(x="log2 (fold change)",y="-log10 (p-value)") +
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
  geom_vline(xintercept = c(-1, 1), lty=4,col="grey",lwd=0.6) +
  theme(legend.position = "none",
        panel.grid=element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
p



##Volcano map(dnd vs wt)
library(ggplot2)
data=res_dnd_wt
head(data)
data$color <- ifelse(data$padj<0.05 & abs(data$log2FoldChange)>= 1,ifelse(data$log2FoldChange > 1,'red','blue'),'gray')
color <- c(red = "red",gray = "gray",blue = "blue")

p <- ggplot(data, aes(log2FoldChange, -log10(padj), col = color)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = color) +
  labs(x="log2 (fold change)",y="-log10 (p-value)") +
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
  geom_vline(xintercept = c(-1, 1), lty=4,col="grey",lwd=0.6) +
  theme(legend.position = "none",
        panel.grid=element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
p
