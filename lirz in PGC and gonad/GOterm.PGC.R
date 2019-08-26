#GO enrich
library(DOSE)
library(org.Dr.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(AnnotationHub)


#Transform ID
resSigup$id <- as.character(resSigup$id)
transID = bitr(resSigup$id,
               fromType="ENTREZID",
               toType=c("ENSEMBL", "ENTREZID"),
               OrgDb="org.Dr.eg.db"
) 

###diff for Buc vs dnd
library(DESeq2)
condition_Buc_dnd<-factor(c("Buc","Buc","dndMo","dndMo"))
Buc_dnd=total[,1:4]
dds_Buc_dnd <- DESeqDataSetFromMatrix(Buc_dnd, DataFrame(condition_Buc_dnd), design= ~ condition_Buc_dnd)
dds_Buc_dnd <- DESeq(dds_Buc_dnd)

resdata <- as.data.frame(counts(dds_Buc_dnd, normalized=TRUE))
res_Buc_dnd <- results(dds_Buc_dnd)
res_Buc_dnd = res_Buc_dnd[order(res_Buc_dnd$padj),]
res_Buc_dnd = as.data.frame(res_Buc_dnd)


###diff Buc vs dnd
##select padj<0.05,foldchang>1(up)
resSig_Buc_dnd_up=res_Buc_dnd[res_Buc_dnd$padj<0.05 & res_Buc_dnd$log2FoldChange>1,]
resSig_Buc_dnd_up<-na.omit(resSig_Buc_dnd_up)
resSig_Buc_dnd_up$id=row.names(resSig_Buc_dnd_up) 


#Biological process
Buc_dnd_up_BP <- enrichGO(resSig_Buc_dnd_up$id,
               "org.Dr.eg.db",
               keyType="ENSEMBL",
               ont="BP"
)

barplot(Buc_dnd_up_BP, showCategory=20,title = "Biological process")


#Molecular Function
Buc_dnd_up_MF <- enrichGO(resSig_Buc_dnd_up$id,
                          "org.Dr.eg.db",
                          keyType="ENSEMBL",
                          ont="MF"
)

barplot(Buc_dnd_up_MF, showCategory=20,title = "Molecular Function")


#Cellular Component
Buc_dnd_up_CC <- enrichGO(resSig_Buc_dnd_up$id,
                          "org.Dr.eg.db",
                          keyType="ENSEMBL",
                          ont="CC"
)

barplot(Buc_dnd_up_CC, showCategory=20,title = "Cellular Component")



##select padj<0.05,foldchang<1(down)
resSig_Buc_dnd_down=res_Buc_dnd[res_Buc_dnd$padj<0.05 & res_Buc_dnd$log2FoldChange<c(-1),]
resSig_Buc_dnd_down<-na.omit(resSig_Buc_dnd_down)
resSig_Buc_dnd_down$id=row.names(resSig_Buc_dnd_down) 


#Biological process
Buc_dnd_down_BP <- enrichGO(resSig_Buc_dnd_down$id,
               "org.Dr.eg.db",
               keyType="ENSEMBL",
               ont="BP"
)

barplot(Buc_dnd_down_BP,title = "Biological process")


#Molecular Function
Buc_dnd_down_MF <- enrichGO(resSig_Buc_dnd_down$id,
                          "org.Dr.eg.db",
                          keyType="ENSEMBL",
                          ont="MF"
)

barplot(Buc_dnd_down_MF, showCategory=20,title = "Molecular Function")


#Cellular Component
Buc_dnd_down_CC <- enrichGO(resSig_Buc_dnd_down$id,
                          "org.Dr.eg.db",
                          keyType="ENSEMBL",
                          ont="CC"
)

barplot(Buc_dnd_down_CC, showCategory=20,title = "Cellular Component")

###diff for Buc vs wt
library(DESeq2)
condition_buc_wt<-factor(c("Buc","Buc","Wt","Wt"))
Buc_wt=total[,c(1,2,5,6)]
dds_Buc_wt <- DESeqDataSetFromMatrix(Buc_wt, DataFrame(condition_buc_wt), design= ~ condition_buc_wt)
dds_Buc_wt <- DESeq(dds_Buc_wt)

resdata <- as.data.frame(counts(dds_Buc_wt, normalized=TRUE))
res_Buc_wt <- results(dds_Buc_wt)
res_Buc_wt = res_Buc_dnd[order(res_Buc_wt$padj),]
res_Buc_wt = as.data.frame(res_Buc_wt)

##select padj<0.05,foldchang>1(up)
resSig_Buc_wt_up=res_Buc_wt[res_Buc_wt$padj<0.05 & res_Buc_wt$log2FoldChange>1,]
resSig_Buc_wt_up<-na.omit(resSig_Buc_wt_up)
resSig_Buc_wt_up$id=row.names(resSig_Buc_wt_up)

#Biological process
Buc_wt_up_BP <- enrichGO(resSig_Buc_wt_up$id,
                          "org.Dr.eg.db",
                          keyType="ENSEMBL",
                          ont="BP"
)

barplot(Buc_wt_up_BP, showCategory=20,title = "Biological process")


#Molecular Function
Buc_wt_up_MF <- enrichGO(resSig_Buc_wt_up$id,
                            "org.Dr.eg.db",
                            keyType="ENSEMBL",
                            ont="MF"
)

barplot(Buc_wt_up_MF, showCategory=20,title = "Molecular Function")


#Cellular Component
Buc_wt_up_CC <- enrichGO(resSig_Buc_wt_up$id,
                            "org.Dr.eg.db",
                            keyType="ENSEMBL",
                            ont="CC"
)

barplot(Buc_wt_up_CC, showCategory=20,title = "Cellular Component")

##select padj<0.05,foldchang<1(down)
resSig_Buc_wt_down=res_Buc_wt[res_Buc_wt$padj<0.05 & res_Buc_wt$log2FoldChange<c(-1),]
resSig_Buc_wt_down<-na.omit(resSig_Buc_wt_down)
resSig_Buc_wt_down$id=row.names(resSig_Buc_wt_down) 

#Biological process
Buc_wt_down_BP <- enrichGO(resSig_Buc_wt_down$id,
                            "org.Dr.eg.db",
                            keyType="ENSEMBL",
                            ont="BP"
)

barplot(Buc_wt_down_BP,title = "Biological process")


#Molecular Function
Buc_wt_down_MF <- enrichGO(resSig_Buc_wt_down$id,
                           "org.Dr.eg.db",
                           keyType="ENSEMBL",
                           ont="MF"
)

barplot(Buc_wt_down_MF, showCategory=20,title = "Molecular Function")


#Cellular Component
Buc_wt_down_CC <- enrichGO(resSig_Buc_wt_down$id,
                           "org.Dr.eg.db",
                           keyType="ENSEMBL",
                           ont="CC"
)

barplot(Buc_wt_down_CC, showCategory=20,title = "Cellular Component")

###diff for dnd vs wt
library(DESeq2)
condition_dnd_wt<-factor(c("dnd","dnd","Wt","Wt"))
dnd_wt=total[,c(3,4,5,6)]
dds_dnd_wt <- DESeqDataSetFromMatrix(dnd_wt, DataFrame(condition_dnd_wt), design= ~ condition_dnd_wt)
dds_dnd_wt <- DESeq(dds_dnd_wt)

resdata <- as.data.frame(counts(dds_dnd_wt, normalized=TRUE))
res_dnd_wt <- results(dds_dnd_wt)
res_dnd_wt = res_dnd_wt[order(res_dnd_wt$padj),]
res_dnd_wt = as.data.frame(res_dnd_wt)

##select padj<0.05,foldchang>1(up)
resSig_dnd_wt_up=res_dnd_wt[res_dnd_wt$padj<0.05 & res_dnd_wt$log2FoldChange>1,]
resSig_dnd_wt_up<-na.omit(resSig_dnd_wt_up)
resSig_dnd_wt_up$id=row.names(resSig_dnd_wt_up)

#Biological process
dnd_wt_up_BP <- enrichGO(resSig_dnd_wt_up$id,
                         "org.Dr.eg.db",
                         keyType="ENSEMBL",
                         ont="BP"
)

barplot(dnd_wt_up_BP, showCategory=20,title = "Biological process")

#Molecular Function
dnd_wt_up_MF <- enrichGO(resSig_dnd_wt_up$id,
                           "org.Dr.eg.db",
                           keyType="ENSEMBL",
                           ont="MF"
)

barplot(dnd_wt_up_MF, showCategory=20,title = "Molecular Function")


#Cellular Component
dnd_wt_up_CC <- enrichGO(resSig_dnd_wt_up$id,
                           "org.Dr.eg.db",
                           keyType="ENSEMBL",
                           ont="CC"
)

barplot(dnd_wt_up_CC, showCategory=20,title = "Cellular Component")

##select padj<0.05,foldchang<1(down)
resSig_dnd_wt_down=res_dnd_wt[res_dnd_wt$padj<0.05 & res_dnd_wt$log2FoldChange<c(-1),]
resSig_dnd_wt_down<-na.omit(resSig_dnd_wt_down)
resSig_dnd_wt_down$id=row.names(resSig_dnd_wt_down) 

#Biological process
dnd_wt_down_BP <- enrichGO(resSig_dnd_wt_down$id,
                           "org.Dr.eg.db",
                           keyType="ENSEMBL",
                           ont="BP"
)

barplot(dnd_wt_down_BP,title = "Biological process")

#Molecular Function
dnd_wt_down_MF <- enrichGO(resSig_dnd_wt_down$id,
                           "org.Dr.eg.db",
                           keyType="ENSEMBL",
                           ont="MF"
)

barplot(dnd_wt_down_MF, showCategory=20,title = "Molecular Function")


#Cellular Component
dnd_wt_down_CC <- enrichGO(resSig_dnd_wt_down$id,
                           "org.Dr.eg.db",
                           keyType="ENSEMBL",
                           ont="CC"
)

barplot(dnd_wt_down_CC, showCategory=20,title = "Cellular Component")
© 2019 GitHub, Inc.
