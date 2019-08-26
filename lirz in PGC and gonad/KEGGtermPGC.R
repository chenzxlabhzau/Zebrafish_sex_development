###diff Buc vs dnd
##select padj<0.05,foldchang>1(up)
##KEGG enrich
library(clusterProfiler)
gene_list_Buc_dnd_up <- mapIds(org.Dr.eg.db, keys = row.names(resSig_Buc_dnd_up),
                    column = "ENTREZID", keytype = "ENSEMBL" )

kegg_Buc_dnd_up <- enrichKEGG(gene_list_Buc_dnd_up, organism="dre",
                 keyType = "ncbi-geneid",
                 pvalueCutoff=0.05, pAdjustMethod="BH",
                 qvalueCutoff=0.1)

kegg_Buc_dnd_up_result=kegg_Buc_dnd_up@result

ggplot(kegg_Buc_dnd_up_result,aes(p.adjust,Description)) +
  geom_point(aes(size=Count,color=-1*log10(p.adjust)))+
  scale_color_gradient(low="green",high ="red")+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
  
   
##select padj<0.05,foldchang<1(down)
##KEGG enrich
library(clusterProfiler)
gene_list_Buc_dnd_down <- mapIds(org.Dr.eg.db, keys = row.names(resSig_Buc_dnd_down),
                               column = "ENTREZID", keytype = "ENSEMBL" )

kegg_Buc_dnd_down <- enrichKEGG(gene_list_Buc_dnd_down, organism="dre",
                              keyType = "ncbi-geneid",
                              pvalueCutoff=0.05, pAdjustMethod="BH",
                              qvalueCutoff=0.1)

kegg_Buc_dnd_down_result=kegg_Buc_dnd_down@result

ggplot(kegg_Buc_dnd_down_result,aes(p.adjust,Description)) +
  geom_point(aes(size=Count,color=-1*log10(p.adjust)))+
  scale_color_gradient(low="green",high ="red")+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))


###diff for Buc vs wt
##select padj<0.05,foldchang>1(up)
##KEGG enrich
library(clusterProfiler)
gene_list_Buc_wt_up <- mapIds(org.Dr.eg.db, keys = row.names(resSig_Buc_dnd_up),
                              column = "ENTREZID", keytype = "ENSEMBL" )

kegg_Buc_wt_up <- enrichKEGG(gene_list_Buc_wt_up, organism="dre",
                             keyType = "ncbi-geneid",
                             pvalueCutoff=0.05, pAdjustMethod="BH",
                             qvalueCutoff=0.1)

kegg_Buc_wt_up_result=kegg_Buc_wt_up@result

ggplot(kegg_Buc_wt_up_result,aes(p.adjust,Description)) +
  geom_point(aes(size=Count,color=-1*log10(p.adjust)))+
  scale_color_gradient(low="green",high ="red")+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
  
##select padj<0.05,foldchang<1(down)
##KEGG enrich
library(clusterProfiler)
gene_list_Buc_wt_down <- mapIds(org.Dr.eg.db, keys = row.names(resSig_Buc_dnd_down),
                              column = "ENTREZID", keytype = "ENSEMBL" )

kegg_Buc_wt_down <- enrichKEGG(gene_list_Buc_wt_down, organism="dre",
                             keyType = "ncbi-geneid",
                             pvalueCutoff=0.05, pAdjustMethod="BH",
                             qvalueCutoff=0.1)

kegg_Buc_wt_down_result=kegg_Buc_wt_down@result

ggplot(kegg_Buc_wt_down_result,aes(p.adjust,Description)) +
  geom_point(aes(size=Count,color=-1*log10(p.adjust)))+
  scale_color_gradient(low="green",high ="red")+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))



###diff for dnd vs wt
##select padj<0.05,foldchang>1(up)
##KEGG enrich
library(clusterProfiler)
gene_list_dnd_wt_up <- mapIds(org.Dr.eg.db, keys = row.names(resSig_dnd_wt_up),
                              column = "ENTREZID", keytype = "ENSEMBL" )

kegg_dnd_wt_up <- enrichKEGG(gene_list_dnd_wt_up, organism="dre",
                             keyType = "ncbi-geneid",
                             pvalueCutoff=0.05, pAdjustMethod="BH",
                             qvalueCutoff=0.1)

kegg_dnd_wt_up_result=kegg_dnd_wt_up@result

ggplot(kegg_dnd_wt_up_result,aes(p.adjust,Description)) +
  geom_point(aes(size=Count,color=-1*log10(p.adjust)))+
  scale_color_gradient(low="green",high ="red")+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
  
  
##select padj<0.05,foldchang<1(down)
##KEGG enrich
library(clusterProfiler)
gene_list_dnd_wt_down <- mapIds(org.Dr.eg.db, keys = row.names(resSig_dnd_wt_down),
                              column = "ENTREZID", keytype = "ENSEMBL" )

kegg_dnd_wt_down <- enrichKEGG(gene_list_dnd_wt_down, organism="dre",
                             keyType = "ncbi-geneid",
                             pvalueCutoff=0.05, pAdjustMethod="BH",
                             qvalueCutoff=0.1)

kegg_dnd_wt_down_result=kegg_dnd_wt_down@result

ggplot(kegg_dnd_wt_down_result,aes(p.adjust,Description)) +
  geom_point(aes(size=Count,color=-1*log10(p.adjust)))+
  scale_color_gradient(low="green",high ="red")+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
