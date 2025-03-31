library(stringr)
library(openxlsx)
library(data.table)
library(reshape2)
library(dplyr)
library(tidydr)
library(tidyr)
library(tidyverse)
library(clusterProfiler)
library(pheatmap)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(fgsea)
library(corrplot)
library(colorspace)
library(survival)
library(survminer)
library(maftools)
library(vegan)
library(forcats)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(rstatix)
library(ggstatsplot)
library(survcomp)
library(customLayout)
#library(ggcor)
library(ggstance)
options(stringsAsFactors = F)

custom_theme <- function() {
  theme_survminer() %+replace%
    theme(text = element_text(family = 'Times'),panel.grid = element_blank())
}
bioForest=function(rt=null,col){
  #
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt[,2])
  hrLow  <- sprintf("%.3f",rt[,3])
  hrHigh <- sprintf("%.3f",rt[,4])
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt[,1]<0.001, "<0.001", sprintf("%.3f", rt[,1]))
  
  #
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
  
  #
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, col[2], col[1])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.5)
  axis(1)
}
getGeneFC=function(gene.exp,group,ulab=ulab,dlab=dlab){
  degs_C1_C3=mg_limma_DEG(gene.exp, 
                          group,
                          ulab=ulab,
                          dlab=dlab)
  
  ## 
  degs_C1_C3_sig<-degs_C1_C3$DEG[which(degs_C1_C3$DEG$adj.P.Val <= 1),]
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
  degs_C1_C3_sig_gene<-rownames(degs_C1_C3_sig)
  
  degs_gene_entz=bitr(degs_C1_C3_sig_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  degs_gene_entz <- dplyr::distinct(degs_gene_entz,SYMBOL,.keep_all=TRUE)
  
  gene_df <- data.frame(logFC=degs_C1_C3_sig$logFC,
                        SYMBOL = rownames(degs_C1_C3_sig))
  gene_df <- merge(gene_df,degs_gene_entz,by="SYMBOL")
  head(gene_df)
  
  geneList<-gene_df$logFC
  names(geneList)=gene_df$ENTREZID 
  head(geneList)
  
  geneList=sort(geneList,decreasing = T)
  head(geneList) ## 
  return(geneList)
}

##############
tcga.pancancer.cli=read.xlsx('TCGA_pancancer_cli_PMID_29625055.xlsx')
head(tcga.pancancer.cli)
tcga.cli=tcga.pancancer.cli[which(tcga.pancancer.cli$type=='OV'),]
head(tcga.cli)
tcga.cli=data.frame(Samples=paste0(tcga.cli$bcr_patient_barcode,'-01'),
                    Age=tcga.cli$age_at_initial_pathologic_diagnosis,
                    # Gender=tcga.cli$gender,
                    Stage=tcga.cli$clinical_stage,
                    Grade=tcga.cli$histological_grade,
                    # event=tcga.cli$new_tumor_event_type ,
                    tcga.cli[,c('OS','OS.time','DSS','DSS.time','DFI','DFI.time','PFI','PFI.time')])
rownames(tcga.cli)=tcga.cli$Samples

#
head(tcga.cli)
tcga.cli$DFI.time
tcga.cli=tcga.cli %>% drop_na(DFI.time)
head(tcga.cli)


table(tcga.cli$Stage)
tcga.cli$Stage[tcga.cli$Stage=='[Not Available]']=NA
tcga.cli$Stage=gsub('[ABC]','',tcga.cli$Stage)
tcga.cli$Stage=gsub('Stage ','',tcga.cli$Stage)
table(tcga.cli$Grade)
tcga.cli$Grade[tcga.cli$Grade%in%c('GB','GX')]=NA
dim(tcga.cli)

tcga.data=read.delim('origin_datas/TCGA/TCGA_OV_TPM.txt',row.names = 1,check.names = F)
tcga.data[1:4,1:4]
table(substr(colnames(tcga.data),14,15))
dim(tcga.data)
# colnames(tcga.data)=substr(colnames(tcga.data),1,12)


com.samples2=intersect(colnames(tcga.data),tcga.cli$Samples)

tcga.data=log2(tcga.data+1)
tcga.exp=tcga.data[,com.samples2]
range(tcga.exp)
dim(tcga.exp)

tcga.cli=tcga.cli[com.samples2,]
head(tcga.cli)
fivenum(tcga.cli$Age)
tcga.cli$Age1=ifelse(tcga.cli$Age>58,'>58','<=58')
tcga.cli$status=ifelse(tcga.cli$DFI==0,'NED','Relapse/Dead')
dim(tcga.cli)
#177


#01.###########
dir.create('results/01.scRNA')
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(dplyr)
library(harmony)

dir_name=list.files('origin_datas/GEO/GSE130000_RAW/')
sample.type=rep(c('primary','relapse'),c(4,2))
datalist=list()
for (i in 1:length(dir_name)){
  # i=1
  files = paste0("origin_datas/GEO/GSE130000_RAW/",dir_name[i])
  counts=fread(file = files,data.table = T,sep = '\t',check.names = F)
  counts[1:5,1:5]
  counts=data.frame(counts)
  rownames(counts)=counts[,1]
  counts=counts[,-1]
  datalist[[i]]<- CreateSeuratObject(counts=counts,project = stringr::str_split_fixed(dir_name[i],'_',3)[,1],min.cells = 3, min.features = 200)
  datalist[[i]]$Samples=stringr::str_split_fixed(dir_name[i],'_',3)[,1]
  datalist[[i]]$type=sample.type[i]
}
names(datalist)=stringr::str_split_fixed(dir_name,'_',3)[,1]

##########
for (i in 1:length(datalist)){
  sce <- datalist[[i]]
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")# 
  sce[["percent.Ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")# 
  datalist[[i]] <- sce
  rm(sce)
}
sce <- merge(x=datalist[[1]], y=datalist[2:length(datalist)])
#

raw_meta=sce@meta.data
raw_count <- table(raw_meta$Samples)
raw_count
sum(raw_count)#  26794
pearplot_befor<-VlnPlot(sce,group.by ='Samples',
                        features = c("nFeature_RNA", "nCount_RNA"),
                        pt.size = 0,
                        ncol = 2)
pearplot_befor
# ggsave('results/06.scRNA/pearplot_befor.pdf',pearplot_befor,height = 6,width = 15)
#
sce=subset(sce, subset=nFeature_RNA>200 & nFeature_RNA<1500 )

clean_meta=sce@meta.data
clean_count <- table(clean_meta$Samples)
clean_count
sum(clean_count)#23079
pearplot_after <- VlnPlot(sce,group.by ='Samples',
                          features = c("nFeature_RNA", "nCount_RNA"),
                          pt.size = 0,
                          ncol = 2)
pearplot_after
# ggsave('results/06.scRNA/pearplot_after.pdf',pearplot_after,height = 6,width = 15)
rm(datalist)
########
# sce <- NormalizeData(sce)
# sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
# sce <- ScaleData(sce, features = rownames(sce))

sce = SCTransform(sce, vars.to.regress = "percent.mt", verbose = F)
#
sce <- RunPCA(sce, features = VariableFeatures(sce))
colnames(sce@meta.data)
##
library(harmony)
sce = RunHarmony(sce, group.by.vars="Samples", max.iter.harmony=50, lambda=0.5,assay.use = "SCT")
#
pca.plot=ElbowPlot(sce,ndims = 50)+theme(text = element_text(family = 'Times',size = 12))
pca.plot
# ggsave('results/06.scRNA/pca.plot.pdf',pca.plot,height = 6,width = 6)

###
sce <- RunUMAP(sce, dims=1:20, reduction="harmony")
DimPlot(sce,group.by='Samples',reduction="umap",label =F,pt.size = 0.2)
# sce <- RunTSNE(sce, dims=1:15, reduction="harmony")
library(clustree)
sce <- FindNeighbors(sce, dims = 1:20, reduction="harmony")
#
sce <- FindClusters(object = sce,resolution = .1)
# DefaultAssay(sce) <- "RNA"
colnames(sce@meta.data)
length(table(sce@meta.data$seurat_clusters))
DimPlot(sce,group.by='seurat_clusters',reduction="umap",label =T,pt.size = 0.2)


marker <- data.frame(cluster = 0:5,cell = 0:5)
marker[marker$cluster %in% c(0,3,4),2] <- 'Epithelial cells'
marker[marker$cluster %in% c(1),2] <- 'Fibroblast'
marker[marker$cluster %in% c(2),2] <- 'Myeloid cells'
marker[marker$cluster %in% c(5),2] <- 'Endothelial cells'


marker
sce@meta.data$cell_type <- sapply(sce@meta.data$seurat_clusters,function(x){marker[x,2]})
my.cols=brewer.pal(11,"Set3")[5:8]
cell_type_umap=DimPlot(sce,group.by='cell_type',reduction="umap",label = F,
                       pt.size = 0.5,cols =my.cols)+
  # theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(text = element_text(family = 'Times',size = 12),legend.position = 'none',
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
cell_type_umap=LabelClusters(plot = cell_type_umap, id = 'cell_type', family = 'Times', size = 4, colour = "black", repel = T)
cell_type_umap
ggsave('results/cell_type_umap.pdf',cell_type_umap,height = 4.5,width = 5)

saveRDS(sce,file = 'results/01.scRNA/sce.rds')




feat_gene<-c("SPARCL1","VWF",
             'WFDC2','CD24',"KRT18","KRT19","EPCAM",
             "COL1A1","COL3A1","DCN",
             "APOE","SPP1","CD74"
)
length(feat_gene)



dotplot_gene_marker=DotPlot(sce, features=feat_gene,group.by = 'cell_type')+
  scale_color_gradientn(colors=c( "dodgerblue", "white", "orange", "firebrick1"))+
  theme_light()+coord_flip()+
  theme(text=element_text(family="Times"), 
        axis.text.x=element_text(angle=30, hjust=1, size=14, face="bold"), 
        axis.text.y=element_text(face="bold", size=14), axis.title.x=element_blank(), axis.title.y=element_blank())
dotplot_gene_marker


############
table(sce$Samples)
meta.data=sce@meta.data
head(meta.data)
table(meta.data$type)
bar = meta.data %>% group_by(type, cell_type) %>% dplyr::count()
state_label = c('primary','relapse')
bar$type = factor(bar$type, levels=state_label)
bar = bar %>% group_by(type) %>% mutate(percent=100*n/sum(n))

cell.prop1=ggplot(data=bar, aes(x=reorder(cell_type,percent), y=percent, fill=type,label = sprintf("%.2f", percent)))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(values=c("#E7298A","skyblue"))+theme_classic()+
  ylab("Percent(%)")+xlab('')+#coord_flip()+
  geom_text(position = position_dodge(width = 1),  size = 2.5)+
  theme(axis.text.x=element_text( size=12, face="bold",angle = 30,hjust = 1),
        legend.text=element_text(family = 'Times', size=12),legend.position = 'top',
        legend.title=element_blank(), text = element_text(family = 'Times'))
cell.prop1


chisq.test(table(sce$cell_type,sce$Samples))
cell_freq1=data.frame(t(prop.table(table(sce$cell_type,sce$Samples),margin=2)),sample.type)
cell_freq1
colnames(cell_freq1)<-c('Samples','cell_type','Freq','type')
cell.prop2=ggplot(cell_freq1,aes(x=Samples,y=Freq,fill=cell_type))+
  scale_fill_manual(values = my.cols)+
  facet_grid(~type,scales = 'free',space='free')+
  geom_bar(position = "fill",stat="identity")+
  xlab('')+ylab('Proportion')+theme_bw()+
  theme(text = element_text(family = 'Times',size=15),
        axis.text.x = element_text(angle = 30,hjust = 1),
        legend.text =element_text(family = 'Times'),
        legend.title = element_text(family = 'Times'))
cell.prop2

###########
table(sce$cell_type)
my.data=subset(sce,cell_type=='Epithelial cells')


ElbowPlot(my.data,ndims = 50)+theme(text = element_text(family = 'Times',size = 12))
my.data <- RunUMAP(my.data, dims=1:20, reduction="harmony")
DimPlot(my.data,group.by='type',reduction="umap",label = F,
        pt.size = 0.5,cols =c("#E7298A","skyblue"))
my.data <- FindNeighbors(my.data, dims = 1:20, reduction="harmony")
# #
my.data <- FindClusters(object = my.data,resolution = .1)
DimPlot(my.data,group.by='seurat_clusters',reduction="umap",label = F,
        pt.size = 0.5)
saveRDS(my.data,'results/01.scRNA/my.data.rds')

Idents(my.data)='type'
table(Idents(my.data))
markers_df <- FindMarkers(object = my.data, ident.1 = 'primary', ident.2 = 'relapse',
                          logfc.threshold = 0.1)
head(markers_df)
write.csv(markers_df,'results/01.scRNA/markers_df.csv')

markers_df <- read_csv("results/01.scRNA/markers_df.csv")
library(org.Hs.eg.db)
gene_list <- markers_df$avg_logFC[order(markers_df$avg_logFC,decreasing = T)]
names(gene_list) <- markers_df$...1[order(markers_df$avg_logFC,decreasing = T)]
head(gene_list)
res <- gseGO(gene_list,    # 
             ont = "BP",  OrgDb = org.Hs.eg.db, keyType = "SYMBOL",    # 
             pvalueCutoff = 0.05, pAdjustMethod = "BH")   # 

fgseaRes=res@result
write.xlsx(fgseaRes,'results/01.scRNA/GSEA_result.xlsx',overwrite = T)
head(fgseaRes)
# rm(res)
fgseaRes$sign<-ifelse(fgseaRes$NES>0,"Activated","Suppressed") ## 
fgseaRes=fgseaRes%>%group_by(sign)%>% slice_max(order_by = abs(NES),n = 5)
GSEA.plot=fgseaRes %>% 
  ggplot( aes(NES, fct_reorder(Description, NES),fill=-log10(pvalue))) + 
  geom_col(orientation='y') + 
  scale_fill_continuous(low="#8787FFFF", high="#FF707FFF", 
                        guide=guide_colorbar(reverse=TRUE)) +   
  ylab(NULL)+  ggtitle("primary vs relapse")+ ##
  theme_bw(base_size = 12)+
  theme(plot.title = element_text(hjust = 0.5),
        text=element_text(family = 'Times',size=14),
        axis.text.y = element_text(family = 'Times',size=18)) +
  scale_y_discrete(labels=function(a)str_wrap(a,width = 40))

GSEA.plot

pdf('results/01.scRNA/Fig1_1.pdf',height = 20,width = 15)
mg_merge_plot(mg_merge_plot(cell_type_umap,dotplot_gene_marker,
                            cell.prop1,cell.prop2,nrow=2,ncol=2,
                            labels = LETTERS[1:4]),
              GSEA.plot,nrow=2,labels = c('','E'))
dev.off()

#02.####################
dir.create('results/02.LASSO_STEP')

table(markers_df$p_val_adj<0.05 & abs(markers_df$avg_logFC)>0.25)
markers_df.fit=markers_df[markers_df$p_val_adj<0.05 & abs(markers_df$avg_logFC)>0.25,]
dim(markers_df.fit)

####
marker.cox=cox_batch(tcga.exp[markers_df$...1,tcga.cli$Samples],time = tcga.cli$DFI.time,event = tcga.cli$DFI)
marker.cox=na.omit(marker.cox)
table(marker.cox$p.value<0.05)
marker.cox.fit=marker.cox[marker.cox$p.value<0.05,]

pdf('results/02.LASSO_STEP/Fig2a.pdf',height = 9,width = 6)
bioForest(rt = marker.cox.fit[order(marker.cox.fit$HR),],col = c('orange','skyblue'))
dev.off()

pre.genes=rownames(marker.cox.fit)
length(pre.genes)
#27
tcga_model_data <- cbind(tcga.cli[, c("DFI.time", "DFI")],
                         t(tcga.exp[pre.genes, tcga.cli$Samples]))
colnames(tcga_model_data) <- gsub('-', '_', colnames(tcga_model_data))


######
library(glmnet)
set.seed(2024)

fit1=glmnet(as.matrix(tcga_model_data[,-c(1:2)])
            ,cbind(time=tcga_model_data$DFI.time,
                   status=tcga_model_data$DFI)
            ,family="cox"
            ,nlambda=100
            , alpha=1)

cv.fit<-cv.glmnet(as.matrix( tcga_model_data[,-c(1:2)])
                  ,cbind(time=tcga_model_data$DFI.time,
                         status=tcga_model_data$DFI)
                  ,family="cox"
                  ,nfolds = 10
                  ,nlambda=100
                  , alpha=1)

sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
print(cv.fit$lambda.min)
length(names(sig.coef))
#11
# 
# 
pdf('results/02.LASSO_STEP/LASSO.pdf',height = 4.5,width = 9,onefile = F)
par(mfrow=c(1,2))
plot(fit1)
plot(cv.fit)
dev.off()



######
fmla <- as.formula(paste0("Surv(DFI.time, DFI) ~"
                          ,paste0(names(sig.coef),collapse = '+')))
cox <- coxph(fmla, data =as.data.frame(tcga_model_data))
cox=step(cox)

lan <- coef(cox)
lan
length(lan)
paste0(round(lan, 3), '*', names(lan),collapse = '+')


gene.coef=data.frame(gene=names(lan),coef=as.numeric(lan))
gene.coef
gene.coef$Type=ifelse(gene.coef$coef>0,'Risk','Protective')
gene.coef.fig=ggplot(gene.coef, aes(x = coef, y = reorder(gene,coef), fill =Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#C75DAA", "#009B9F")) +
  labs(x = 'coefficient', y = "") +
  geom_text(aes(label = round(coef,2),hjust =2), data = subset(gene.coef, coef > 0))+ 
  geom_text(aes(label = round(coef,2), hjust = -1), data = subset(gene.coef, coef < 0))+  
  theme_bw()+ theme(text = element_text(family = 'Times',size = 14,face = 'bold'),legend.position = 'top')
gene.coef.fig
ggsave('results/02.LASSO_STEP/gene.coef.fig.pdf',gene.coef.fig,height = 4,width = 4)

gene.forest=ggforest(cox, data = tcga_model_data, 
                     main = "Hazardratio", fontsize =1.0, 
                     noDigits = 2)
gene.forest
ggsave('results/02.LASSO_STEP/gene.forest.pdf',gene.forest,height = 4,width = 6)


#03.###########
dir.create('results/03.model')
risktype.col=c("orange", "#009B9F")
###########
risk.tcga=as.numeric(lan%*%as.matrix(t(tcga_model_data[tcga.cli$Samples,names(lan)])))
tcga.risktype.cli=data.frame(tcga.cli,Riskscore=risk.tcga)
tcga.risktype.cli$Risktype=ifelse(tcga.risktype.cli$Riskscore>median(risk.tcga),'High','Low')

tcga.roc=ggplotTimeROC(time = tcga.risktype.cli$DFI.time,
                       status = tcga.risktype.cli$DFI,
                       score = tcga.risktype.cli$Riskscore,mks = c(1:5))
tcga.roc

tcga.km.DFI=ggsurvplot(fit=survfit( Surv(DFI.time/365, DFI) ~ Risktype,
                                    data = tcga.risktype.cli),
                       data=tcga.risktype.cli,
                       conf.int = F,pval = T,fun = "pct",risk.table =T, size = 0.7,
                       surv.median.line = 'hv',title='TCGA',
                       linetype = c("solid", "dashed","strata")[1],
                       palette = risktype.col,ggtheme = custom_theme(),
                       legend = c(0.8,0.85), # 
                       legend.title = "Risktype",legend.labs=c('High','Low'))
tcga.km.DFI=mg_merge_plot(tcga.km.DFI$plot,tcga.km.DFI$table,nrow=2,heights = c(2.5,1),align = 'v')
tcga.km.DFI


tcga.km.PFI=ggsurvplot(fit=survfit( Surv(PFI.time/365, PFI) ~ Risktype,
                                    data = tcga.risktype.cli),
                       data=tcga.risktype.cli,
                       conf.int = F,pval = T,fun = "pct",risk.table =T, size = 0.7,
                       surv.median.line = 'hv',title='TCGA',
                       linetype = c("solid", "dashed","strata")[1],
                       palette = risktype.col,ggtheme = custom_theme(),
                       legend = c(0.8,0.85), # 
                       legend.title = "Risktype",legend.labs=c('High','Low'))
tcga.km.PFI=mg_merge_plot(tcga.km.PFI$plot,tcga.km.PFI$table,nrow=2,heights = c(2.5,1),align = 'v')
tcga.km.PFI


tcga.km.OS=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                                   data = tcga.risktype.cli),
                      data=tcga.risktype.cli,
                      conf.int = F,pval = T,fun = "pct",risk.table =T, size = 0.7,
                      surv.median.line = 'hv',title='TCGA',
                      linetype = c("solid", "dashed","strata")[1],
                      palette = risktype.col,ggtheme = custom_theme(),
                      legend = c(0.8,0.85), # 
                      legend.title = "Risktype",legend.labs=c('High','Low'))
tcga.km.OS=mg_merge_plot(tcga.km.OS$plot,tcga.km.OS$table,nrow=2,heights = c(2.5,1),align = 'v')
tcga.km.OS






##GSE63885【OK】#######
# GSE63885=getGEOExpData('GSE63885')
# save(GSE63885,file = 'origin_datas/GEO/GSE63885.RData')
load('origin_datas/GEO/GSE63885.RData')
GSE63885.cli=GSE63885$Sample
GSE63885.cli=data.frame(Samples=GSE63885.cli$Acc,
                        event=GSE63885.cli$`clinical status at last follow-up (awd - alive with disease; dod - dead from disease; ned - no evidence of disease)`,
                        DFI.time=GSE63885.cli$`dfs - disease-free survival [days]`,
                        OS.time=GSE63885.cli$`os - overall survival [days]`)
table(GSE63885.cli$event)
GSE63885.cli$DFI.time
GSE63885.cli=GSE63885.cli%>%drop_na(DFI.time)
GSE63885.cli$DFI=ifelse(GSE63885.cli$event=='DOD',1,0)
GSE63885.cli$OS=ifelse(GSE63885.cli$event=='DOD',1,0)

GSE63885.exp=GSE63885$Exp$GPL570_54675_Data_col1
GSE63885.exp[1:5,1:5]
GSE63885.exp=exp_probe2symbol_v2(datExpr = GSE63885.exp,GPL = 'GPL570')
range(GSE63885.exp)

GSE63885_model_data <- cbind(GSE63885.cli[, c("DFI.time", "DFI")],
                             t(GSE63885.exp[names(lan), GSE63885.cli$Samples]))
colnames(GSE63885_model_data) <- gsub('-', '_', colnames(GSE63885_model_data))
GSE63885_model_data


fmla1=as.formula(paste0("Surv(DFI.time, DFI) ~"
                        ,paste0(names(lan),collapse = '+')))
GSE63885.cox <- coxph(fmla1, data =as.data.frame(GSE63885_model_data))
GSE63885.lan=coef(GSE63885.cox)


risk.GSE63885=as.numeric(GSE63885.lan%*%as.matrix(t(GSE63885_model_data[GSE63885.cli$Samples,names(GSE63885.lan)])))
GSE63885.risktype.cli=data.frame(GSE63885.cli,Riskscore=risk.GSE63885)
GSE63885.risktype.cli$Risktype=ifelse(GSE63885.risktype.cli$Riskscore>median(risk.GSE63885),'High','Low')

GSE63885.roc=ggplotTimeROC(GSE63885.risktype.cli$DFI.time,
                           GSE63885.risktype.cli$DFI,
                           GSE63885.risktype.cli$Riskscore,mks = c(1:5))
GSE63885.roc



GSE63885.km.DFI=ggsurvplot(fit=survfit( Surv(DFI.time/365, DFI) ~ Risktype,
                                        data = GSE63885.risktype.cli),
                           data=GSE63885.risktype.cli,
                           conf.int = F,pval = T,fun = "pct",risk.table =T, size = 0.7,
                           surv.median.line = 'hv',title='GSE63885',
                           linetype = c("solid", "dashed","strata")[1],
                           palette = risktype.col,ggtheme = custom_theme(),
                           legend = c(0.8,0.85), # 
                           legend.title = "Risktype",legend.labs=c('High','Low'))
GSE63885.km.DFI=mg_merge_plot(GSE63885.km.DFI$plot,GSE63885.km.DFI$table,nrow=2,heights = c(2.5,1),align = 'v')
GSE63885.km.DFI


GSE63885.km.OS=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                                        data = GSE63885.risktype.cli),
                           data=GSE63885.risktype.cli,
                           conf.int = F,pval = T,fun = "pct",risk.table =T, size = 0.7,
                           surv.median.line = 'hv',title='GSE63885',
                           linetype = c("solid", "dashed","strata")[1],
                           palette = risktype.col,ggtheme = custom_theme(),
                           legend = c(0.8,0.85), # 
                           legend.title = "Risktype",legend.labs=c('High','Low'))
GSE63885.km.OS=mg_merge_plot(GSE63885.km.OS$plot,GSE63885.km.OS$table,nrow=2,heights = c(2.5,1),align = 'v')
GSE63885.km.OS

pdf('results/03.model/Fig3.pdf',height = 10,width = 15)
mg_merge_plot(tcga.roc,tcga.km.DFI,tcga.km.OS,
              GSE63885.roc,GSE63885.km.DFI,ncol=3,nrow=2,labels = LETTERS[1:5])
dev.off()

#04.#########
dir.create('results/04.gene_expression')
# plotMutiBar(table(tcga.risktype.cli$status,tcga.risktype.cli$Risktype))
tcga.gene.df=data.frame(Risktype=tcga.risktype.cli$Risktype,
                        t(tcga.exp[names(lan),tcga.risktype.cli$Samples]))
tcga.gene.df=melt(tcga.gene.df)
head(tcga.gene.df)
fig4a=ggplot(tcga.gene.df,aes(x=variable,y=value,fill=Risktype))+geom_boxplot()+
  scale_fill_manual(values = risktype.col)+
  stat_compare_means(aes(group=Risktype), label = "p.signif", method = 'wilcox.test')+
  ylab('Expression')+xlab('')+ggtitle('TCGA')+theme_bw()+
  theme(text = element_text(family = 'Times',size=12))
fig4a

# plotMutiBar(table(GSE63885.risktype.cli$event,GSE63885.risktype.cli$Risktype))
GSE63885.gene.df=data.frame(Risktype=GSE63885.risktype.cli$Risktype,
                            t(GSE63885.exp[names(lan),GSE63885.risktype.cli$Samples]))

GSE63885.gene.df=melt(GSE63885.gene.df)
head(GSE63885.gene.df)
fig4b=ggplot(GSE63885.gene.df,aes(x=variable,y=value,fill=Risktype))+geom_boxplot()+
  scale_fill_manual(values = risktype.col)+
  stat_compare_means(aes(group=Risktype), label = "p.signif", method = 'wilcox.test')+
  ylab('Expression')+xlab('')+ggtitle('GSE63885')+theme_bw()+
  theme(text = element_text(family = 'Times',size=12))
fig4b


fig4c=DotPlot(sce, features=names(lan),group.by = 'cell_type')+
  scale_color_gradientn(colors=c( "dodgerblue", "white", "orange", "firebrick1"))+
  theme_light()+coord_flip()+
  theme(text=element_text(family="Times"), 
        axis.text.x=element_text(angle=45, hjust=1, size=14, face="bold"), 
        axis.text.y=element_text(face="bold", size=14), axis.title.x=element_blank(), axis.title.y=element_blank())
fig4c

fig4d=DotPlot(sce, features=names(lan),group.by = 'type')+
  scale_color_gradientn(colors=c( "dodgerblue", "white", "orange", "firebrick1"))+
  theme_light()+coord_flip()+
  theme(text=element_text(family="Times"), 
        axis.text.x=element_text(size=14, face="bold"),
        axis.text.y=element_text(face="bold", size=14), axis.title.x=element_blank(), axis.title.y=element_blank())
fig4d

fig4=mg_merge_plot(mg_merge_plot(fig4a,fig4b,common.legend = T,labels = LETTERS[1:2]),
                   mg_merge_plot(fig4c,fig4d,align = 'h',labels = LETTERS[3:4]),
                   nrow = 2,ncol=1)
ggsave('results/04.gene_expression/Fig4.pdf',fig4,height = 10,width = 9)


#05.##############
dir.create('results/05.clinical_KM')
head(tcga.risktype.cli)
tcga_cox_datas=tcga.risktype.cli
colnames(tcga_cox_datas)
table(tcga_cox_datas$Age1)

table(tcga_cox_datas$Stage)
tcga_cox_datas$Stage[tcga_cox_datas$Stage=='I'|tcga_cox_datas$Stage=='II']<-'I+II'
# tcga_cox_datas$Stage[tcga_cox_datas$Stage=='III'|tcga_cox_datas$Stage=='IV']<-'III+IV'

table(tcga_cox_datas$Grade)
# tcga_cox_datas$Grade[tcga_cox_datas$Grade=='G1'|tcga_cox_datas$Grade=='G2']<-'G1+G2'
tcga_cox_datas$Grade[tcga_cox_datas$Grade=='G3'|tcga_cox_datas$Grade=='G4']<-'G3+G4'



head(tcga_cox_datas)
table(tcga_cox_datas$Age1)
fig5a=ggsurvplot(fit=survfit( Surv(DFI.time/365, DFI) ~ Risktype,
                        data = tcga_cox_datas[which(tcga_cox_datas$Age1=='<=58'),]),
           data=tcga_cox_datas[which(tcga_cox_datas$Age1=='<=58'),],
           conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,
           title='Age<=58',ggtheme=custom_theme(),
           linetype = c("solid", "dashed","strata")[1],
           legend = c('top', 'bottom', 'left', 'right', 'none')[1],
           legend.title = "Risktype",
           legend.labs = c("High","Low"))


fig5b=ggsurvplot(fit=survfit( Surv(DFI.time/365, DFI) ~ Risktype,
                        data = tcga_cox_datas[which(tcga_cox_datas$Age1=='>58'),]),
           data=tcga_cox_datas[which(tcga_cox_datas$Age1=='>58'),],
           conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,
           title='Age>58',ggtheme=custom_theme(),
           linetype = c("solid", "dashed","strata")[1],
           legend = c('top', 'bottom', 'left', 'right', 'none')[1],
           legend.title = "Risktype",
           legend.labs = c("High","Low"))


fig5c=ggsurvplot(fit=survfit( Surv(DFI.time/365, DFI) ~ Risktype,
                        data = tcga_cox_datas[which(tcga_cox_datas$Stage=='I+II'),]),
           data=tcga_cox_datas[which(tcga_cox_datas$Stage=='I+II'),],
           conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,
           title='Stage I+II',ggtheme=custom_theme(),
           linetype = c("solid", "dashed","strata")[1],
           legend = c('top', 'bottom', 'left', 'right', 'none')[1],
           legend.title = "Risktype",
           legend.labs = c("High","Low"))


fig5d=ggsurvplot(fit=survfit( Surv(DFI.time/365, DFI) ~ Risktype,
                        data = tcga_cox_datas[which(tcga_cox_datas$Stage=='III'),]),
           data=tcga_cox_datas[which(tcga_cox_datas$Stage=='III'),],
           conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,
           title='Stage III',ggtheme=custom_theme(),
           linetype = c("solid", "dashed","strata")[1],
           legend = c('top', 'bottom', 'left', 'right', 'none')[1],
           legend.title = "Risktype",
           legend.labs = c("High","Low"))


fig5e=ggsurvplot(fit=survfit( Surv(DFI.time/365, DFI) ~ Risktype,
                        data = tcga_cox_datas[which(tcga_cox_datas$Grade=='G2'),]),
           data=tcga_cox_datas[which(tcga_cox_datas$Grade=='G2'),],
           conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,
           title='G2',ggtheme=custom_theme(),
           linetype = c("solid", "dashed","strata")[1],
           legend = c('top', 'bottom', 'left', 'right', 'none')[1],
           legend.title = "Risktype",
           legend.labs = c("High","Low"))

fig5f=ggsurvplot(fit=survfit( Surv(DFI.time/365, DFI) ~ Risktype,
                        data = tcga_cox_datas[which(tcga_cox_datas$Grade=='G3+G4'),]),
           data=tcga_cox_datas[which(tcga_cox_datas$Grade=='G3+G4'),],
           conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,
           title='G3+G4',ggtheme=custom_theme(),
           linetype = c("solid", "dashed","strata")[1],
           legend = c('top', 'bottom', 'left', 'right', 'none')[1],
           legend.title = "Risktype",
           legend.labs = c("High","Low"))

fig5a=mg_merge_plot(fig5a$plot,fig5a$table,nrow=2,heights = c(2.5,1),align = 'v')
fig5b=mg_merge_plot(fig5b$plot,fig5b$table,nrow=2,heights = c(2.5,1),align = 'v')
fig5c=mg_merge_plot(fig5c$plot,fig5c$table,nrow=2,heights = c(2.5,1),align = 'v')
fig5d=mg_merge_plot(fig5d$plot,fig5d$table,nrow=2,heights = c(2.5,1),align = 'v')
fig5e=mg_merge_plot(fig5e$plot,fig5e$table,nrow=2,heights = c(2.5,1),align = 'v')
fig5f=mg_merge_plot(fig5f$plot,fig5f$table,nrow=2,heights = c(2.5,1),align = 'v')

fig5=mg_merge_plot(fig5a,fig5c,fig5e,fig5b,fig5d,fig5f,
                   ncol=3,nrow=2,labels = LETTERS[1:6])
ggsave('results/05.clinical_KM/Fig5.pdf',fig5,height = 10,width = 15)


#06.############
dir.create('results/06.pathway')
#############
tcga.geneList.risktype=getGeneFC(gene.exp=tcga.exp[,tcga.risktype.cli$Samples],
                                 group=tcga.risktype.cli$Risktype,ulab='High',dlab='Low')
h.all.gmt<-read.gmt("h.all.v7.5.1.entrez.gmt")
set.seed(777)
tcga.geneList.risktype.gsea<-GSEA(tcga.geneList.risktype,TERM2GENE = h.all.gmt,seed=T)
tcga.geneList.risktype.gsea.res=tcga.geneList.risktype.gsea@result
write.xlsx(tcga.geneList.risktype.gsea.res,'results/06.pathway/TCGA.risktypeGSEA.res.xlsx',overwrite = T)

library(dotplotGsea)
risktype.gsea.dotplot=dotplotGsea(data = tcga.geneList.risktype.gsea,order.by = 'NES')

fig6a=risktype.gsea.dotplot$plot+
  theme(text = element_text(family = 'Times'))+
  ggtitle('High risk VS Low risk')+
  scale_y_discrete(labels=function(x) str_remove(x,"HALLMARK_"))
fig6a



save.image(file = 'project.RData')
##IC50############
load('results/tcga_durg_ic50_res.RData')
tcga_durg_ic50_res[1:5,1:5]
colnames(tcga_durg_ic50_res)[1]='Cisplatin'

IC50.mat=cbind.data.frame(Riskscore=tcga.risktype.cli$Riskscore,
                          tcga_durg_ic50_res[tcga.risktype.cli$Samples,])

library(ggcorrplot)
library(psych)
IC50_RS_cor <- corr.test(x =IC50.mat$Riskscore,
                         y = IC50.mat[,-1],
                         method = "spearman",adjust = "BH",ci = F)


IC50_RS_cor_res=data.frame(drugs=colnames( IC50.mat[,-1]))
IC50_RS_cor_res$cor<-as.numeric(IC50_RS_cor$r)
IC50_RS_cor_res$p.adj<-as.numeric(IC50_RS_cor$p.adj)
head(IC50_RS_cor_res)
table(IC50_RS_cor_res$p.adj<0.05,abs(IC50_RS_cor_res$cor)>0.3)
IC50_RS_cor_res=IC50_RS_cor_res[IC50_RS_cor_res$p.adj<0.05 & abs(IC50_RS_cor_res$cor)>0.3,]
IC50_RS_cor_res=IC50_RS_cor_res[order(IC50_RS_cor_res$cor),]
head(IC50_RS_cor_res)

library(rcartocolor)
fig6b=ggplot(data=IC50_RS_cor_res,aes(x=cor,y=reorder(drugs,cor), color = -log10(p.adj))) +
  geom_point(aes(size=abs(cor)),show.legend = F) +
  scale_colour_gradient(low = "#80B1D3",high = 'orange')+
  geom_segment(aes(yend=drugs,xend=0),size=.5) +
  labs(x='spearman Correlation',y='Drugs')+theme_bw()+
  theme(text = element_text(family = 'Times',size=15),legend.position = 'top',
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


IC50.df=data.frame(Risktype=tcga.risktype.cli$Risktype,
                   tcga_durg_ic50_res[tcga.risktype.cli$Samples,IC50_RS_cor_res$drugs])
head(IC50.df)
IC50.df= reshape2::melt(IC50.df)
head(IC50.df)


fig6c=ggplot(IC50.df,aes(x=Risktype,y=value,fill=Risktype))+
  geom_boxplot()+facet_wrap(~variable,scales = 'free',nrow = 2,ncol = 3)+
  stat_compare_means(aes(group=Risktype), label = "p.signif",
                     method = 'wilcox.test')+
  scale_fill_manual(values =risktype.col)+
  xlab('Risktype')+ylab('IC50')+
  theme_bw()+theme(text = element_text(family = 'Times',size = 14),
                   legend.position = 'none',
                   axis.text.x = element_text(color = "black", size = 10),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())

fig6=mg_merge_plot(fig6a,mg_merge_plot(fig6b,fig6c,labels = c('B','C')),
                   nrow=2,ncol=1,labels = c('A',''),heights = c(1.5,1))
ggsave('results/06.pathway/Fig6.pdf',fig6,height = 12,width = 12)


