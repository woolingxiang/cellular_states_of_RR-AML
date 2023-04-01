###==================================================================================
### Integration of scRNA-seq datasets of AML patients
setwd('/public/workspace/AML_Leukemia')
options(stringsAsFactors = F)
source('~/R/functions.R')
source('~/R/requiredPcg.R')

path = Sys.glob('/public/workspace/AML_Leukemia/00-DATA/AMLData/*.rds')
aml.list = lapply(path,function(i) {
	id = gsub('*.*/','',gsub('.label.rds','',i))
	tmp = readRDS(i)
	tmp$'source' = rep(id,times=ncol(tmp))
	tmp = NormalizeData(tmp, verbose = FALSE)
	tmp = FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 4000, verbose = FALSE) 
	return(tmp)})
anchors = FindIntegrationAnchors(object.list = aml.list,dims = 1:30,anchor.features=4000)
integrated = IntegrateData(anchorset = anchors, dims = 1:30)
integrated = ScaleData(integrated, verbose = FALSE,features=rownames(integrated)) 
integrated = RunPCA(integrated, npcs = 50, verbose = FALSE,features=VariableFeatures(object = integrated))
integrated = JackStraw(integrated,num.replicate = 100)
integrated = ScoreJackStraw(integrated,dims=1:20)
integrated = FindNeighbors(integrated, dims = 1:40)
integrated = FindClusters(integrated, resolution = 3) 
integrated = RunTSNE(integrated,dims=1:40)
integrated = RunUMAP(integrated,dims=1:40)
saveRDS(integrated,file='AML.integrated.rds')




###==================================================================================
### Integration of  scRNA-seq datasets of healthy donors
setwd('/public/workspace/AML_Leukemia')
options(stringsAsFactors = F)
source('~/R/functions.R')
source('~/R/requiredPcg.R')
path = Sys.glob('/public/workspace/AML_Leukemia/00-DATA/HCData/*.rds')
dat.list = lapply(path,function(i) {
	id = gsub('*.*/','',gsub('.rds','',i))
	tmp = readRDS(i)
	tmp$'source' = rep(id,times=ncol(tmp))
	tmp = NormalizeData(tmp, verbose = FALSE)
	tmp = FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 4000, verbose = FALSE) 
	return(tmp) })
anchors = FindIntegrationAnchors(object.list = dat.list,dims = 1:30,anchor.features=4000)
integrated = IntegrateData(anchorset = anchors, dims = 1:30)
integrated = ScaleData(integrated, verbose = FALSE,features=rownames(integrated))
integrated = RunPCA(integrated, npcs = 50, verbose = FALSE,features=VariableFeatures(object = integrated))
integrated = JackStraw(integrated,num.replicate = 100)
integrated = ScoreJackStraw(integrated,dims=1:20)
integrated = FindNeighbors(integrated, dims = 1:35)
integrated = FindClusters(integrated, resolution = 3) 
integrated = RunTSNE(integrated,dims=1:35)
integrated = RunUMAP(integrated,dims=1:35)
saveRDS(integrated,file='HC.integrated.rds')




###==================================================================================
### Integration of scRNA-seq datasets of healthy donors and aml patients
setwd('/public/workspace/AML_Leukemia')
options(stringsAsFactors = F)
source('~/R/functions.R')
source('~/R/requiredPcg.R')
path = c('AML.integrated.rds','HC.integrated.rds')
dat.list = lapply(path,function(i) {
	tmp = readRDS(i)
	tmp = FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 4000, verbose = FALSE) 
	return(tmp) })
integrated = ScaleData(integrated, verbose = FALSE,features=rownames(integrated))
integrated = RunPCA(integrated, npcs = 50, verbose = FALSE,features=VariableFeatures(object = integrated))
integrated = JackStraw(integrated,num.replicate = 100)
integrated = ScoreJackStraw(integrated,dims=1:20)
integrated = FindNeighbors(integrated, dims = 1:35)
integrated = FindClusters(integrated, resolution = 5) 
integrated = RunTSNE(integrated,dims=1:35)
integrated = RunUMAP(integrated,dims=1:35)
integrated$source2 = ifelse(substr(integrated$source,1,3)=='AML', integrated$source, 'HC')
integrated$source3 = ifelse(substr(integrated$source,1,3)=='AML', 'AML', 'HC')
saveRDS(integrated,file='/public/workspace/AML_Leukemia/AML_HC.integrated.rds')







###==================================================================================
### Identification of Leukemia-like Cells States
path = '/public/workspace/AML_Leukemia/Figure1'
if(!file.exists(path)) dir.create(path)
setwd(path)

options(stringsAsFactors = F)
source('~/R/functions.R')
source('~/R/requiredPcg.R')
dat = readRDS('/public/workspace/AML_Leukemia/AML_HC.integrated.rds')
sum = dat@meta.data %>% mutate(sourceBin = ifelse(dat$source2 == 'AML',1,0)) %>% group_by(seurat_clusters) %>% 
		summarise(AML.pct = sum(sourceBin)/n()*100, HC.pct = 1-sum(sourceBin)/n()*100)
m.cluster = sum$seurat_clusters[which(sum$HC.pct<60)]
dat$Malignant_cluster = ifelse(dat$seurat_clusters %in% m.cluster,'Leukemia','NormalLike')
dat$Malignant_cluster = ifelse(dat$source2=='HC' | dat$CellType %in% c('B','NaiveT','CTL','CTL_NK','Ery'),'NormalLike',dat$Malignant_cluster)

saveRDS(dat, file = '/public/workspace/AML_Leukemia/AML_HC.integrated.rds')




###==================================================================================
### Identification of cellular states of leukemia-like cells 
options(stringsAsFactors = F)
source('~/R/functions.R')
source('~/R/requiredPcg.R')

dat = readRDS('/public/workspace/AML_LeukemiaAML_HC.integrated.rds')
dat = subset(dat,subset=(Malignant_cluster=='Leukemia'))
monocle = newMonocleAnalysis(dat,'Monocle.pdf',get_color_scheme('clusters'))
monocle$State2 = ifelse(monocle$State==9, 'QSC', monocle$State==1 | monocle$State==2,'PSP',
				ifelse(monocle$State==3, 'GMP', ifelse(monocle$State=='8','PG',
					ifelse(monocle$State==4, 'ProMono', 'Mono'))))
dat$State2 = monocle$State
save(monocle,file='monocle.RData')
saveRDS(dat,file='/public/workspace/AML_Leukemia/AML_Leukemia.integrated.rds')



###==================================================================================
### Labeltransform reference
order = c('QSC','PSP','GMP','PG','ProMono','Mono')
tree = readRDS('/public/workspace/AML_Leukemia/AML_HC.integrated.rds')
mye = readRDS('/public/workspace/AML_Leukemia/AML_Leukemia.integrated.rds')
tmp = merge(data.frame(tree@meta.data,ID=colnames(tree)),data.frame(mye@meta.data,ID=colnames(mye)),by='ID',all.x=TRUE) 
rownames(tmp) = tmp$ID

tree$State = tmp[intersect(colnames(tree),rownames(tmp)),'State2']
tree$State = ifelse(!is.na(tree$State),tree$State,ifelse(is.na(tree$State) & tree$source3 == 'AML','Unknown',tree$pred.fine))
tree$State2 = ifelse(tree$State %in% c("Hematopoietic stem cells_CD133+ CD34dim","Common myeloid progenitors",
										"Granulocyte/monocyte progenitors","Hematopoietic stem cells_CD38- CD34+"),'HSC_Prog',
				ifelse(tree$State %in% c("Monocytes","Colony Forming Unit-Granulocytes","Colony Forming Unit-Monocytes",
										 "Granulocytes (Neutrophilic Metamyelocytes)","Granulocytes (Neutrophils)"),'Myeloid',tree$State))
tree = subset(tree,subset=(State2 %in% c('HSC_Prog','Myeloid',order)))
saveRDS(tree,file='/public/workspace/AML_Leukemia/AML_HC_labeltransform_tree.rds')




###==================================================================================
### DEGS between QSC and PSP
options(stringsAsFactors = F)
source('~/R/functions.R')
source('~/R/requiredPcg.R')

dat=readRDS("/public/workspace/AML_Leukemia/AML_Leukemia.integrated.rds")
DefaultAssay(dat)="RNA"
Idents(dat) = factor(dat$State2, levels = unique(dat$State2))
dat.markers = FindAllMarkers(dat, only.pos = TRUE,test.use="MAST",min.pct = 0.25, logfc.threshold = 0.25)
dat.markers$mean.1=NA
dat.markers$mean.2=NA
dat.markers$median.1=NA
dat.markers$median.2=NA
for(i in 1:nrow(dat.markers)){
    case.label=rownames(dat@meta.data)[dat$State2==as.character(dat.markers$cluster[i])]
    control.label=rownames(dat@meta.data)[dat$State2!=as.character(dat.markers$cluster[i])]
    dat.markers$mean.1[i]=mean(dat@assays$RNA@data[dat.markers$gene[i],case.label])
    dat.markers$mean.2[i]=mean(dat@assays$RNA@data[dat.markers$gene[i],control.label])
    dat.markers$median.1[i]=median(dat@assays$RNA@data[dat.markers$gene[i],case.label])
    dat.markers$median.2[i]=median(dat@assays$RNA@data[dat.markers$gene[i],control.label]) }
write.table(dat.markers,"allmarker.mean.median.txt",sep="\t",quote=F)

marker = read.table("allmarker.mean.median.txt",sep="\t",stringsAsFactors=F,header=T)
marker$avg_logFC = round(marker$avg_logFC,2)
marker$pct.ratio = marker$pct.1/marker$pct.2
marker9 = marker%>%filter(cluster%in%c("QSC"),avg_logFC>=0.5,pct.1>=0.4,pct.ratio>=2)
marker1_2 = marker%>%filter(cluster%in%c("PSP"),avg_logFC>=0.5,pct.1>=0.4,pct.ratio>=2)
del = intersect(marker9$gene,marker1_2$gene)

select = setdiff(marker9$gene,del)
final=as.matrix(dat@assays$RNA@data[select,])
a = rownames(dat@meta.data)[dat$State2=="QSC"]
mean=apply(final[,a],2,function(x){mean(x)})
oa=rev(order(mean))

select = setdiff(marker1_2$gene,del)
final=as.matrix(dat@assays$RNA@data[select,])
b = rownames(dat@meta.data)[dat$State2=="PSP"]
mean=apply(final[,b],2,function(x){mean(x)})
ob = order(mean)
select = c(del,setdiff(marker9$gene,del),setdiff(marker1_2$gene,del))
cm = cbind(as.matrix(dat@assays$RNA@data[select,a[oa]]),as.matrix(dat@assays$RNA@data[select,b[ob]]))

pdf("Figure3_A.pdf",useDingbats=F)
annotation_col=data.frame(Type=c(rep("QSC",length(a)),rep("PSP",length(b))))
rownames(annotation_col)=colnames(cm)
ann_colors = list(Type= c('QSC'="darkred",'PSP'='lightblue')) 
bk = c(seq(-6,6,by=0.02))
pheatmap(cm,breaks=bk,scale="row",color = c(colorRampPalette(colors = c("lightblue","white"))(length(bk)/2),
    colorRampPalette(colors = c("white","darkred"))(length(bk)/2)),cluster_cols=FALSE,cluster_rows=FALSE,border=NA,
    annotation=annotation_col,show_rownames=TRUE,show_colnames=FALSE,annotation_colors=ann_colors,fontsize=5,
    gaps_col=c(4009),gaps_row=c(5,39))
dev.off()




###==================================================================================
### TF of QSC and PSP
dat = read.table("../TF/myeloid.auc_mtx.tsv",row.names=1,header=T,sep="\t",stringsAsFactors=F) 
dat2 = readRDS("/public/workspace/AML_Leukemia/AML_Leukemia.integrated.rds")

a1 = colnames(dat2)[dat2$State2 %in% c("QSC")]
a2 = colnames(dat2)[dat2$State2 %in% c("PSP")]


group1 = dat[which(rownames(dat) %in% a1),]
group2 = dat[which(rownames(dat) %in% a2),]
cm = rbind(group1,group2)
p = apply(cm,2,function(x){wilcox.test(x[1:nrow(group1)],x[(nrow(group1)+1):nrow(cm)])[[3]]})
fc = apply(cm,2,function(x){mean(as.numeric(x[1:nrow(group1)]))/mean(as.numeric(x[(nrow(group1)+1):nrow(cm)]))})
mean.9 = apply(cm,2,function(x){mean(as.numeric(x[1:nrow(group1)]))})
mean.1_2 = apply(cm,2,function(x){mean(as.numeric(x[(nrow(group1)+1):nrow(cm)]))})

x = data.frame(label=names(fc),Mean.9=mean.9,Mean.1_2=mean.1_2, FC=as.numeric(fc), P.Value=as.numeric(p))
x$P.Value[x$P.Value==0]=1.381878e-300
rownames(x)=gsub("\\...","",rownames(x))

DefaultAssay(dat2) = "RNA"

tmp = dat2@assays$RNA@data[intersect(rownames(x),rownames(dat2)),a1]
tmp.n9 = apply(tmp,1,function(x){length(which(x>0))/ncol(tmp)})

tmp = dat2@assays$RNA@data[intersect(rownames(x),rownames(dat2)),a2]
tmp.n1_2 = apply(tmp,1,function(x){length(which(x>0))/ncol(tmp)})

select=intersect(names(which(tmp.n9<0.05)),names(which(tmp.n1_2<0.05)))
select=setdiff(names(tmp.n9),select)
x=x[select,]

logFCcut = 0.25
pvalCut = 1.30103 
logFCcut2 = 0.38
pvalCut2 = 12
logFCcut3=1
pvalCut3=20

n1 = length(x[, 1])
cols = rep("grey", n1)
names(cols)= rownames(x)
cols[-log10(x$P.Value) > pvalCut & log2(x$FC) >logFCcut]= "#9C9C9C"
cols[-log10(x$P.Value) > pvalCut2 & log2(x$FC) > logFCcut2]= "#ED4F4F"
cols[-log10(x$P.Value) > pvalCut & log2(x$FC) < -logFCcut]= "#B2DF8A"
cols[-log10(x$P.Value) > pvalCut2 & log2(x$FC) < -logFCcut2]= "#329E3F"
color_transparent = adjustcolor(cols, alpha.f = 0.5)
x$color_transparent = color_transparent

n1 = length(x[, 1])
size = rep(1, n1)
size[-log10(x$P.Value) > pvalCut & log2(x$FC) > logFCcut]= 2
size[-log10(x$P.Value) > pvalCut2 & log2(x$FC) > logFCcut2]= 4
size[-log10(x$P.Value) > pvalCut3 & log2(x$FC) > logFCcut3]= 6
size[-log10(x$P.Value) > pvalCut & log2(x$FC) < -logFCcut]= 2
size[-log10(x$P.Value) > pvalCut2 & log2(x$FC) < -logFCcut2]= 4
size[-log10(x$P.Value) > pvalCut3 & log2(x$FC) < -logFCcut3]= 6

p1 = ggplot(data=x, aes(log2(FC), -log10(P.Value), label = label)) +
  		geom_point(alpha = 0.6, size = size, colour = x$color_transparent) +
		labs(x=bquote(~Log[2]~"(fold change)"), y=bquote(~-Log[10]~italic("P-value")), title="") + 
  		scale_x_continuous(
		    breaks = c( -1, -0.38, -logFCcut, 0, 0.38, logFCcut, 1),
		    labels = c( -1, -0.38, -logFCcut, 0, 0.38, logFCcut, 1),
		    limits = c(-1.5, 1.5) ) +
  		geom_vline(xintercept = c(-logFCcut, logFCcut), color="grey91", linetype="longdash", lwd = 0.5) + 
  		geom_hline(yintercept = pvalCut, color="grey91", linetype="longdash", lwd = 0.5) +
  		geom_vline(xintercept = c(-logFCcut2, logFCcut2), color="grey91", linetype="longdash", lwd = 0.5) +
  		geom_hline(yintercept = pvalCut2, color="grey91", linetype="longdash", lwd = 0.5)+
  		theme_bw(base_size = 12) +
  		theme(panel.grid=element_blank())


plot = p1 + geom_text_repel(aes(x = log2(FC), y = -log10(P.Value), label = ifelse(log2(FC) > log2(1.3) & -log10(P.Value) > 12 , rownames(x),"")),
        					colour="darkred", size = 3, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"))+ 
    		geom_text_repel(aes(x = log2(FC), y = -log10(P.Value), label = ifelse(log2(FC) < (-log2(1.3)) & -log10(P.Value) > 12 , rownames(x),"")),
        					colour="darkgreen", size = 3, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"))

pdf("Figure3_B.pdf")
print(plot)
dev.off()




###==================================================================================
### metabolism analysis of QSC and PSP
options(stringsAsFactors = F)
source('~/R/functions.R')
source('~/R/requiredPcg.R')
source('/public/workspace/lily/software/SingleCellMetabolic/utils.R')
source('/public/workspace/lily/software/SingleCellMetabolic/runGSEA_preRank.R')

datpath = '/public/workspace/AML_Leukemia/AML_Leukemia.integrated.rds'

### cluster
tmp_data = readRDS(datpath) %>% FindNeighbors(dims=1:20) %>% FindClusters(resolution=2.5) %>%
		RunUMAP(dims=1:20) %>% RunTSNE(dims=1:20)
saveRDS(tmp_data,file = datpath)
all_data = as.matrix(tmp_data@assays$RNA@data)

info = data.frame()
a = names(table(tmp_data$seurat_clusters))
for(i in 1:length(a)){
    n = table(tmp_data$State2[tmp_data$seurat_clusters==a[i]])/sum(table(tmp_data$State2[tmp_data$seurat_clusters==a[i]]))
    label = names(which(n>=0.5))
    if(length(label)!=0){
        tmp = data.frame(cluster=a[i],label=label)
        info = rbind(info,tmp)
    }
}

info2 = data.frame(cluster=tmp_data$seurat_clusters)
info2$cluster = as.character(info2$cluster)
info2 = left_join(info2,info,by="cluster")
tmp_data$new=info2$label

cell_type = tmp_data$seurat_clusters# because all cell is tumor ,so cell type is seurat cluster 
tumor = unname(tmp_data$seurat_clusters)
col_data = data.frame(tumor=tumor,cellType=as.character(cell_type),row.names=names(cell_type))
pathways = gmtPathways("/public/workspace/lily/software/SingleCellMetabolic/Data/KEGG_metabolism.gmt")
metabolics = unique(as.vector(unname(unlist(pathways))))
row_data = data.frame(metabolic=rep(FALSE,nrow(all_data)),row.names = rownames(all_data))
row_data[rownames(row_data)%in%metabolics,"metabolic"]=TRUE
sce = SingleCellExperiment(assays = all_data,colData = col_data,rowData = row_data)
selected_tumor_sce = sce 
selected_tumor_metabolic_sce = sce[rowData(sce)$metabolic,] 


### scRNA_pathway_activity 
pathway_file = "/public/workspace/lily/software/SingleCellMetabolic/Data/KEGG_metabolism.gmt"
pathways = gmtPathways(pathway_file)
pathway_names = names(pathways)
all_cell_types = as.vector(selected_tumor_metabolic_sce$cellType)
cell_types = unique(all_cell_types)

gene_pathway_number = num_of_pathways(pathway_file,rownames(selected_tumor_metabolic_sce)[rowData(selected_tumor_metabolic_sce)$metabolic])
set.seed(123)
normalization_method = "Deconvolution"

##Calculate the pathway activities
#mean ratio of genes in each pathway for each cell type
mean_expression_shuffle = matrix(NA,nrow=length(pathway_names),ncol=length(cell_types),dimnames = list(pathway_names,cell_types))
mean_expression_noshuffle = matrix(NA,nrow=length(pathway_names),ncol=length(cell_types),dimnames = list(pathway_names,cell_types))
###calculate the pvalues using shuffle method
pvalues_mat = matrix(NA,nrow=length(pathway_names),ncol=length(cell_types),dimnames = (list(pathway_names, cell_types)))
norm_tpm = all_data

for(p in pathway_names){
  	genes = pathways[[p]]
  	genes_comm = intersect(genes, rownames(norm_tpm))
  	if(length(genes_comm) < 5) next

  	pathway_metabolic_tpm = norm_tpm[genes_comm, ]
  	pathway_metabolic_tpm = pathway_metabolic_tpm[rowSums(pathway_metabolic_tpm)>0,]
  
  	mean_exp_eachCellType = apply(pathway_metabolic_tpm, 1, function(x)by(x, all_cell_types, mean))

  	#remove genes which are zeros in any celltype to avoid extreme ratio value
  	keep = colnames(mean_exp_eachCellType)[colAlls(mean_exp_eachCellType>0.001)]

  	if(length(keep)<3) next
  
  	#using the loweset value to replace zeros for avoiding extreme ratio value
  	pathway_metabolic_tpm = pathway_metabolic_tpm[keep,]
  	pathway_metabolic_tpm = t( apply(pathway_metabolic_tpm,1,function(x) {x[x<=0] = min(x[x>0]);x} ))

  
  	pathway_number_weight = 1 / gene_pathway_number[keep,]
  	mean_exp_eachCellType = apply(pathway_metabolic_tpm, 1, function(x)by(x, all_cell_types, mean))
  	ratio_exp_eachCellType = t(mean_exp_eachCellType) / colMeans(mean_exp_eachCellType)
  	#exclude the extreme ratios
  	col_quantile = apply(ratio_exp_eachCellType,2,function(x) quantile(x,na.rm=T))
  	col_q1 = col_quantile["25%",]
  	col_q3 = col_quantile["75%",]
  	col_upper = col_q3 * 3
  	col_lower = col_q1 / 3
  	outliers = apply(ratio_exp_eachCellType,1,function(x) {any( (x>col_upper)|(x<col_lower) )} )
  
  	if(sum(!outliers) < 3) next
  
  	keep = names(outliers)[!outliers]
  	pathway_metabolic_tpm = pathway_metabolic_tpm[keep,]
  	pathway_number_weight = 1 / gene_pathway_number[keep,]
  	mean_exp_eachCellType = apply(pathway_metabolic_tpm, 1, function(x)by(x, all_cell_types, mean))
  	ratio_exp_eachCellType = t(mean_exp_eachCellType) / colMeans(mean_exp_eachCellType)
  	mean_exp_pathway = apply(ratio_exp_eachCellType,2, function(x) weighted.mean(x, pathway_number_weight/sum(pathway_number_weight)))
  	mean_expression_shuffle[p, ] =  mean_exp_pathway[cell_types]
  	mean_expression_noshuffle[p, ] =  mean_exp_pathway[cell_types]
    
  	##shuffle 5000 times:  
  	##define the functions 
  	group_mean = function(x){
    	sapply(cell_types,function(y) rowMeans(pathway_metabolic_tpm[,shuffle_cell_types_list[[x]]==y,drop=F])) }
  	column_weigth_mean = function(x){
    	apply(ratio_exp_eachCellType_list[[x]],2, function(y) weighted.mean(y, weight_values)) }
  	#####  
  	times = 1:5000
 	weight_values = pathway_number_weight/sum(pathway_number_weight)
  	shuffle_cell_types_list = lapply(times,function(x) sample(all_cell_types)) 
  	names(shuffle_cell_types_list) = times
  	mean_exp_eachCellType_list = lapply(times,function(x) group_mean(x))
  	ratio_exp_eachCellType_list = lapply(times,function(x) mean_exp_eachCellType_list[[x]] / rowMeans(mean_exp_eachCellType_list[[x]]))
  	mean_exp_pathway_list = lapply(times,function(x) column_weigth_mean(x))
  
  	shuffle_results = matrix(unlist(mean_exp_pathway_list),ncol=length(cell_types),byrow = T) 
  	rownames(shuffle_results) = times
  	colnames(shuffle_results) = cell_types
  	for(c in cell_types){
    	if(is.na(mean_expression_shuffle[p,c])) next
    	if(mean_expression_shuffle[p,c]>1){
      		pval = sum(shuffle_results[,c] > mean_expression_shuffle[p,c]) / 5000 
    	}else if(mean_expression_shuffle[p,c]<1){
      		pval = sum(shuffle_results[,c] < mean_expression_shuffle[p,c]) / 5000 }
    	if(pval>0.01) mean_expression_shuffle[p, c] = NA  ### NA is  blank in heatmap
    	pvalues_mat[p,c] = pval
  	}
}
all_NA = rowAlls(is.na(mean_expression_shuffle))
mean_expression_shuffle = mean_expression_shuffle[!all_NA,]
dat = mean_expression_shuffle
sort_row = c()
sort_column = c()

for(i in colnames(dat)){
  	select_row = which(rowMaxs(dat,na.rm = T) == dat[,i])
  	tmp = rownames(dat)[select_row][order(dat[select_row,i],decreasing = T)]
  	sort_row = c(sort_row,tmp) }
sort_column = apply(dat[sort_row,],2,function(x) order(x)[nrow(dat)])
sort_column = names(sort_column)
dat[is.na(dat)] = 1
pdf("KEGGpathway_activity_heatmap.pdf",onefile=T,width=6,height=9)
mybreaks = c( seq(0.8, 0.85, length.out=23), seq(0.86, 1.05, length.out=43), seq(1.06, max(dat),length.out=34)) 
color = colorRampPalette(c("blue","white","red"))(100)
pheatmap(dat[sort_row,sort_column],cluster_cols = F,cluster_rows = F,
color=color,breaks=mybreaks,fontsize_row=6,fontsize_col=6,cellwidth=9,cellheight=6)
dev.off()

write.table(mean_expression_noshuffle,file="KEGGpathway_activity_noshuffle.txt",row.names=T,col.names=T,quote=F,sep="\t")
write.table(mean_expression_shuffle,file="KEGGpathway_activity_shuffle.txt",row.names=T,col.names=T,quote=F,sep="\t")
write.table(pvalues_mat,file="KEGGpathway_activity_shuffle_pvalue.txt",row.names=T,col.names=T,quote=F,sep="\t")


datpath = '/public/workspace/AML_Leukemia/AML_Leukemia.integrated.rds'
tmp_data=readRDS(datpath)
info = data.frame()
a = names(table(tmp_data$seurat_clusters))
for(i in 1:length(a)){
    n = table(tmp_data$State2[tmp_data$seurat_clusters==a[i]])/sum(table(tmp_data$State2[tmp_data$seurat_clusters==a[i]]))
    label = names(which(n>=0.5))
    if(length(label)!=0){
        tmp = data.frame(cluster=a[i],label=label)
        info = rbind(info,tmp)
    } }

dat = read.table("KEGGpathway_activity_shuffle.txt",row.names=1,header=T,sep="\t",stringsAsFactors=F)
dat[is.na(dat)] = 1
cluster.case = paste("X",info$cluster[info$label=="QSC"],sep="")#9
cluster.control = paste("X",info$cluster[info$label=="PSP"],sep="")#11
cm=cbind(dat[,cluster.case],dat[,cluster.control])

p=apply(cm,1,function(x){wilcox.test(as.numeric(x[1:9]),as.numeric(x[10:20]))[[3]]})
filter1=names(which(p<=0.05))

fc=apply(cm,1,function(x){mean(as.numeric(x[1:9]))/mean(as.numeric(x[10:20]))})

cm=cm[filter1,]
colnames(cm)=gsub("X","C",colnames(cm))

pdf("Figure3_C.pdf",useDingbats=F)
annotation_col=data.frame(Type=c(rep("QSC",9),
    rep("PSP",11)))
rownames(annotation_col)=gsub("X","C",colnames(cm))
ann_colors = list(Type= c('QSC'="darkred",
    'PSP'='lightblue')) 
bk = c(seq(-2,2,by=0.02))
pheatmap(cm,
    breaks=bk,scale="row",
    color = c(colorRampPalette(colors = c("lightblue","white"))(length(bk)/2),
    		  colorRampPalette(colors = c("white","darkred"))(length(bk)/2)),
    cluster_cols=TRUE,cluster_rows=FALSE,border=NA,
    annotation=annotation_col,show_rownames=TRUE,border_color = "white",
    show_colnames=TRUE,annotation_colors=ann_colors,fontsize=7)
dev.off()



###==================================================================================
### Labeltransform for matched scRNA-seq datasets of Pt3
options(stringsAsFactors = F)
source('~/R/functions.R')
source('~/R/requiredPcg.R')

predicted.cols = c(`QSC`="#10507F",`PSP`="#ED8C25",
		`GMP`="#582566",`PG`="#E71F19",`ProMono`="#FFD92F",`Mono`="#6BB82D",
		Myeloid="olivedrab2",HSC_Prog="#4169E1",T_NK="lightgrey",B="#737373",Ery="#FADBDF")
tree = readRDS('/public/workspace/AML_Leukemia/AML_HC_labeltransform_tree.rds')
dat = readRDS('/public/workspace/AML_Leukemia/Pt3/integrated2Sample.rds')
anchors = FindTransferAnchors(reference = tree, query = dat, reference.assay = 'RNA', features = VariableFeatures(object = tree), query.assay = "RNA",  reduction = "cca")
predictions = TransferData(anchorset = anchors, refdata = tree$State2, dims = 1:30, weight.reduction="cca") 
predictions = data.frame(predictions,source=dat$source)
save(predictions,file = 'Pt3.labeltransform.RData')



###==================================================================================
### Labeltransform for matched scRNA-seq datasets of Pt9
options(stringsAsFactors = F)
source('~/R/functions.R')
source('~/R/requiredPcg.R')

predicted.cols = c(`QSC`="#10507F",`PSP`="#ED8C25",
		`GMP`="#582566",`PG`="#E71F19",`ProMono`="#FFD92F",`Mono`="#6BB82D",
		Myeloid="olivedrab2",HSC_Prog="#4169E1",T_NK="lightgrey",B="#737373",Ery="#FADBDF")
tree = readRDS('/public/workspace/AML_Leukemia/AML_HC_labeltransform_tree.rds')
dat = readRDS('/public/workspace/AML_Leukemia/Pt9/integrated2Sample.rds')
anchors = FindTransferAnchors(reference = tree, query = dat, reference.assay = 'RNA', features = VariableFeatures(object = tree), query.assay = "RNA",  reduction = "cca")
predictions = TransferData(anchorset = anchors, refdata = tree$State2, dims = 1:30, weight.reduction="cca") 
predictions = data.frame(predictions,source=dat$source)
save(predictions,file = 'Pt9.labeltransform.RData')




###==================================================================================
### Labeltransform for matched scRNA-seq datasets of Pt10
options(stringsAsFactors = F)
source('~/R/functions.R')
source('~/R/requiredPcg.R')

predicted.cols = c(`QSC`="#10507F",`PSP`="#ED8C25",
		`GMP`="#582566",`PG`="#E71F19",`ProMono`="#FFD92F",`Mono`="#6BB82D",
		Myeloid="olivedrab2",HSC_Prog="#4169E1",T_NK="lightgrey",B="#737373",Ery="#FADBDF")
tree = readRDS('/public/workspace/AML_Leukemia/AML_HC_labeltransform_tree.rds')
dat = readRDS('/public/workspace/AML_Leukemia/Pt10/integrated2Sample.rds')
anchors = FindTransferAnchors(reference = tree, query = dat, reference.assay = 'RNA', features = VariableFeatures(object = tree), query.assay = "RNA",  reduction = "cca")
predictions = TransferData(anchorset = anchors, refdata = tree$State2, dims = 1:30, weight.reduction="cca") 
predictions = data.frame(predictions,source=dat$source)
save(predictions,file = 'Pt10.labeltransform.RData')



###=============================================================================================
### DEGs of each group of leukemia-like cells in Pt3
dat = readRDS('/public/workspace/AML_Leukemia/Pt3/integrated2Sample.rds')
load('Pt3.labeltransform.RData')
dat = AddMetaData(dat,predictions)
dat.tmp = subset(dat,subset=(Leukemia_Normal=='Leukemia_Like'))
DefaultAssay(dat.tmp) = 'RNA'

Idents(dat.tmp) = factor(dat.tmp$predicted.id2, levels = unique(dat.tmp$predicted.id2))
sapply(as.character(unique(dat.tmp$predicted.id2)),function(x){
	tmp = FindMarkers(dat.tmp, ident.1 = "NT", group.by = 'source', subset.ident = x,min.pct = 0.25, logfc.threshold = 0.25)
	write.table(tmp,paste0('Pt3_',x,'.txt'),row.names=T,col.names=T,sep='\t',quote=F)})

protein = read.table('/public/workspace/AML_Leukemia/protein_coding_gene.txt')
list = lapply(Sys.glob('*.txt'),function(x){
		tmp = read.table(x,header=T,row.names=1,sep='\t')
		tmp = tmp[intersect(setdiff(rownames(tmp),rownames(tmp)[grep('^RPL|^RPS|^MT',rownames(tmp))]),protein$V1),]
		tmp = tmp[which(tmp$pct.1-tmp$pct.2>=0.2),]
		tmp = tmp[order(tmp$avg_logFC,decreasing=T),]
		if(nrow(tmp)>15) tmp = tmp[1:15,]
		return(tmp)})
names(list) = gsub('.txt','',Sys.glob('*.txt'))



###=============================================================================================
### DEGs of each group of leukemia-like cells in Pt9
dat = readRDS('/public/workspace/AML_Leukemia/Pt9/integrated2Sample.rds')
load('Pt9.labeltransform.RData')
dat = AddMetaData(dat,predictions)
dat.tmp = subset(dat,subset=(Leukemia_Normal=='Leukemia_Like'))
DefaultAssay(dat.tmp) = 'RNA'

Idents(dat.tmp) = factor(dat.tmp$predicted.id2, levels = unique(dat.tmp$predicted.id2))
sapply(as.character(unique(dat.tmp$predicted.id2)),function(x){
	tmp = FindMarkers(dat.tmp, ident.1 = "NT", group.by = 'source', subset.ident = x,min.pct = 0.25, logfc.threshold = 0.25)
	write.table(tmp,paste0(gsub(':','.',x),'.txt'),row.names=T,col.names=T,sep='\t',quote=F)})

protein = read.table('/public/workspace/AML_Leukemia/protein_coding_gene.txt')
list = lapply(Sys.glob('*.txt'),function(x){
		tmp = read.table(x,header=T,row.names=1,sep='\t')
		tmp = tmp[intersect(setdiff(rownames(tmp),rownames(tmp)[grep('^RPL|^RPS|^MT',rownames(tmp))]),protein$V1),]
		tmp = tmp[order(tmp$avg_logFC,decreasing=T),]
		if(nrow(tmp)>20) tmp = tmp[1:20,]
		return(tmp)})
names(list) = gsub('.txt','',Sys.glob('*.txt'))


