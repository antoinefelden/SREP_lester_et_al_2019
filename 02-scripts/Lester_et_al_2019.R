##### Different bacterial and viral pathogens trigger distinct immune responses in a globally invasive ant
##################################################################################################################################

##### Working directories ########################################################################################################

rm(list = ls())
DataDir = "/Users/antoinefelden/Documents/Research/Manuscripts/5-ViralLoads/Lester_et_al_2019/01-data"
FigDir = "/Users/antoinefelden/Documents/Research/Manuscripts/5-ViralLoads/Lester_et_al_2019/03-figures"

##### Libraries and functions ####################################################################################################

library("limma")
library("edgeR")
library("RColorBrewer")
library("gplots")
library("rtracklayer")
library("ggplot2")
library("gridExtra")
library("reshape")
library("multcomp")
library("mixOmics")
library("WGCNA")
library("pgirmess")
library("car")
library("PerformanceAnalytics")

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

my_palette <- colorRampPalette(c("blue", "white", "red3"))(n = 24)

##### RNA-seq data ###############################################################################################################

### Load data

gene_set = read.csv(paste(DataDir,"/Lhum_RefSeq.csv",sep=""))

pheno_data <- read.csv(paste(DataDir,"/phenotypes_northland.csv",sep=""))
data_rnaseq = read.csv(paste(DataDir,"/gene_count_matrix.csv",sep=""))

ref_immune <- read.csv(paste(DataDir,"/File_S1_lhum_immune_genes.csv",sep=""),na=TRUE,h=T); ref_immune <- unique(ref_immune[,c("LOC","Name","Pathway","Type","Comment")])
ref_candidates <- ref_immune
ref_candidates$ref_gene_name <- paste("LOC",ref_candidates$LOC,sep="")

# virus name database
virus_names <- data.frame("full_name" = c("Deformed wing virus", "Kashmir bee virus", "Linepithema humile bunyan-like virus 1", "Linepithema humile C virus 1", "Linepithema humile narna-like virus 1", "Linepithema humile partiti-like virus 1", "Linepithema humile picorna-like virus 1", "Linepithema humile polycipivirus 2", "Linepithema humile rhabdo-like virus 1", "Linepithema humile toti-like virus 1", "Linepithema humile virus 1"),
                          "point_name" = c("Deformed.wing.virus", "Kashmir.bee.virus", "Linepithema.humile.bunyan.like.virus.1", "Linepithema.humile.C.virus.1", "Linepithema.humile.narna.like.virus.1", "Linepithema.humile.partiti.like.virus.1", "Linepithema.humile.picorna.like.virus.1", "Linepithema.humile.polycipivirus.2", "Linepithema.humile.rhabdo.like.virus.1", "Linepithema.humile.toti.like.virus.1", "Linepithema.humile.virus.1"),
                          "symbol" =  c("DWV","KBV","LhuBLV1","LhuCV1","LhuNLV1","LhuPLV1","LhuPiLV1","LhuPcV2","LhuRLV1","LhuTLV1","LHUV-1"),
                          "family" = c("+","+","-","+","+","ds","+","+","-","ds","+"))

# creates DGElist object
x <- DGEList(as.matrix(data_rnaseq[,c(2:length(data_rnaseq))]),remove.zeros=F)

### Organising sample info

# Match the order of samples in pheno_data with the order of samples in data_rnaseq
pheno_data <- pheno_data[order(match(pheno_data$ids,colnames(data_rnaseq)[2:length(colnames(data_rnaseq))])),]

pheno_data$group <- ifelse(pheno_data$site == "WAI","bees","no_bees")

x$samples$group <- pheno_data$group
x$samples$site <- pheno_data$site

### Organising gene info

gtf <- rtracklayer::import(paste(DataDir,'/nth_CS323.gtf',sep=""))
gtf_df=as.data.frame(gtf)
gtf_df_transcripts <- droplevels(subset(gtf_df, type == "transcript"))

unique_gtf_df_transcripts <- unique(gtf_df_transcripts[,c("gene_id","ref_gene_name")])

# Calculate gene-level length (i.e. average from isoforms)
mean_length <- aggregate(Gene...Transcripts...Length~Gene...Symbol, data=unique(gene_set[,c("Gene...Symbol","Gene...Transcripts...Length")]),FUN="mean")
unique_description <- unique(gene_set[,c("Gene...Symbol","Gene...Description")])
unique_description_length <- merge(unique_description,mean_length,by="Gene...Symbol")

# Merge all infos
all_genes_info <- na.omit(merge(unique_gtf_df_transcripts,unique_description_length,by.x = "ref_gene_name", by.y = "Gene...Symbol",all.x=TRUE))
summary(all_genes_info)

x$genes <- data.frame("gene_id"=data_rnaseq[,1],"line"=seq(1:nrow(data_rnaseq)))

### Data pre-processing

cpm <- cpm(x)
log_cpm <- cpm(x,log=TRUE)

# Remove low-expressed transcripts
table(rowSums(x$counts==0)==ncol(x$counts))
keep.exprs <- rowSums(cpm>1) >= 3
x_filt <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)
dim(x_filt)

cpm_filt <- cpm(x_filt)
log_cpm_filt <- cpm(x_filt,log=TRUE)

samplenames <- colnames(x)

nsamples <- ncol(x)
colpal <- colorRampPalette(brewer.pal(8, "Set1"))
col <- colpal(nsamples)
par(mfrow=c(1,2))
plot(density(log_cpm[,1]), col=col[1], lwd=2, las=2,
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(log_cpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}

plot(density(log_cpm_filt[,1]), col=col[1], lwd=2, las=2,
          main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(log_cpm_filt[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}

dev.copy2pdf(file=paste(FigDir,"/filtering_low_expr.pdf",sep=""), width=8, height=6)


### Normalising gene expression distribution

x2 <- x_filt
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5

par(mfrow=c(1,2))
cpm_x2 <- cpm(x2,log=TRUE)
boxplot(cpm_x2, las=2, col=col, main="")
title(main="A. Unnormalised data",ylab="Log-cpm")
x2 <- calcNormFactors(x2)
x2$samples$norm.factors
cpm_x2 <- cpm(x2,log=TRUE)
boxplot(cpm_x2, las=2, col=col, main="")
title(main="B. Normalised data",ylab="Log-cpm")
dev.copy2pdf(file=paste(FigDir,"/normalisation.pdf",sep=""), width=8, height=6)

# Go ahead with TMM normalisation
x_filt_tmm = calcNormFactors(x_filt, method='TMM')
cpm_filt_tmm <- cpm(x_filt_tmm,normalized.lib.sizes=TRUE)
log_cpm_filt_tmm <- cpm(x_filt_tmm,log=TRUE)

### Unsupervised clustering of samples

group <- as.factor(pheno_data$group)

par(mfrow=c(1,2),mar=c(4,4,4,4))
col.group <- group
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
plotMDS(log_cpm_filt_tmm, labels=group, col=col.group, dim=c(1,2))
plotMDS(log_cpm_filt_tmm, labels=group, col=col.group, dim=c(2,3))

dev.copy2pdf(file=paste(FigDir,"/unsupervised_clustering.pdf",sep=""), width=8, height=6)

### Calculate FPKM values for immune genes of interest

ref_candidates_length <- merge(ref_candidates,all_genes_info,by.x="ref_gene_name",by.y="ref_gene_name")

cpm_filt_tmm <- data.frame(cpm(x_filt_tmm,normalized.lib.sizes=TRUE))

selected_ref_candidates <- subset(ref_candidates,ref_candidates$Pathway == "Imd"|ref_candidates$Pathway == "JAK-STAT"|ref_candidates$Pathway == "JNK"|ref_candidates$Pathway == "RNAi"|ref_candidates$Pathway == "TOLL")

# calculate fpkm values
cpm_filt_tmm$gene_id <- x_filt_tmm$genes$gene_id
candidates_cpm_filt_tmm <- merge(ref_candidates_length,cpm_filt_tmm,by.x="gene_id",by.y="gene_id")

rpkm_data <- rpkm(candidates_cpm_filt_tmm[,10:ncol(candidates_cpm_filt_tmm)],normalized.lib.sizes=FALSE,gene.length=candidates_cpm_filt_tmm$Gene...Transcripts...Length,log=TRUE)

t_data_num_all <- data.frame(scale(t(rpkm_data),scale=TRUE))
colnames(t_data_num_all) <- paste("LOC",candidates_cpm_filt_tmm$LOC,sep="")
t_data_num_candidates <- t_data_num_all[,na.omit(match(selected_ref_candidates$ref_gene_name,colnames(t_data_num_all)))]

# Preferred alternative: compute FPKM from x_filt_TMM on all genes, and then subset (not implemented in manuscript)
#x_filt_tmm = calcNormFactors(x_filt, method='TMM')
#x_filt_tmm$genes$length <- all_genes_info$Gene...Transcripts...Length[match(x_filt_tmm$genes$gene_id,all_genes_info$gene_id)]
#x_filt_tmm$genes$ref_gene_name <- all_genes_info[match(droplevels(x_filt_tmm$genes$gene_id),all_genes_info$gene_id),1]
#rpkm_all <- data.frame(rpkm(x_filt_tmm,normalized.lib.sizes=TRUE,gene.length=x_filt_tmm$genes$length,log=TRUE))
#t_data_num_all <- data.frame(scale(t(rpkm_all),scale=TRUE))
#colnames(t_data_num_all) <- x_filt_tmm$genes$ref_gene_name
#selected_ref_candidates <- subset(ref_candidates,ref_candidates$Pathway == "Imd"|ref_candidates$Pathway == "JAK-STAT"|ref_candidates$Pathway == "JNK"|ref_candidates$Pathway == "RNAi"|ref_candidates$Pathway == "TOLL")
#t_data_num_candidates <- t_data_num_all[,na.omit(match(selected_ref_candidates$ref_gene_name,colnames(t_data_num_all)))]

### Viral loads

BLASTx_all <- rbind(read.csv(paste(DataDir,"/viruses_nth_refseq.blastx.outfmt6_clean.csv",sep="")),
                            read.csv(paste(DataDir,"/viljakainen_viruses_nthland.blastx.outfmt6_clean.csv",sep="")))
BLASTx_all$sseqid <- gsub("[|]","",gsub(".*gb[|]","",BLASTx_all$sseqid))
                   
expression_matrix <- read.table(paste(DataDir,"/genes.TMM.EXPR.matrix",sep=""))                  
expression_matrix <- expression_matrix[-which(rowSums(expression_matrix)==0),]

# Use library sizes of reads mapped to the Argentine ant genome to normalise viral counts -> host-library-size-standardised TMM-normalised TPMs
mapped_reads_lib_size <- data.frame("sample" = rownames(x$samples), "alias" = colnames(expression_matrix), "lib.size" = x$samples$lib.size)
for(n in seq(1:ncol(expression_matrix))){
  expression_matrix[,n] <- (expression_matrix[,n]/mapped_reads_lib_size[n,3])*(mean(mapped_reads_lib_size[,3]))
}

expression_matrix$gene_id <- rownames(expression_matrix)

virus_hits <- read.csv(paste(DataDir,"/virus_hits_nthland.csv",sep=""))

list_TOI_BLASTx <- unique(BLASTx_all[,1:2])
TOI_expr_BLASTx <- merge(list_TOI_BLASTx,expression_matrix,by.x="qseqid",by.y="gene_id")[,-1]
TOI_expr_BLASTx <- aggregate(TOI_expr_BLASTx[,-1],by=list(sseqid=TOI_expr_BLASTx$sseqid),sum)

TOI_expr_BLASTx$virus <- virus_hits$organism[match(TOI_expr_BLASTx$sseqid,virus_hits$accession)]
rownames(TOI_expr_BLASTx) <- TOI_expr_BLASTx$sseqid; TOI_expr_BLASTx <- TOI_expr_BLASTx[,-1]

t_data_num_viruses <- data.frame(scale(t(TOI_expr_BLASTx[,-ncol(TOI_expr_BLASTx)]),center = FALSE,scale=TRUE))


data_microbe=data.frame("sample"= pheno_data$ids)
  for (microbe in as.character(unique(virus_names$full_name))) {
  print(microbe)
  subset_accessions <- subset(virus_hits,virus_hits$organism == microbe); subset_accessions <- as.character(subset_accessions$accession)
  subset_loads <- data.frame(t_data_num_viruses[,match(subset_accessions,colnames(t_data_num_viruses))])
  print(subset_loads)
  if (length(subset_loads)>1) {sum_loads <- data.frame(as.numeric(rowSums(subset_loads)))} else {sum_loads <- subset_loads}
  colnames(sum_loads) <- microbe
  print(sum_loads)
  data_microbe <- cbind(data_microbe,sum_loads)
  }
data_microbe

# immune genes
X <- as.matrix(t_data_num_candidates)#; colnames(X) <- paste("LOC",colnames(X),sep="")
pathway_genes <- data.frame("ref_gene_name" = colnames(X), "description" = ref_candidates[match(colnames(X),ref_candidates$ref_gene_name),2],
                            "pathway" = ref_candidates[match(colnames(X),ref_candidates$ref_gene_name),3])

# pathogens
Y <- data.frame(data_microbe[,2:length(data_microbe)]); rownames(Y) <- data_microbe[,1]
colnames(Y) <- virus_names[match(colnames(Y),virus_names$point_name),2]

# Compile data for PLS
PLS_data <- list("genes"=data.frame(X),"pathogens"=data.frame(Y),"phenotype"=pheno_data,"gene_names"=pathway_genes, "pathogen_names" = virus_names)
col_genes <- c("darkblue", "darkorange1", "purple","red3","dodgerblue")
col_pathogens <- c("darkgreen", "darkolivegreen1", "darkseagreen3","dimgrey")

# PLS
PLS <- pls(Y=PLS_data$genes,X=PLS_data$pathogens,ncomp=8, mode = "regression")

tune.PLS <- perf(PLS, validation = "Mfold", folds = 10, progressBar = FALSE, nrepeat = 100)
plot(tune.PLS$Q2.total,main="Selection of the number of components")
abline(h=0.0975,col="red")
dev.copy2pdf(file=paste(FigDir,"/perf_plot.pdf",sep=""), width=10, height=6)

par(mfrow=c(1,2))
plot_cim <- cim(PLS, comp=c(1,4),margins=c(7,38), row.cex = 1.2, col.cex = 1.5, transpose = TRUE,
                col.names=PLS_data$gene_names[,2], row.names=PLS_data$pathogen_names[,3],
                row.sideColors = col_pathogens[factor(PLS_data$pathogen_names[,4])],
                col.sideColors = col_genes[factor(PLS_data$gene_names[,3])])
legend_genes <- legend(x=0.31,y=0.3,legend=sort(unique(PLS_data$gene_names[,3])),fill=col_genes,ncol=1,bty="o",cex=1.5)
legend_viruses <- legend(x=0.31,y=0.6,legend=sort(unique(PLS_data$pathogen_names[,3])),fill=col_pathogens,ncol=1,bty="o",cex=1.5)
dev.copy2pdf(file=paste(FigDir,"/heatmap_PLS_selected_cp1_4.pdf",sep=""), width=16, height=20)

##### Correlations and principal components analysis ###################################################################################

### Data that isn't logged:
data.df <- read.csv(paste(DataDir,"/File_S2_Immune_qPCR.csv",sep=""),T)
head(data.df)

### The Box-Cox family of scaled power transformations for variable x equals (x^(lambda)-1)/lambda for lambda not equal to 0, and log(x) if lambda = 0. 
powerTransform(data.df+1)
summary(powerTransform(data.df+1))
fit_powers = powerTransform(data.df+1)$roundlam


trans_data = bcPower(data.df+1,fit_powers)

names(trans_data) = paste0(names(data.df), "_t")

#### Correlations 

# Plot of values (for Fig. 3 below diagonal)
chart.Correlation(trans_data, histogram = TRUE)


# Table of Spearman correlation values (for Fig. 3 above diagonal)
res <- cor(trans_data, method="spearman")
round(res, 2)

### Information for sequential Bonferroni corrections (for Fig. 3)
attach(trans_data)

y = I_J_HOP_t
cor.test(y,I_R_DCR_t, method="spearman")
cor.test(y,I_R_AGO_t, method="spearman")
cor.test(y,I_T_TOLL_t, method="spearman")
cor.test(y,I_T_DEF2_t, method="spearman")
cor.test(y,I_T_PELLE_t, method="spearman")
cor.test(y,I_I_FADD_t, method="spearman")
cor.test(y,I_I_JNK_t, method="spearman")
cor.test(y,I_I_IKB_t, method="spearman")
cor.test(y,V_LHUV_t, method="spearman")
cor.test(y,V_KBV_t, method="spearman")
cor.test(y,V_DWV_t, method="spearman")
cor.test(y,V_BQCV_t, method="spearman")
cor.test(y,B_FA_t, method="spearman")

y = I_R_DCR_t
cor.test(y,I_R_AGO_t, method="spearman")
cor.test(y,I_T_TOLL_t, method="spearman")
cor.test(y,I_T_DEF2_t, method="spearman")
cor.test(y,I_T_PELLE_t, method="spearman")
cor.test(y,I_I_FADD_t, method="spearman")
cor.test(y,I_I_JNK_t, method="spearman")
cor.test(y,I_I_IKB_t, method="spearman")
cor.test(y,V_LHUV_t, method="spearman")
cor.test(y,V_KBV_t, method="spearman")
cor.test(y,V_DWV_t, method="spearman")
cor.test(y,V_BQCV_t, method="spearman")
cor.test(y,B_FA_t, method="spearman")

y = I_R_AGO_t
cor.test(y,I_T_TOLL_t, method="spearman")
cor.test(y,I_T_DEF2_t, method="spearman")
cor.test(y,I_T_PELLE_t, method="spearman")
cor.test(y,I_I_FADD_t, method="spearman")
cor.test(y,I_I_JNK_t, method="spearman")
cor.test(y,I_I_IKB_t, method="spearman")
cor.test(y,V_LHUV_t, method="spearman")
cor.test(y,V_KBV_t, method="spearman")
cor.test(y,V_DWV_t, method="spearman")
cor.test(y,V_BQCV_t, method="spearman")
cor.test(y,B_FA_t, method="spearman")

y = I_T_TOLL_t
cor.test(y,I_T_DEF2_t, method="spearman")
cor.test(y,I_T_PELLE_t, method="spearman")
cor.test(y,I_I_FADD_t, method="spearman")
cor.test(y,I_I_JNK_t, method="spearman")
cor.test(y,I_I_IKB_t, method="spearman")
cor.test(y,V_LHUV_t, method="spearman")
cor.test(y,V_KBV_t, method="spearman")
cor.test(y,V_DWV_t, method="spearman")
cor.test(y,V_BQCV_t, method="spearman")
cor.test(y,B_FA_t, method="spearman")

y = I_T_DEF2_t
cor.test(y,I_T_PELLE_t, method="spearman")
cor.test(y,I_I_FADD_t, method="spearman")
cor.test(y,I_I_JNK_t, method="spearman")
cor.test(y,I_I_IKB_t, method="spearman")
cor.test(y,V_LHUV_t, method="spearman")
cor.test(y,V_KBV_t, method="spearman")
cor.test(y,V_DWV_t, method="spearman")
cor.test(y,V_BQCV_t, method="spearman")
cor.test(y,B_FA_t, method="spearman")

y = I_T_PELLE_t
cor.test(y,I_I_FADD_t, method="spearman")
cor.test(y,I_I_JNK_t, method="spearman")
cor.test(y,I_I_IKB_t, method="spearman")
cor.test(y,V_LHUV_t, method="spearman")
cor.test(y,V_KBV_t, method="spearman")
cor.test(y,V_DWV_t, method="spearman")
cor.test(y,V_BQCV_t, method="spearman")
cor.test(y,B_FA_t, method="spearman")

y = I_I_FADD_t
cor.test(y,I_I_JNK_t, method="spearman")
cor.test(y,I_I_IKB_t, method="spearman")
cor.test(y,V_LHUV_t, method="spearman")
cor.test(y,V_KBV_t, method="spearman")
cor.test(y,V_DWV_t, method="spearman")
cor.test(y,V_BQCV_t, method="spearman")
cor.test(y,B_FA_t, method="spearman")

y = I_I_JNK_t
cor.test(y,I_I_IKB_t, method="spearman")
cor.test(y,V_LHUV_t, method="spearman")
cor.test(y,V_KBV_t, method="spearman")
cor.test(y,V_DWV_t, method="spearman")
cor.test(y,V_BQCV_t, method="spearman")
cor.test(y,B_FA_t, method="spearman")

y = I_I_IKB_t
cor.test(y,V_LHUV_t, method="spearman")
cor.test(y,V_KBV_t, method="spearman")
cor.test(y,V_DWV_t, method="spearman")
cor.test(y,V_BQCV_t, method="spearman")
cor.test(y,B_FA_t, method="spearman")

y = V_LHUV_t
cor.test(y,V_KBV_t, method="spearman")
cor.test(y,V_DWV_t, method="spearman")
cor.test(y,V_BQCV_t, method="spearman")
cor.test(y,B_FA_t, method="spearman")

y = V_KBV_t
cor.test(y,V_DWV_t, method="spearman")
cor.test(y,V_BQCV_t, method="spearman")
cor.test(y,B_FA_t, method="spearman")

y = V_DWV_t
cor.test(y,V_BQCV_t, method="spearman")
cor.test(y,B_FA_t, method="spearman")

y = V_BQCV_t
cor.test(y,B_FA_t, method="spearman")


### PCA ##########################################################################################################################

pc_pathogen = prcomp(trans_data[,10:14])
## no scaling - best in this situation, since all expressions are relative to the same reference genes


### Table 2 - with reordering of components
pc_pathogen
summary(pc_pathogen)


#######################
### Linear regression of each immune gene on all five pathogen PCs

### Shows how all 9 immume genes load on the 5 PCs. 

### Gives statistical significance of effect of each PC, 
### which shows which pathogens are most important for each immune gene

### Definitely have no dependence between explanatory vars - all PCs linearly independent by construction

#######################

imm_reg_mod = lm(as.matrix(trans_data[,1:9])~pc_pathogen$x)

### Table 3 info
summary(imm_reg_mod)
