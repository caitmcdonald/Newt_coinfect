
#######################################################################################################################

### Convert transcript-level abundance to gene counts using tximport, then continue with edgeR differential gene expression analysis ###

### README: It is not appropriate to just sum transcript counts to gene level because that violates a negative binomial assumption of DGE programs. Need to filter through tximport first #####
### README 2: performing DGE at gene level instead of transcript level is more biologically relevant #####

#######################################################################################################################



#######################################################################################################################
### RUN TXIMPORT TO GENERATE GENE COUNTS FROM SALMON TRANSCRIPT COUNTS ###
#######################################################################################################################

# load packages, setwd
library(tximport)
library(edgeR)
library(dplyr)
library(data.table)
setwd("~/nvir/nvir_tximport_edger/")

# create named vector pointing to quantification files:
## 1. create vector of filenames by reading in a table containing sample IDs
# should have row for each library, one column must match the names of the directories in which the quant files are held.
samples <- read.table("nvir_samples_file.txt", header=TRUE)

## 2. combine this with file.path and "quant.sf"
files <- file.path("~/nvir/nvir_tximport_edger", samples$rep, "quant.sf")
names(files) <- paste0(samples$rep)
all(file.exists(files))

# Read in a transcript-to-gene mapping file. (Can get this from various model org packages, or create your own. Must have first column=transcript ID, second column=gene ID)
tx2gene <- read.table("tx2gene.txt", header=TRUE)

# import transcript-level estimates, convert to gene counts
txi <- tximport(files, type="salmon", tx2gene = tx2gene)
#head(txi$length) #this is 'effective length' output from salmon (which itself is computed by taking into account factors that affect probability of sampling fragments from the transcript, e.g. fragment length distribution, sequence-specific bias, gc-bias)



#######################################################################################################################
### SET UP GENE COUNTS DATA FOR EDGER ###
#######################################################################################################################

# the gene-level estimate file above can now be used by edgeR

# identify gene counts
cts <- txi$counts[,c(1:4,13:16,25:28,37:40,5:8,17:20,29:32,41:44,9:12,21:24,33:36,45:48)] #re-order liver/skin/spleen
cts_data <- as.data.frame(cts, optional=TRUE)
cts_data <- tibble::rownames_to_column(cts_data, "gene")

# add annotation data to counts data
annotations <- fread("Newt_gene_swissprot_all_annot.txt") #for some reason, fread works while read.table doesn't!

nvir_counts <- left_join(cts_data,annotations)

# add effective gene length to annotations dataframe
eff_gene_length <- txi$length[1:49377]
nvir_counts1 <- cbind(nvir_counts, eff_gene_length)
nvir_counts1 <- nvir_counts1[,c(1,54,60,55,64,2:49)] #annotation data first
names(nvir_counts1)
#write.table(nvir_counts1, "~/nvir_final/Nvir_tximport_genes_and_annot_nofilter.txt", col.names=T, sep="\t", row.names=F, quote=F) #write raw counts e.g. to be used in volcano plots, heatmaps

# assign treatments for DGE experiment grouping
treatments <- factor(c(rep("liv.bd", 4), rep("liv.both", 4), rep("liv.bsal", 4), rep("liv.ctl", 4), rep("sk.bd", 4), rep("sk.both", 4), rep("sk.bsal", 4), rep("sk.ctl", 4), rep("spl.bd", 4), rep("spl.both", 4),  rep("spl.bsal", 4), rep("spl.ctl", 4)))



#######################################################################################################################
### FILTER, NORMALIZE AS SUGGESTED BY TXIMPORT VIGNETTE ###
#######################################################################################################################

# normalize by gene length
normMat <- txi$length[,c(1:4,13:16,25:28,37:40,5:8,17:20,29:32,41:44,9:12,21:24,33:36,45:48)] #reorder
normMat <- normMat/exp(rowMeans(log(normMat))) #normalize by gene length


# calculate offset, separated by tissue
o <- log(calcNormFactors(nvir_counts1[,c(6,7,9:23,25,27:39,41:48,50:52)]/normMat[,c(1,2,4:18,20,22:34,36:43,45:47)])) + log(colSums(nvir_counts1[,c(6,7,9:23,25,27:39,41:48,50:52)]/normMat[,c(1,2,4:18,20,22:34,36:43,45:47)])) #all samples, all tissues included
o.liv <- log(calcNormFactors(nvir_counts1[,c(6,7,9:21)]/normMat[c(1,2,4:16)])) + log(colSums(nvir_counts1[,c(6,7,9:21)]/normMat[c(1,2,4:16)])) #remove outlier bd_liv_rep3
o.skin <- log(calcNormFactors(nvir_counts1[,c(22,23,25,27:37)]/normMat[c(17,18,20,22:32)])) + log(colSums(nvir_counts1[,c(22,23,25,27:37)])/normMat[c(17,18,20,22:32)]) #remove bd_skin_rep3 and both_skin_rep1
o.spleen <- log(calcNormFactors(nvir_counts1[,c(38,39,41:48,50:52)]/normMat[c(33,34,36:43,45:47)])) + log(colSums(nvir_counts1[,c(38,39,41:48,50:52)]/normMat[c(33,34,36:43,45:47)])) #remove bd_spleen_rep3, bsal_spleen_rep4, control_spleen_rep4

# generate DGEList
y <- DGEList(counts=nvir_counts1[,c(6,7,9:23,25,27:39,41:48,50:52)], group=treatments[c(1,2,4:18,20,22:34,36:43,45:47)], genes=nvir_counts1[,1:4])
y.liv <- DGEList(counts=nvir_counts1[,c(6,7,9:21)], group=treatments[c(1,2,4:16)], genes=nvir_counts1[,1:4])
y.skin <- DGEList(counts=nvir_counts1[,c(22,23,25,27:37)], group=treatments[c(17,18,20,22:32)], genes=nvir_counts1[,1:4])
y.spleen <- DGEList(counts=nvir_counts1[,c(38,39,41:48,50:52)], group=treatments[c(33,34,36:43,45:47)], genes=nvir_counts1[,1:4])

# scale libraries by offset
y <- scaleOffset(y, t(t(log(normMat[,c(1,2,4:18,20,22:34,36:43,45:47)])) + o))
y.liv <- scaleOffset(y.liv, t(t(log(normMat[,c(1,2,4:16)])) + o.liv))
y.skin <- scaleOffset(y.skin, t(t(log(normMat[,c(17,18,20,22:32)])) + o.skin))
y.spleen <- scaleOffset(y.spleen, t(t(log(normMat[,c(33,34,36:43,45:47)])) + o.spleen))

# filter counts
keep <- rowSums(cpm(y)>2)>=3
keep.liv <- rowSums(cpm(y.liv)>2)>=3 #tximport recommends filterByExpr(y), but this is too relaxed imo and results in BCV>0.55
keep.skin <- rowSums(cpm(y.skin)>2)>=3
keep.spleen <- rowSums(cpm(y.spleen)>2)>=3

ykeep <- y[keep, , keep.lib.sizes=FALSE]
ykeep.liv <- y.liv[keep.liv, , keep.lib.sizes=FALSE]
ykeep.skin <- y.skin[keep.skin, , keep.lib.sizes=FALSE]
ykeep.spleen <- y.spleen[keep.spleen, , keep.lib.sizes=FALSE]

dim(ykeep)
#dim(ykeep.liv)
#dim(ykeep.skin) 
#dim(ykeep.spleen) 

allMDS <- plotMDS(ykeep)
#livMDS <- plotMDS(ykeep.liv)
#skinMDS <- plotMDS(ykeep.skin)
#spleenMDS <- plotMDS(ykeep.spleen)

###calculate RPKM for boxplots ----> need to have a dataframe with gene length included.
rpkm_l <- rpkm(ykeep.liv)
rpkm_sk <- rpkm(ykeep.skin)
rpkm_spl <- rpkm(ykeep.spleen)

genes_l <- ykeep.liv$genes
genes_sk <- ykeep.skin$genes
genes_spl <- ykeep.spleen$genes

rpkm_liv <- cbind(genes_l, rpkm_l)
rpkm_skin <- cbind(genes_sk, rpkm_sk)
rpkm_spleen <- cbind(genes_spl, rpkm_spl)

#write.table(rpkm_liv, "~/nvir_final/plots/boxplots/rpkm_liver.txt", sep="\t", quote=F, row.names=F)
#write.table(rpkm_skin, "~/nvir_final/plots/boxplots/rpkm_skin.txt", sep="\t", quote=F, row.names=F)
#write.table(rpkm_spleen, "~/nvir_final/plots/boxplots/rpkm_spleen.txt", sep="\t", quote=F, row.names=F)

# generate filtered and log-transformed counts matrix
#filt_trans <- log2(ykeep$counts+0.25)
#nvir_filt <- cbind(ykeep$genes,filt_trans)
#write.table(nvir_filt, "~/nvir_final/Nvir_counts_filtered_and_log.txt", sep="\t", quote=F)



#######################################################################################################################
### DESIGN MODEL MATRIX & ESTIMATE DISPERSION ###
#######################################################################################################################

# design models
design.mat.liv <- model.matrix(~ 0 + ykeep.liv$samples$group) #specify no intercept so pairwise compare all groups
colnames(design.mat.liv) <- levels(ykeep.liv$samples$group)

design.mat.skin <- model.matrix(~ 0 + ykeep.skin$samples$group)
colnames(design.mat.skin) <- levels(ykeep.skin$samples$group)

design.mat.spleen <- model.matrix(~ 0 + ykeep.spleen$samples$group)
colnames(design.mat.spleen) <- levels(ykeep.spleen$samples$group)

# estimate dispersions
ydisp.liv <- estimateGLMCommonDisp(ykeep.liv,design.mat.liv, verbose=TRUE)
ydisp.liv <- estimateGLMTrendedDisp(ydisp.liv,design.mat.liv, method="power") #model is fit with trended dispersion
ydisp.liv <- estimateGLMTagwiseDisp(ydisp.liv,design.mat.liv) #tagwise shown for plotting purposes
#plotBCV(ydisp.liv)

ydisp.skin <- estimateGLMCommonDisp(ykeep.skin,design.mat.skin, verbose=TRUE)
ydisp.skin <- estimateGLMTrendedDisp(ydisp.skin,design.mat.skin, method="power")
ydisp.skin <- estimateGLMTagwiseDisp(ydisp.skin,design.mat.skin)
#plotBCV(ydisp.skin)

ydisp.spleen <- estimateGLMCommonDisp(ykeep.spleen,design.mat.spleen, verbose=TRUE)
ydisp.spleen <- estimateGLMTrendedDisp(ydisp.spleen,design.mat.spleen, method="power")
ydisp.spleen <- estimateGLMTagwiseDisp(ydisp.spleen,design.mat.spleen)
#plotBCV(ydisp.spleen)

# fit model
fit.liv <- glmFit(ydisp.liv, design.mat.liv, robust=TRUE)
fit.skin <- glmFit(ydisp.skin, design.mat.skin, robust=TRUE)
fit.spleen <- glmFit(ydisp.spleen, design.mat.spleen, robust=TRUE)

# specify contrasts
my.contrasts.liv <- makeContrasts(  
  liv.bd.ctl=liv.bd-liv.ctl,
  liv.bsal.ctl=liv.bsal-liv.ctl,
  liv.both.ctl=liv.both-liv.ctl,
  liv.both.bd=(liv.both-liv.ctl)-(liv.bd-liv.ctl),
  liv.both.bsal=(liv.both-liv.ctl)-(liv.bsal-liv.ctl),
  liv.bd.bsal=(liv.bd-liv.ctl)-(liv.bsal-liv.ctl),
  liv.both.vs.single=(liv.both-liv.ctl)-((liv.bsal+liv.bd)/2-liv.ctl),
  liv.inf.ctl=(liv.bd+liv.bsal+liv.both)/3-liv.ctl,
  levels=design.mat.liv)

my.contrasts.skin <- makeContrasts(
  sk.bd.ctl=sk.bd-sk.ctl,
  sk.bsal.ctl=sk.bsal-sk.ctl,
  sk.both.ctl=sk.both-sk.ctl,
  sk.both.bd=(sk.both-sk.ctl)-(sk.bd-sk.ctl),
  sk.both.bsal=(sk.both-sk.ctl)-(sk.bsal-sk.ctl),
  sk.bd.bsal=(sk.bd-sk.ctl)-(sk.bsal-sk.ctl),
  sk.both.vs.single=(sk.both-sk.ctl)-((sk.bsal+sk.bd)/2-sk.ctl),
  sk.inf.ctl=(sk.bd+sk.bsal+sk.both)/3-sk.ctl,
  levels=design.mat.skin)

my.contrasts.spleen <- makeContrasts( 
  spl.bd.ctl=spl.bd-spl.ctl,
  spl.bsal.ctl=spl.bsal-spl.ctl,
  spl.both.ctl=spl.both-spl.ctl,
  spl.both.bd=(spl.both-spl.ctl)-(spl.bd-spl.ctl),
  spl.both.bsal=(spl.both-spl.ctl)-(spl.bsal-spl.ctl),
  spl.bd.bsal=(spl.bd-spl.ctl)-(spl.bsal-spl.ctl),
  spl.both.vs.single=(spl.both-spl.ctl)-((spl.bsal+spl.bd)/2-spl.ctl),
  spl.inf.ctl=(spl.bd+spl.bsal+spl.both)/3-spl.ctl,
  levels=design.mat.spleen)



#######################################################################################################################
### RUN GLM LRT FOR DGE ###
#######################################################################################################################

# liver
liv.inf.ctl <- glmLRT(fit.liv, contrast=my.contrasts.liv[,"liv.inf.ctl"])
liv.both.vs.single <- glmLRT(fit.liv, contrast=my.contrasts.liv[,"liv.both.vs.single"])
liv.both.bd <- glmLRT(fit.liv, contrast=my.contrasts.liv[,"liv.both.bd"])
liv.both.bsal <- glmLRT(fit.liv, contrast=my.contrasts.liv[,"liv.both.bsal"])
liv.bd.ctl <- glmLRT(fit.liv, contrast=my.contrasts.liv[,"liv.bd.ctl"])
liv.bsal.ctl <- glmLRT(fit.liv, contrast=my.contrasts.liv[,"liv.bsal.ctl"])
liv.both.ctl <- glmLRT(fit.liv, contrast=my.contrasts.liv[,"liv.both.ctl"])
liv.bd.bsal <- glmLRT(fit.liv, contrast=my.contrasts.liv[,"liv.bd.bsal"])
#summary(decideTests(liv.both.bsal, p.value = 0.05))


# skin
sk.inf.ctl <- glmLRT(fit.skin, contrast=my.contrasts.skin[,"sk.inf.ctl"])
sk.both.vs.single <- glmLRT(fit.skin, contrast=my.contrasts.skin[,"sk.both.vs.single"])
sk.both.bd <- glmLRT(fit.skin, contrast=my.contrasts.skin[,"sk.both.bd"])
sk.both.bsal <- glmLRT(fit.skin, contrast=my.contrasts.skin[,"sk.both.bsal"])
sk.bd.ctl <- glmLRT(fit.skin, contrast=my.contrasts.skin[,"sk.bd.ctl"])
sk.bsal.ctl <- glmLRT(fit.skin, contrast=my.contrasts.skin[,"sk.bsal.ctl"])
sk.both.ctl <- glmLRT(fit.skin, contrast=my.contrasts.skin[,"sk.both.ctl"])
sk.bd.bsal <- glmLRT(fit.skin, contrast=my.contrasts.skin[,"sk.bd.bsal"])
#summary(decideTests(sk.both.bsal, p.value = 0.05))


# spleen
spl.inf.ctl <- glmLRT(fit.spleen, contrast=my.contrasts.spleen[,"spl.inf.ctl"])
spl.both.vs.single <- glmLRT(fit.spleen, contrast=my.contrasts.spleen[,"spl.both.vs.single"])
spl.both.bd <- glmLRT(fit.spleen, contrast=my.contrasts.spleen[,"spl.both.bd"])
spl.both.bsal <- glmLRT(fit.spleen, contrast=my.contrasts.spleen[,"spl.both.bsal"])
spl.bd.ctl <- glmLRT(fit.spleen, contrast=my.contrasts.spleen[,"spl.bd.ctl"])
spl.bsal.ctl <- glmLRT(fit.spleen, contrast=my.contrasts.spleen[,"spl.bsal.ctl"])
spl.both.ctl <- glmLRT(fit.spleen, contrast=my.contrasts.spleen[,"spl.both.ctl"])
spl.bd.bsal <- glmLRT(fit.spleen, contrast=my.contrasts.spleen[,"spl.bd.bsal"])
#summary(decideTests(spl.both.bsal, p.value = 0.05))



#######################################################################################################################
### MAKE MD PLOTS OF DGE ###
#######################################################################################################################

# skin
#pdf(file="~/nvir_final/DGE/MD_sk_all.pdf", 9,7)
par(mfrow=c(2,4),oma=c(0,0,2,0))
plotMD(sk.inf.ctl, legend=FALSE, main="Infected vs. Control", ylim=c(-12,12),adjust.method = "BH", p.value = 0.05, cex=0.6)
abline(h = c(-1, 1), col = "grey")
plotMD(sk.bd.ctl, legend=FALSE, main="Bd vs. Control", ylim=c(-12,12),adjust.method = "BH", p.value = 0.05, cex=0.6)
abline(h = c(-1, 1), col = "grey")
plotMD(sk.bsal.ctl, legend=FALSE, main="Bsal vs. Control", ylim=c(-12,12),adjust.method = "BH", p.value = 0.05, cex=0.6)
abline(h = c(-1, 1), col = "grey")
plotMD(sk.both.ctl, legend=FALSE, main="Coinfected vs. Control", ylim=c(-12,12),adjust.method = "BH", p.value = 0.05, cex=0.6)
abline(h = c(-1, 1), col = "grey")
plotMD(sk.both.vs.single, legend=FALSE, main="Coinfected vs. Singly infected", ylim=c(-12,12),adjust.method = "BH", p.value = 0.05, cex=0.6)
abline(h = c(-1, 1), col = "grey")
plotMD(sk.both.bd, legend=FALSE, main="Coinfected vs. Bd", ylim=c(-12,12),adjust.method = "BH", p.value = 0.05, cex=0.6)
abline(h = c(-1, 1), col = "grey")
plotMD(sk.both.bsal, legend=FALSE, main="Coinfected vs. Bsal", ylim=c(-12,12),adjust.method = "BH", p.value = 0.05, cex=0.6)
abline(h = c(-1, 1), col = "grey")
plotMD(sk.bd.bsal, legend=FALSE, main="Bd vs. Bsal", ylim=c(-12,12),adjust.method = "BH", p.value = 0.05, cex=0.6)
abline(h = c(-1, 1), col = "grey")
title("Differential expression: SKIN", outer=TRUE) 
#dev.off()

# spleen
#pdf(file="~/nvir_final/DGE/MD_spl_all.pdf", 9,7)
par(mfrow=c(2,4),oma=c(0,0,2,0))
plotMD(spl.inf.ctl, legend=FALSE, main="Infected vs. Control", ylim=c(-12,12), adjust.method = "BH", p.value = 0.05, cex=0.6)
abline(h = c(-1, 1), col = "grey")
plotMD(spl.bd.ctl, legend=FALSE, main="Bd vs. Control", ylim=c(-12,12), adjust.method = "BH", p.value = 0.05, cex=0.6)
abline(h = c(-1, 1), col = "grey")
plotMD(spl.bsal.ctl, legend=FALSE, main="Bsal vs. Control", ylim=c(-12,12), adjust.method = "BH", p.value = 0.05, cex=0.6)
abline(h = c(-1, 1), col = "grey")
plotMD(spl.both.ctl, legend=FALSE, main="Coinfected vs. Control", ylim=c(-12,12), adjust.method = "BH", p.value = 0.05, cex=0.6)
abline(h = c(-1, 1), col = "grey")
plotMD(spl.both.vs.single, legend=FALSE, main="Coinfected vs. Singly infected", ylim=c(-12,12), adjust.method = "BH", p.value = 0.05, cex=0.6)
abline(h = c(-1, 1), col = "grey")
plotMD(spl.both.bd, legend=FALSE, main="Coinfected vs. Bd", ylim=c(-12,12), adjust.method = "BH", p.value = 0.05, cex=0.6)
abline(h = c(-1, 1), col = "grey")
plotMD(spl.both.bsal, legend=FALSE, main="Coinfected vs. Bsal", ylim=c(-12,12), adjust.method = "BH", p.value = 0.05, cex=0.6)
abline(h = c(-1, 1), col = "grey")
plotMD(spl.bd.bsal, legend=FALSE, main="Bd vs. Bsal", ylim=c(-12,12), adjust.method = "BH", p.value = 0.05, cex=0.6)
abline(h = c(-1, 1), col = "grey")
title("Differential expression: SPLEEN", outer=TRUE) 
#dev.off()

# liver
#pdf(file="~/nvir_final/DGE/MD_liv_all.pdf", 9,7)
par(mfrow=c(2,4),oma=c(0,0,2,0))
plotMD(liv.inf.ctl, legend=FALSE, main="Infected vs. Control", ylim=c(-12,12), adjust.method = "BH", p.value = 0.05, cex=0.6)
abline(h = c(-1, 1), col = "grey")
plotMD(liv.bd.ctl, legend=FALSE, main="Bd vs. Control", ylim=c(-12,12),adjust.method = "BH", p.value = 0.05, cex=0.6)
abline(h = c(-1, 1), col = "grey")
plotMD(liv.bsal.ctl, legend=FALSE, main="Bsal vs. Control", ylim=c(-12,12),adjust.method = "BH", p.value = 0.05, cex=0.6)
abline(h = c(-1, 1), col = "grey")
plotMD(liv.both.ctl, legend=FALSE, main="Coinfected vs. Control", ylim=c(-12,12),adjust.method = "BH", p.value = 0.05, cex=0.6)
abline(h = c(-1, 1), col = "grey")
plotMD(liv.both.vs.single, legend=FALSE, main="Coinfected vs. Singly infected", ylim=c(-12,12),adjust.method = "BH", p.value = 0.05, cex=0.6)
abline(h = c(-1, 1), col = "grey")
plotMD(liv.both.bd, legend=FALSE, main="Coinfected vs. Bd", ylim=c(-12,12),adjust.method = "BH", p.value = 0.05, cex=0.6)
abline(h = c(-1, 1), col = "grey")
plotMD(liv.both.bsal, legend=FALSE, main="Coinfected vs. Bsal", ylim=c(-12,12),adjust.method = "BH", p.value = 0.05, cex=0.6)
abline(h = c(-1, 1), col = "grey")
plotMD(liv.bd.bsal, legend=FALSE, main="Bd vs. Bsal", ylim=c(-12,12),adjust.method = "BH", p.value = 0.05, cex=0.6)
abline(h = c(-1, 1), col = "grey")
title("Differential expression: LIVER", outer=TRUE) 
#dev.off()



#######################################################################################################################
### GENERATE ALL DGE TABLES OF ALL COMPARISONS ###
#######################################################################################################################

# specify test inputs
input <- list(sk.inf.ctl, sk.both.vs.single, sk.both.bd,sk.both.bsal,sk.bd.ctl,sk.bsal.ctl,sk.both.ctl,sk.bd.bsal)
#input <- list(spl.inf.ctl, spl.both.vs.single, spl.both.bd,spl.both.bsal,spl.bd.ctl,spl.bsal.ctl,spl.both.ctl,spl.bd.bsal)
#input <- list(liv.inf.ctl, liv.both.vs.single, liv.both.bd,liv.both.bsal,liv.bd.ctl,liv.bsal.ctl,liv.both.ctl,liv.bd.bsal)

# create empty list to hold output
topDE.volc <- list() #for volcano plots
topDE.down <- list()
topDE.up <- list()

# generate tables of DGE separated by increased or decreased expression
for (i in 1:length(input))
{
  topDE <- topTags(input[[i]], adjust.method="BH", n=NULL, p.value=0.05)
  topDE.volc[[i]] <- topDE[order(topDE$table$logFC),] #for volcano plots
  topDE.down[[i]] <- topDE[which(topDE$table$logFC <= 0),] #output is not ordered by logFC
  topDE.up[[i]] <- topDE[which(topDE$table$logFC>=0),]
}

# create file names vector
names(topDE.up) <- c("sk_inf_ctl_up", "sk_both_vs_single_up", "sk_both_bd_up","sk_both_bsal_up","sk_bd_ctl_up","sk_bsal_ctl_up","sk_both_ctl_up","sk_bd_bsal_up")
names(topDE.down) <- c("sk_inf_ctl_down", "sk_both_vs_single_down", "sk_both_bd_down","sk_both_bsal_down","sk_bd_ctl_down","sk_bsal_ctl_down","sk_both_ctl_down","sk_bd_bsal_down")

#names(topDE.up) <- c("spl_inf_ctl_up", "spl_both_vs_single_up", "spl_both_bd_up","spl_both_bsal_up","spl_bd_ctl_up","spl_bsal_ctl_up","spl_both_ctl_up","spl_bd_bsal_up")
#names(topDE.down) <- c("spl_inf_ctl_down", "spl_both_vs_single_down", "spl_both_bd_down","spl_both_bsal_down","spl_bd_ctl_down","spl_bsal_ctl_down","spl_both_ctl_down","spl_bd_bsal_down")

#names(topDE.up) <- c("liv_inf_ctl_up", "liv_both_vs_single_up", "liv_both_bd_up","liv_both_bsal_up","liv_bd_ctl_up","liv_bsal_ctl_up","liv_both_ctl_up","liv_bd_bsal_up")
#names(topDE.down) <- c("liv_inf_ctl_down", "liv_both_vs_single_down", "liv_both_bd_down","liv_both_bsal_down","liv_bd_ctl_down","liv_bsal_ctl_down","liv_both_ctl_down","liv_bd_bsal_down")

# create file names for volcano plots
names(topDE.volc) <- c("sk.inf.ctl", "sk.both.vs.single", "sk.both.bd","sk.both.bsal","sk.bd.ctl","sk.bsal.ctl","sk.both.ctl","sk.bd.bsal")
#names(topDE.volc) <- c("spl.inf.ctl", "spl.both.vs.single", "spl.both.bd","spl.both.bsal","spl.bd.ctl","spl.bsal.ctl","spl.both.ctl","spl.bd.bsal")
#names(topDE.volc) <- c("liv.inf.ctl", "liv.both.vs.single", "liv.both.bd","liv.both.bsal","liv.bd.ctl","liv.bsal.ctl","liv.both.ctl","liv.bd.bsal")

# write tables of all DEGs for volcano plots
setwd("~/nvir_final/plots/volcano/")
lapply(1:length(topDE.volc),
       function(i) write.table(topDE.volc[[i]],
                               file = paste0(names(topDE.volc[i]),".txt"),
                               row.names=F, quote=F, sep="\t"))

# write tables of all DEGs, separated up and down-regulated
setwd("~/nvir_final/DGE/")
lapply(1:length(topDE.down),
       function(i) write.table(topDE.down[[i]],
                               file = paste0(names(topDE.down[i]),".txt"),
                               row.names=F, quote=F, sep="\t"))
lapply(1:length(topDE.up),
       function(i) write.table(topDE.up[[i]],
                               file = paste0(names(topDE.up[i]),".txt"),
                               row.names=F, quote=F, sep="\t"))


#######################################################################################################################