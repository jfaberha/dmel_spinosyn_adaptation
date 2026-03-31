### In R ###
#configure r environment
setwd("/scratch/user/jfaberha/20260325_052109/admera/gp_analysis/rudflies_2023_redo/r")
library(emmeans)
library(matrixStats)
library(ACER)
library(poolSeq)
library(ggplot2)
library(gplots)
library(ggpubr)
library(stats)
library(ggfortify)
library(dplyr)
library(zoo)
library(rstatix)
library(slider)
library(RColorBrewer)
library(geiger)
library(RRHO)
library(RRHO2)
library(stringr)
library(scales)

## set treatment group color scheme upfront for plotting
#treatment group order: E, PA, SE, SP
sample_cols <- c("#D26183","#495184","#848556","#D9B851")

#################################
### Load required input files ###
#################################

## Variant Effect Predictor (VEP) annotation file with appended FLYCADD scores
vep <- read.table("filtered-all.annot.vcf.FLYCADD.tsv", header=TRUE) 
## Sample metadata table
haf.meta <- read.table("rudflies_2023_meta.tsv", header=TRUE)
## Hafpipe imputed allele frequency table
haf.freq <- read.delim("rudflies_2023_hafpipe.csv", header=TRUE, sep = ",")

## Do a bit of reformatting for the frequency table header
names(haf.freq) <- gsub("[.]af","",names(haf.freq))
names(haf.freq) <- gsub("X","",names(haf.freq))
names(haf.freq)[1]="CHROM"
names(haf.freq)[2]="POS"

## Filter for founder samples for PA, S, and E
haf.meta.filt <- haf.meta[haf.meta$batch == "a" & haf.meta$experiment == "spino" & (haf.meta$treat.fix == "PA" | haf.meta$treat.fix == "S" | haf.meta$treat.fix == "E"),]
## Use this next line if you want to run aCMH with balanced pops, since all pops represented in TPT3
#haf.meta.filt <- haf.meta.filt[which((haf.meta.filt$cage.fix %in% haf.meta.filt[haf.meta.filt$tpt==3,]$cage.fix)==TRUE),]
## Find matching sample names in the frequency table, post metadata filtering
haf.freq.filt <- haf.freq[, which((names(haf.freq) %in% haf.meta.filt$samp)==TRUE)]
## And run the reciprocal filter
haf.meta.filt <- haf.meta.filt[which((haf.meta.filt$samp %in% names(haf.freq.filt))==TRUE),]

## additional filter to remove low variance loci
haf.sites.filt <- na.omit(haf.freq[rowVars(as.matrix(haf.freq.filt))>0.001,c(1:2)])
haf.freq.filt <- na.omit(haf.freq.filt[rowVars(as.matrix(haf.freq.filt))>0.001,])

## tables need a bit of reformatting for some downstream plotting
haf.meta.filt$treat.fix <- as.factor(haf.meta.filt$treat.fix)
haf.meta.filt$tpt <- as.factor(haf.meta.filt$tpt)


######################################
### For downstream analysis create ###
### a TPT1-only metadata table     ###
######################################

## Create a new column, "condition", which contains the population outcomes for S samples
## For example; SE - extinct S populations, SP - persistent S populations
filt.table <- cbind(cage=c("3","7","11","15","21","27","33","37","41","45","2","6","10","14","20","24","31","36","40","44","1","5","9","13","19","23","29","35","39","43","1","5","9","13","19","23","29","35","39","43"),tpt=c("1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1"),condition=c("SP","SP","SE","SE","SE","SE","SP","SP","SE","SE","PA","PA","PA","PA","PA","PA","PA","PA","PA","PA","E","E","E","E","E","E","E","E","E","E","E","E","E","E","E","E","E","E","E","E"))
haf.meta.T1filt <- merge(haf.meta, filt.table, by=c("cage","tpt"))

haf.meta.T1filt <- unique(haf.meta.T1filt[haf.meta.T1filt$batch == "a" & (haf.meta.T1filt$treat == "S" | haf.meta.T1filt$treat == "PA" | haf.meta.T1filt$treat == "E"),])
haf.meta.T1filt <- haf.meta.T1filt[which((haf.meta.T1filt$samp %in% c(1:9,30:171))==TRUE),]
haf.freq.T1filt <- haf.freq[, which((names(haf.freq) %in% haf.meta.T1filt$samp)==TRUE)]

#additional filter to remove extreme low variance loci
haf.sites.T1filt <- haf.freq[rowVars(as.matrix(haf.freq.T1filt))>0.001,c(1:2)]
haf.freq.T1filt <- haf.freq.T1filt[rowVars(as.matrix(haf.freq.T1filt))>0.001,]

#tables need a bit of reformatting for glm to run
haf.meta.T1filt$treat <- as.factor(haf.meta.T1filt$treat)
haf.meta.T1filt$tpt <- as.factor(haf.meta.T1filt$tpt)
haf.meta.T1filt$condition <- as.factor(haf.meta.T1filt$condition)
haf.meta.T1filt <- haf.meta.T1filt[order(haf.meta.T1filt$samp),]

##########################################
### Create master table of GLM results ###
##########################################

### Start with PA and S contrasts ###
## Load PA vs S, all samples combined
contrast.PAvS.treat.all <- read.table("rudflies_2023_PAvS_treat.all.table.GLMcontrast.txt", header=TRUE)
contrast.PAvS.treat.all <- contrast.PAvS.treat.all[,-3] #remove AF mean column
names(contrast.PAvS.treat.all)[3] <- "PAvS" #label glm results column
## Load PA vs S for each timepoint
contrast.PAvS.treat <- read.table("rudflies_2023_PAvS_treat.GLMcontrast.txt", header=TRUE)
names(contrast.PAvS.treat) <- c("PAvS.T1","PAvS.T2","PAvS.T3","PAvS.T4")
contrast.PAvS.tpt <- read.table("rudflies_2023_PAvS_tpt.GLMcontrast.txt", header=TRUE)
names(contrast.PAvS.tpt) <- c("PAvS.T1vT2","PAvS.T1vT3","PAvS.T1vT4","PAvS.T2vT3","PAvS.T2vT4","PAvS.T3vT4") #label glm results columns
## Load PA and S "timepoint:treatment" interaction
contrast.PAvS.int <- read.table("rudflies_2023_PAvS_int.GLMcontrast.txt", header=TRUE)
names(contrast.PAvS.int) <- c("int.PAvS.T1vT2","int.PAvS.T1vT3","int.PAvS.T1vT4","int.PAvS.T2vT3","int.PAvS.T2vT4","int.PAvS.T3vT4") #label glm results columns
## Join all these glm results together, and include loci columns
haf.freq.PAvS.filt.loci <- read.table("rudflies_2023_PAvS_haf.freq.filt.loci.txt", header=TRUE) #load loci columns
contrast.PAvS <- cbind(haf.freq.PAvS.filt.loci, contrast.PAvS.treat, contrast.PAvS.tpt, contrast.PAvS.int)
contrast.PAvS <- merge(contrast.PAvS.treat.all, contrast.PAvS, by=c("CHROM","POS"))
contrast.PAvS <- contrast.PAvS[order(contrast.PAvS[,1], contrast.PAvS[,2]), ] #sort by locus
## Save all PA and S contrasts for posterity
#write.table(contrast.PAvS, file="rudflies_2023_redo.PAvS.multiGLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)

## Look at PA-only and S-only time-point contrast GLM results ##
# Load S-only results
contrast.tpt.S.table <- read.table("rudflies_2023_S_tpt.table.GLMcontrast.txt", header=TRUE)
names(contrast.tpt.S.table) <- c("CHROM","POS","S.af.mean","S.T1vT2","S.T1vT3","S.T1vT4","S.T2vT3","S.T2vT4","S.T3vT4") #label glm results columns
# Load PA-only results
contrast.tpt.PA.table <- read.table("rudflies_2023_PA_tpt.table.GLMcontrast.txt", header=TRUE)
names(contrast.tpt.PA.table) <- c("CHROM","POS","PA.af.mean","PA.T1vT2","PA.T1vT3","PA.T1vT4","PA.T2vT3","PA.T2vT4","PA.T3vT4") #label glm results columns
# Merge the two
contrast.tpt.S_only.PA_only <- merge(contrast.tpt.S.table, contrast.tpt.PA.table, by=c("CHROM","POS"))
contrast.tpt.S_only.PA_only <- contrast.tpt.S_only.PA_only[,-c(3,10)] #remove AF mean cols
# Save relabeled tables for posterity
#write.table(contrast.tpt.S.table, file="rudflies_2023_S_tpt.table.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)
#write.table(contrast.tpt.PA.table, file="rudflies_2023_PA_tpt.table.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)

## TPT1-only GLM contrasts: PA vs S vs SE vs E - all pairwise ##
contrast.PAvSvSEvE.table <- read.table("rudflies_2023_PAvSvSEvE.wLoci.GLMcontrast.txt", header=TRUE)
names(contrast.PAvSvSEvE.table) <- c("CHROM","POS","EvPA.T1","EvSE.T1","EvSP.T1","PAvSE.T1","PAvSP.T1","SEvSP.T1") #label glm results columns
# Save relabeled table for posterity
# write.table(contrast.PAvSvSEvE.table, file="rudflies_2023_PAvSvSEvE.wLoci.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)

## Load PA and E contrasts ##
# PA vs E, all time-points combined
contrastout.PAvE.table <- read.table("rudflies_2023_PAvE_treat.all.table.GLMcontrast.txt", header=TRUE)
# PA vs E, at individual time-points
contrastout.treat.PA.E.table <- read.table("rudflies_2023_PAvE_treat.table.GLMcontrast2.txt", header=TRUE)
# Combine these two tables
contrastout.PAvE.table <- merge(contrastout.PAvE.table, contrastout.treat.PA.E.table, by=c("CHROM","POS"))
names(contrastout.PAvE.table) <- c("CHROM","POS","PAvE","PAvE.T1","PAvE.T2","PAvE.T3","PAvE.T4") #label glm results columns
# Save relabeled table for posterity
#write.table(contrastout.PAvE.table, file="rudflies_2023_PAvE_treat.all.table.wLoci.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)

## Load SE and E contrasts ##
# S vs E, all time-points combined
contrastout.SvE.table <- read.table("rudflies_2023_SvE_treat.all.table.GLMcontrast.txt", header=TRUE)[,c(1,2,4)]
# S vs E, at individual time-points
contrastout.treat.S.E.table <- read.table("rudflies_2023_SvE_treat.GLMcontrast.table.txt", header=TRUE)
# Combine these two tables
contrastout.SvE.table <- merge(contrastout.SvE.table, contrastout.treat.S.E.table, by=c("CHROM","POS"))
# S and E: "treatment:time-point" interaction
contrastout.int.S.E.table <- read.table("rudflies_2023_SvE_int.GLMcontrast.table.txt", header=TRUE)
# Add these results
contrastout.SvE.table <- merge(contrastout.SvE.table, contrastout.int.S.E.table, by=c("CHROM","POS"))
names(contrastout.SvE.table) <- c("CHROM","POS","SvE","SvE.T1","SvE.T2","SvE.T3","SvE.T4","int.SvE.T1vT2","int.SvE.T1vT3","int.SvE.T1vT4","int.SvE.T2vT3","int.SvE.T2vT4","int.SvE.T3vT4") #label glm results columns
# Save full table for posterity
#write.table(contrastout.SvE.table, file="rudflies_2023_SvE_treat.all.table.wLoci.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)


### Merge all glm results ###
glm.all <- merge(contrast.PAvS, contrast.tpt.S_only.PA_only, by=c("CHROM","POS"))
glm.all <- merge(glm.all, contrast.PAvSvSEvE.table, by=c("CHROM","POS"))
glm.all <- merge(glm.all, contrastout.PAvE.table, by=c("CHROM","POS"))
glm.all <- merge(glm.all, contrastout.SvE.table, by=c("CHROM","POS"))
glm.all <- glm.all[glm.all$CHROM != "4",] #remove chromsome 4
glm.all <- glm.all %>%
  		arrange(CHROM,POS) #sort by locus

## How many sites remain after merging all filtered GLM results tables?
dim(glm.all)
#[1] 1298959      53

## P-value correction
for(i in c(3:53)) { #start after loci columns and append FDR cols at the end of the table
  glm.all <- cbind(glm.all, p.adjust(glm.all[,i], method = "fdr"))
  colnames(glm.all)[i+51] <- paste(names(glm.all)[i], ".fdr", sep="")
}

## Double check number of cols after adding FDR cols
dim(glm.all)
#[1] 1298959     104

## Calculate -log10 for p-values
for(i in 3:53) { #start after loci columns and append logp cols at the end of the table
  glm.all <- cbind(glm.all, as.numeric(-log10(glm.all[,i])))
  colnames(glm.all)[i+102] <- paste(names(glm.all)[i], ".logp", sep="")
}

## Only write table if needed, these files are huge!
#write.table(glm.all, file=paste("rudflies_2023_redo.glm.masterfile.txt", sep=""), quote = FALSE, row.names = F)

#################################
### Calculate rolling windows ###
#################################

### 20 consecutive SNP window ###
## Create indexes for each chromosome to slide over and calculate mean over 20-SNP windows
## If we don't do this, the average will be calculated across consecutive chromosomes

## 2L
glm_2L <- unique(glm.all[glm.all$CHROM=="2L",]) 
glm_2L <- glm_2L %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
glm_2L <- glm_2L[is.finite(rowSums(glm_2L[,-c(1:2)])),]
## create new empty file for filtered rows
glm_2L$POS <- as.integer(glm_2L$POS)
glm_2L.rolwin20 <- glm_2L[,c(1,2)]

for(i in 3:ncol(glm_2L)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_2L[is.finite(glm_2L[,i])=="TRUE",i])+100
	glm_2L[is.infinite(glm_2L[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_2L.rolwin20 <- cbind(glm_2L.rolwin20,rollapply(glm_2L[,i], width = 20, FUN = mean, align = "center", fill = NA))
    colnames(glm_2L.rolwin20)[i] <- paste(names(glm_2L)[i], ".rolwin20", sep="")
}

## 2R
glm_2R <- unique(glm.all[glm.all$CHROM=="2R",])
glm_2R <- glm_2R %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
glm_2R <- glm_2R[is.finite(rowSums(glm_2R[,-c(1:2)])),]
## create new empty file for filtered rows
glm_2R$POS <- as.integer(glm_2R$POS)
glm_2R.rolwin20 <- glm_2R[,c(1,2)]

for(i in 3:ncol(glm_2R)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_2R[is.finite(glm_2R[,i])=="TRUE",i])+100
	glm_2R[is.infinite(glm_2R[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_2R.rolwin20 <- cbind(glm_2R.rolwin20,rollapply(glm_2R[,i], width = 20, FUN = mean, align = "center", fill = NA))
    colnames(glm_2R.rolwin20)[i] <- paste(names(glm_2R)[i], ".rolwin20", sep="")
}

## 3L
glm_3L <- unique(glm.all[glm.all$CHROM=="3L",])
glm_3L <- glm_3L %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
glm_3L <- glm_3L[is.finite(rowSums(glm_3L[,-c(1:2)])),]
## create new empty file for filtered rows
glm_3L$POS <- as.integer(glm_3L$POS)
glm_3L.rolwin20 <- glm_3L[,c(1,2)]

for(i in 3:ncol(glm_3L)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_3L[is.finite(glm_3L[,i])=="TRUE",i])+100
	glm_3L[is.infinite(glm_3L[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_3L.rolwin20 <- cbind(glm_3L.rolwin20,rollapply(glm_3L[,i], width = 20, FUN = mean, align = "center", fill = NA))
    colnames(glm_3L.rolwin20)[i] <- paste(names(glm_3L)[i], ".rolwin20", sep="")
}

## 3R
glm_3R <- unique(glm.all[glm.all$CHROM=="3R",])
glm_3R <- glm_3R %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
glm_3R <- glm_3R[is.finite(rowSums(glm_3R[,-c(1:2)])),]
## create new empty file for filtered rows
glm_3R$POS <- as.integer(glm_3R$POS)
glm_3R.rolwin20 <- glm_3R[,c(1,2)]

for(i in 3:ncol(glm_3R)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_3R[is.finite(glm_3R[,i])=="TRUE",i])+100
	glm_3R[is.infinite(glm_3R[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_3R.rolwin20 <- cbind(glm_3R.rolwin20,rollapply(glm_3R[,i], width = 20, FUN = mean, align = "center", fill = NA))
    colnames(glm_3R.rolwin20)[i] <- paste(names(glm_3R)[i], ".rolwin20", sep="")
}

## X
glm_X <- unique(glm.all[glm.all$CHROM=="X",])
glm_X <- glm_X %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
glm_X <- glm_X[is.finite(rowSums(glm_X[,-c(1:2)])),]
## create new empty file for filtered rows
glm_X$POS <- as.integer(glm_X$POS) #reformat
glm_X.rolwin20 <- glm_X[,c(1,2)]

for(i in 3:ncol(glm_X)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_X[is.finite(glm_X[,i])=="TRUE",i])+100
	glm_X[is.infinite(glm_X[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_X.rolwin20 <- cbind(glm_X.rolwin20,rollapply(glm_X[,i], width = 20, FUN = mean, align = "center", fill = NA))
    colnames(glm_X.rolwin20)[i] <- paste(names(glm_X)[i], ".rolwin20", sep="")
}

## Now, join the the new rolling window tables for all chromosomes
glm.all.rolwin20 <- rbind(glm_2L.rolwin20,glm_2R.rolwin20,glm_3L.rolwin20,glm_3R.rolwin20,glm_X.rolwin20)

## Make a rolling window table with just p-value, remove logp and fdr columns 
glm.all.rolwin20.pval <- select(glm.all.rolwin20,-contains("logp"),-contains("fdr"))
## Make a rolling window table with just logp, selects them and append loci
glm.all.rolwin20.logp <- cbind(glm.all.rolwin20.pval[,c(1:2)],select(glm.all.rolwin20,contains("logp")))


## Create indexes for each chromosome to slide over and calculate SD over 20-SNP windows
## If we don't do this, the SD will be calculated across consecutive chromosomes

## 2L
glm_2L_sd <- unique(glm.all[glm.all$CHROM=="2L",]) 
glm_2L_sd <- glm_2L_sd %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
glm_2L_sd <- glm_2L_sd[is.finite(rowSums(glm_2L_sd[,-c(1:2)])),]
## create new empty file for filtered rows
glm_2L_sd$POS <- as.integer(glm_2L_sd$POS)
glm_2L_sd.rolwin20 <- glm_2L_sd[,c(1,2)]

for(i in 3:ncol(glm_2L_sd)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_2L_sd[is.finite(glm_2L_sd[,i])=="TRUE",i])+100
	glm_2L_sd[is.infinite(glm_2L_sd[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_2L_sd.rolwin20 <- cbind(glm_2L_sd.rolwin20,rollapply(glm_2L_sd[,i], width = 20, FUN = mean, align = "center", fill = NA))
    colnames(glm_2L_sd.rolwin20)[i] <- paste(names(glm_2L_sd)[i], ".rolwin20", sep="")
}

## 2R
glm_2R_sd <- unique(glm.all[glm.all$CHROM=="2R",])
glm_2R_sd <- glm_2R_sd %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
glm_2R_sd <- glm_2R_sd[is.finite(rowSums(glm_2R_sd[,-c(1:2)])),]
## create new empty file for filtered rows
glm_2R_sd$POS <- as.integer(glm_2R_sd$POS)
glm_2R_sd.rolwin20 <- glm_2R_sd[,c(1,2)]

for(i in 3:ncol(glm_2R_sd)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_2R_sd[is.finite(glm_2R_sd[,i])=="TRUE",i])+100
	glm_2R_sd[is.infinite(glm_2R_sd[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_2R_sd.rolwin20 <- cbind(glm_2R_sd.rolwin20,rollapply(glm_2R_sd[,i], width = 20, FUN = mean, align = "center", fill = NA))
    colnames(glm_2R_sd.rolwin20)[i] <- paste(names(glm_2R_sd)[i], ".rolwin20", sep="")
}

## 3L
glm_3L_sd <- unique(glm.all[glm.all$CHROM=="3L",])
glm_3L_sd <- glm_3L_sd %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
glm_3L_sd <- glm_3L_sd[is.finite(rowSums(glm_3L_sd[,-c(1:2)])),]
## create new empty file for filtered rows
glm_3L_sd$POS <- as.integer(glm_3L_sd$POS)
glm_3L_sd.rolwin20 <- glm_3L_sd[,c(1,2)]

for(i in 3:ncol(glm_3L_sd)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_3L_sd[is.finite(glm_3L_sd[,i])=="TRUE",i])+100
	glm_3L_sd[is.infinite(glm_3L_sd[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_3L_sd.rolwin20 <- cbind(glm_3L_sd.rolwin20,rollapply(glm_3L_sd[,i], width = 20, FUN = mean, align = "center", fill = NA))
    colnames(glm_3L_sd.rolwin20)[i] <- paste(names(glm_3L_sd)[i], ".rolwin20", sep="")
}

## 3R
glm_3R_sd <- unique(glm.all[glm.all$CHROM=="3R",])
glm_3R_sd <- glm_3R_sd %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
glm_3R_sd <- glm_3R_sd[is.finite(rowSums(glm_3R_sd[,-c(1:2)])),]
## create new empty file for filtered rows
glm_3R_sd$POS <- as.integer(glm_3R_sd$POS)
glm_3R_sd.rolwin20 <- glm_3R_sd[,c(1,2)]

for(i in 3:ncol(glm_3R_sd)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_3R_sd[is.finite(glm_3R_sd[,i])=="TRUE",i])+100
	glm_3R_sd[is.infinite(glm_3R_sd[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_3R_sd.rolwin20 <- cbind(glm_3R_sd.rolwin20,rollapply(glm_3R_sd[,i], width = 20, FUN = mean, align = "center", fill = NA))
    colnames(glm_3R_sd.rolwin20)[i] <- paste(names(glm_3R_sd)[i], ".rolwin20", sep="")
}

## X
glm_X_sd <- unique(glm.all[glm.all$CHROM=="X",])
glm_X_sd <- glm_X_sd %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
glm_X_sd <- glm_X_sd[is.finite(rowSums(glm_X_sd[,-c(1:2)])),]
## create new empty file for filtered rows
glm_X_sd$POS <- as.integer(glm_X_sd$POS) #reformat
glm_X_sd.rolwin20 <- glm_X_sd[,c(1,2)]

for(i in 3:ncol(glm_X_sd)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_X_sd[is.finite(glm_X_sd[,i])=="TRUE",i])+100
	glm_X_sd[is.infinite(glm_X_sd[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_X_sd.rolwin20 <- cbind(glm_X_sd.rolwin20,rollapply(glm_X_sd[,i], width = 20, FUN = mean, align = "center", fill = NA))
    colnames(glm_X_sd.rolwin20)[i] <- paste(names(glm_X_sd)[i], ".rolwin20", sep="")
}

## Now, join the the new rolling window tables for all chromosomes
glm.all.rolwin20_sd <- rbind(glm_2L_sd.rolwin20,glm_2R_sd.rolwin20,glm_3L_sd.rolwin20,glm_3R_sd.rolwin20,glm_X_sd.rolwin20)

## Make a rolling window table with just p-value, remove logp and fdr columns 
glm.all.rolwin20_sd.pval <- select(glm.all.rolwin20_sd,-contains("logp"),-contains("fdr"))
## Make a rolling window table with just logp, selects them and append loci
glm.all.rolwin20_sd.logp <- cbind(glm.all.rolwin20_sd.pval[,c(1:2)],select(glm.all.rolwin20_sd,contains("logp")))

## Merge annotations for filtering and highlighting
glm.all.rolwin20.annot <- merge(glm.all.rolwin20, vep, by=c("CHROM","POS"))

## Load spinosad-resistance candidate gene list
spino.cand <- read.table("spino.cand.list.txt", header=FALSE)
names(spino.cand) <- "Gene"

## Now merge to find All SNPs in and around candidate genes
glm.all.rolwin20.annot.spino <- merge(spino.cand, glm.all.rolwin20.annot, by="Gene")




################
### Figure 2 ###
################

## Panel D
# This manhattan plot shows the contrast between S and E samples combined across all TPTs
manh.SvE.all.rolwin20 <- ggplot(glm.all.rolwin20, aes(POS, SvE.logp.rolwin20)) +
	geom_line(alpha = 0.5, colour = "gray80") +  
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though 
	geom_hline(yintercept=min(na.omit(glm.all.rolwin20[glm.all.rolwin20$SvE.fdr.rolwin20 < 0.01,]$SvE.logp.rolwin20)),color="black",linetype="dashed",linewidth=.25) +
	## Highlight points significant in S vs E GLM contrast
	geom_point(data=na.omit(glm.all.rolwin20[glm.all.rolwin20$SvE.logp.rolwin20 > min(na.omit(glm.all.rolwin20[glm.all.rolwin20$SvE.fdr.rolwin20<0.01,]$SvE.logp.rolwin20)),]), aes(POS, SvE.logp.rolwin20), color = "#39568A", size = 1) +
	## Highlight spinosad resistance candidate SNPs
	geom_point(data=glm.all.rolwin20.annot.spino,aes(POS, SvE.logp.rolwin20), color = "#C45C5C", size = 1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm.all.rolwin20$SvE.logp.rolwin20)))) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("E vs S: 20-SNP sliding window")
  
## Panel E
# This manhattan plot shows the contrast between PA and E samples combined across all TPTs
manh.PAvE.all.rolwin20 <- ggplot(glm.all.rolwin20, aes(POS, PAvE.logp.rolwin20)) + 
	geom_line(alpha = 0.5, colour = "gray80") +  
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	geom_hline(yintercept=min(na.omit(glm.all.rolwin20[glm.all.rolwin20$PAvE.fdr.rolwin20 < 0.01,]$PAvE.logp.rolwin20)),color="black",linetype="dashed",linewidth=.25) +
	## Highlight points significant in S vs E GLM contrast
	geom_point(data=na.omit(glm.all.rolwin20[glm.all.rolwin20$SvE.logp.rolwin20 > min(na.omit(glm.all.rolwin20[glm.all.rolwin20$SvE.fdr.rolwin20<0.01,]$SvE.logp.rolwin20)),]),aes(POS, PAvE.logp.rolwin20), color = "#39568A", size = 1) +
	## Highlight spinosad resistance candidate SNPs
	geom_point(data=glm.all.rolwin20.annot.spino,aes(POS, PAvE.logp.rolwin20), color = "#C45C5C", size = 1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm.all.rolwin20$PAvE.logp.rolwin20)))) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("E vs PA: 20-SNP sliding window")
  
## Save them together
pdf(file = "rudflies_2023_redo.SvE.all.PAvE.all.rolwin20.glm.manh.pdf", width=7.5, height=6)
	ggarrange(manh.SvE.all.rolwin20, manh.PAvE.all.rolwin20,
              ncol = 1, nrow = 2)
dev.off()

### AF change side panels ###

## General idea is that for spinosad resistance candidate SNPs we'll be comparing absolute 
## AF differences in "S vs E" and "PA vs E" contrasts and compare those AF differences to 
## those of matched background SNPs. We'll check whether they we see higher AF differences 
## for spinosad candidates than for matched SNPs. This would tell us whether candidate AFs
## are changing in spino-exposed populations more than expected by chance, suggesting 
## possible adaptative response. We'll do the same analysis for "S vs E" and "PA vs E" 
## glm contrast top outliers.

## We'll also check AF means from "S vs E" contrast significant SNPs in the "PA vs E" 
## contrast to see if we can detect some adaptive convergence.

## Start by calculating mean allele frequencies for each treat by TPT and across TPTs
freq_means <- cbind(haf.sites.filt, 
  E.af=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="E"]),
  S.af=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="S"]),
  PA.af=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="PA"]),
  E.T1.af=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="E"&haf.meta.filt$tpt=="1"]),
  E.T2.af=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="E"&haf.meta.filt$tpt=="2"]),
  E.T3.af=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="E"&haf.meta.filt$tpt=="3"]),
  E.T4.af=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="E"&haf.meta.filt$tpt=="4"]),
  S.T1.af=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="S"&haf.meta.filt$tpt=="1"]),
  S.T2.af=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="S"&haf.meta.filt$tpt=="2"]),
  S.T3.af=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="S"&haf.meta.filt$tpt=="3"]),
  S.T4.af=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="S"&haf.meta.filt$tpt=="4"]),
  PA.T1.af=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="PA"&haf.meta.filt$tpt=="1"]),
  PA.T2.af=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="PA"&haf.meta.filt$tpt=="2"]),
  PA.T3.af=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="PA"&haf.meta.filt$tpt=="3"]),
  PA.T4.af=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="PA"&haf.meta.filt$tpt=="4"]))
  
## Merge AF means only with annotation info
freq_means_annot <- merge(freq_means, vep, by=c("CHROM","POS"))

## Make sure we don't have duplicate loci, generated by overlapping annotations
freq_means_annot_nr <- freq_means_annot[!duplicated(freq_means_annot[,c(1:2)]), ]

## Now filter for spino candidates only
freq_means_annot_nr_spino <- freq_means_annot_nr[freq_means_annot_nr$Gene %in% spino.cand$Gene,]

## Now filter for non-spino candidates as a background SNP list
freq_means_annot_nr_bg <- freq_means_annot_nr[!freq_means_annot_nr$Gene %in% spino.cand$Gene,]

## Merge raw GLM results table (pre-rolling window average) with annotation info 
glm.all.annot <- merge(glm.all, vep, by=c("CHROM","POS"))
## remove duplicate
glm.all.annot.nr <- glm.all.annot[!duplicated(glm.all.annot[,c(1:2)]), ]
## Since there are more significant "S vs E" SNPs than all N spinosad candidate SNPs, 
## make same-sized list of the top N SNPs from "S vs E" contrast.
glm.all.annot.nr.topEvS <- head(glm.all.annot.nr[order(glm.all.annot.nr$SvE.fdr), ],
	nrow(freq_means_annot_nr_spino))
## Make AF mean frequency table for top "S vs E" SNP list
freq_means_annot_nr_EvS <- merge(glm.all.annot.nr.topEvS[,c(1:2)], 
								freq_means_annot_nr, 
								by=c("CHROM","POS"))
## Make AF mean frequency table for "S vs E" background SNPs (not top outliers)
freq_means_annot_nr_EvS_bg <- freq_means_annot_nr[!freq_means_annot_nr$Uploaded_variation %in% glm.all.annot.nr.topEvS$Uploaded_variation,]

## Since there are more significant "PA vs E" SNPs than all N spinosad candidate SNPs, 
## make same-sized list of the top N SNPs from "PA vs E" contrast.
glm.all.annot.nr.topEvPA <- head(glm.all.annot.nr[order(glm.all.annot.nr$SvE.fdr), ],nrow(freq_means_annot_nr_spino))
## Make AF mean frequency table for top "PA vs E" SNP list
freq_means_annot_nr_EvPA <- merge(glm.all.annot.nr.topEvPA[,c(1:2)],
								freq_means_annot_nr, 
								by=c("CHROM","POS"))
## Make AF mean frequency table for "PA vs E" background SNPs (not top outliers)
freq_means_annot_nr_EvPA_bg <- freq_means_annot_nr[!freq_means_annot_nr$Uploaded_variation %in% glm.all.annot.nr.topEvPA$Uploaded_variation,]

## The following loop will go through the list of spino candidates, find all matching
## background background SNPs, sample 100 of them, then calculate absolute frequency
## differences for treatment contrasts of interest.

## Instead of running full loop, you can just reload the table from a previous run
freq_means_results <- read.table("rudflies_2023_redo.freq_means_results.txt", header=TRUE)

### DON'T RUN CODE BLOCK IF YOU JUST LOADED TABLE FROM PREVIOUS RUN ###
freq_means_results <- c() #create empty object
## Loop through all cols pval values for candidate genes
for(g in c(1:nrow(freq_means_annot_nr_spino))) { #number of sites in spino candidate genes
  tryCatch({
  	temp_chrom <- freq_means_annot_nr_spino[g,1] #find chrom
  	temp_e <- freq_means_annot_nr_spino[g,3] #find E mean frequency
  	temp_s <- freq_means_annot_nr_spino[g,4] #find S mean frequency
  	temp_pa <- freq_means_annot_nr_spino[g,5] #find PA mean frequency
  	temp_cons <- freq_means_annot_nr_spino[g,24] #find allele consequence
  	#select bg SNPs by matching chromosome, consequence, and approximate size
  	temp_bg <- freq_means_annot_nr_bg[freq_means_annot_nr_bg$CHROM %in% temp_chrom[1] &
  		freq_means_annot_nr_bg$Consequence %in% temp_cons &
  		freq_means_annot_nr_bg$E.af>(temp_e - (0.25 * min(temp_e,1-temp_e))) &
  		freq_means_annot_nr_bg$E.af<(temp_e + (0.25 * min(temp_e,1-temp_e))),]
  	#randomly sample 100 SNPs, or use all samples if less than 100 matches
  	temp_bg <- temp_bg[sample(nrow(temp_bg), min(c(nrow(temp_bg),100))), ] 
  	#calculate bg allele differences
  	temp_af <- cbind(bg.E2S=abs(temp_bg$E.af-temp_bg$S.af), #raw diff for S vs E
  		bg.E2PA=abs(temp_bg$E.af-temp_bg$PA.af), #raw diff for PA vs E
  		bg.E2S.perc=abs((temp_bg$E.af-temp_bg$S.af)/temp_bg$E.af), #perc diff for S vs E
  		bg.E2PA.perc=abs((temp_bg$E.af-temp_bg$PA.af)/temp_bg$E.af))#perc diff for PA vs E
  	#calculate means for all bg allele differences
  	temp_means <- colMeans(temp_af) #find means for all bg SNP differences
  	names(temp_means) <- paste0(names(temp_means), ".mean") 
  	#calculate variances for all bg allele differences
  	temp_vars <- colVars(temp_af)
  	names(temp_vars) <- paste0(names(temp_vars), ".var")
  	#add row to table
  	freq_means_results <- rbind(freq_means_results,append(temp_means,temp_vars))
  }, error=function(e){})
}  

#save table so we don't need to run this loop every time
colnames(freq_means_results) <- paste0(colnames(freq_means_results), ".spino")
write.table(freq_means_results, file="rudflies_2023_redo.freq_means_results.txt",sep = "\t", quote = FALSE, row.names = F)


## The following loop will go through the list of "S vs E" top outliers, find all matching
## background background SNPs, sample 100 of them, then calculate absolute frequency
## differences for treatment contrasts of interest.

## Instead of running full loop, you can just reload the table from a previous run
freq_means_results_EvS <- read.table("rudflies_2023_redo.freq_means_results_EvS.txt", header=TRUE)

### DON'T RUN CODE BLOCK IF YOU JUST LOADED TABLE FROM PREVIOUS RUN ###
freq_means_results_EvS <- c() #create empty object
## Loop through all cols pval values for candidate genes
for(g in c(1:nrow(freq_means_annot_nr_spino))) { #number of sites in spino candidate genes
  tryCatch({
  	temp_chrom <- freq_means_annot_nr_EvS[g,1] #find chrom
  	temp_e <- freq_means_annot_nr_EvS[g,3] #find E mean frequency
  	temp_s <- freq_means_annot_nr_EvS[g,4] #find S mean frequency
  	temp_pa <- freq_means_annot_nr_EvS[g,5] #find PA mean frequency
  	temp_cons <- freq_means_annot_nr_EvS[g,24] #find allele consequence
  	#select bg SNPs by matching chromosome, consequence, and approximate size
  	temp_bg <- freq_means_annot_nr_EvS_bg[freq_means_annot_nr_EvS_bg$CHROM %in% temp_chrom[1] &
  		freq_means_annot_nr_EvS_bg$Consequence %in% temp_cons &
  		freq_means_annot_nr_EvS_bg$E.af>(temp_e - (0.25 * min(temp_e,1-temp_e))) &
  		freq_means_annot_nr_EvS_bg$E.af<(temp_e + (0.25 * min(temp_e,1-temp_e))),]
  	#randomly sample 100 SNPs, or use all samples if less than 100 matches
  	temp_bg <- temp_bg[sample(nrow(temp_bg), min(c(nrow(temp_bg),100))), ] 
  	#calculate bg allele differences
  	temp_af <- cbind(bg.E2S=abs(temp_bg$E.af-temp_bg$S.af), #raw diff for S vs E
  		bg.E2PA=abs(temp_bg$E.af-temp_bg$PA.af), #raw diff for PA vs E
  		bg.E2S.perc=abs((temp_bg$E.af-temp_bg$S.af)/temp_bg$E.af), #perc diff for S vs E
  		bg.E2PA.perc=abs((temp_bg$E.af-temp_bg$PA.af)/temp_bg$E.af))#perc diff for PA vs E
  	#calculate means for all bg allele differences
  	temp_means <- colMeans(temp_af)
  	names(temp_means) <- paste0(names(temp_means), ".mean")
  	#calculate variances for all bg allele differences
  	temp_vars <- colVars(temp_af)
  	names(temp_vars) <- paste0(names(temp_vars), ".var")
  	#add row to table
  	freq_means_results_EvS <- rbind(freq_means_results_EvS,append(temp_means,temp_vars))
  }, error=function(e){})
}  

#save table so we don't need to run this loop every time
colnames(freq_means_results_EvS) <- paste0(colnames(freq_means_results_EvS), ".EvS")
write.table(freq_means_results_EvS, file="rudflies_2023_redo.freq_means_results_EvS.txt",sep = "\t", quote = FALSE, row.names = F)


## The following loop will go through the list of "PA vs E" top outliers, find all 
## matching background background SNPs, sample 100 of them, then calculate absolute 
## frequency differences for treatment contrasts of interest.

## Instead of running full loop, you can just reload the table from a previous run
freq_means_results_EvPA <- read.table("rudflies_2023_redo.freq_means_results_EvPA.txt", header=TRUE)

### DON'T RUN CODE BLOCK IF YOU JUST LOADED TABLE FROM PREVIOUS RUN ###
freq_means_results_EvPA <- c() #create empty object
## Loop through all cols pval values for candidate genes
for(g in c(1:nrow(freq_means_annot_nr_spino))) { #number of sites in spino candidate genes
  tryCatch({
  	temp_chrom <- freq_means_annot_nr_EvPA[g,1] #find chrom
  	temp_e <- freq_means_annot_nr_EvPA[g,3] #find E mean frequency
  	temp_s <- freq_means_annot_nr_EvPA[g,4] #find S mean frequency
  	temp_pa <- freq_means_annot_nr_EvPA[g,5] #find PA mean frequency
  	temp_cons <- freq_means_annot_nr_EvPA[g,24] #find allele consequence
  	#select bg SNPs by matching chromosome, consequence, and approximate size
  	temp_bg <- freq_means_annot_nr_EvPA_bg[freq_means_annot_nr_EvPA_bg$CHROM %in% temp_chrom[1] &
  		freq_means_annot_nr_EvPA_bg$Consequence %in% temp_cons &
  		freq_means_annot_nr_EvPA_bg$E.af>(temp_e - (0.25 * min(temp_e,1-temp_e))) &
  		freq_means_annot_nr_EvPA_bg$E.af<(temp_e + (0.25 * min(temp_e,1-temp_e))),] 
  	#randomly sample 100 SNPs, or use all samples if less than 100 matches
  	temp_bg <- temp_bg[sample(nrow(temp_bg), min(c(nrow(temp_bg),100))), ]
  	#calculate bg allele differences
  	temp_af <- cbind(bg.E2S=abs(temp_bg$E.af-temp_bg$S.af),
  		bg.E2PA=abs(temp_bg$E.af-temp_bg$PA.af),
  		bg.E2S.perc=abs((temp_bg$E.af-temp_bg$S.af)/temp_bg$E.af),
  		bg.E2PA.perc=abs((temp_bg$E.af-temp_bg$PA.af)/temp_bg$E.af))
  	#calculate means for all bg allele differences
  	temp_means <- colMeans(temp_af)
  	names(temp_means) <- paste0(names(temp_means), ".mean")
  	#calculate variances for all bg allele differences
  	temp_vars <- colVars(temp_af)
  	names(temp_vars) <- paste0(names(temp_vars), ".var")
  	#add row to table
  	freq_means_results_EvPA <- rbind(freq_means_results_EvPA,append(temp_means,temp_vars))
  }, error=function(e){})
}  

#save table so we don't need to run this loop every time
colnames(freq_means_results_EvPA) <- paste0(colnames(freq_means_results_EvPA), ".EvPA")
write.table(freq_means_results_EvPA, file="rudflies_2023_redo.freq_means_results_EvPA.txt",sep = "\t", quote = FALSE, row.names = F)


### Now that we've calculated background AF difference, time to calculate the same for 
### actual spino candidates and GLM outlier SNPs, then make a master table for each.

## First, make spino candidate table
freq_means_annot_nr_spino_afdiff <- cbind(freq_means_annot_nr_spino,
	cand.E2S.spino=abs(freq_means_annot_nr_spino$E.af-freq_means_annot_nr_spino$S.af),
	cand.E2PA.spino=abs(freq_means_annot_nr_spino$E.af-freq_means_annot_nr_spino$PA.af),
	cand.E2S.perc.spino=abs((freq_means_annot_nr_spino$E.af - freq_means_annot_nr_spino$S.af)/freq_means_annot_nr_spino$E.af),
	cand.E2PA.perc.spino=abs((freq_means_annot_nr_spino$E.af - freq_means_annot_nr_spino$PA.af)/freq_means_annot_nr_spino$E.af),
	freq_means_results)

## Second, make "E vs S" GLM outlier table
freq_means_annot_nr_EvS_afdiff <- cbind(freq_means_annot_nr_EvS,
	cand.E2S.EvS=abs(freq_means_annot_nr_EvS$E.af - freq_means_annot_nr_EvS$S.af),
	cand.E2PA.EvS=abs(freq_means_annot_nr_EvS$E.af - freq_means_annot_nr_EvS$PA.af), cand.E2S.perc.EvS=abs((freq_means_annot_nr_EvS$E.af - freq_means_annot_nr_EvS$S.af)/freq_means_annot_nr_EvS$E.af),
	cand.E2PA.perc.EvS=abs((freq_means_annot_nr_EvS$E.af - freq_means_annot_nr_EvS$PA.af)/freq_means_annot_nr_EvS$E.af),
	freq_means_results_EvS)

## Third, make "E vs PA" GLM outlier table
freq_means_annot_nr_EvPA_afdiff <- cbind(freq_means_annot_nr_EvPA,
	cand.E2S.EvPA=abs(freq_means_annot_nr_EvPA$E.af-freq_means_annot_nr_EvPA$S.af),
	cand.E2PA.EvPA=abs(freq_means_annot_nr_EvPA$E.af-freq_means_annot_nr_EvPA$PA.af),
	cand.E2S.perc.EvPA=abs((freq_means_annot_nr_EvPA$E.af - freq_means_annot_nr_EvPA$S.af)/freq_means_annot_nr_EvPA$E.af),
	cand.E2PA.perc.EvPA=abs((freq_means_annot_nr_EvPA$E.af - freq_means_annot_nr_EvPA$PA.af)/freq_means_annot_nr_EvPA$E.af),
	freq_means_results_EvPA)


## Combine all relevant AF difference stats into a master table
freq_means_annot_nr_afdiff <- cbind(freq_means_annot_nr_spino_afdiff[,c(33:ncol(freq_means_annot_nr_spino_afdiff))], freq_means_annot_nr_EvS_afdiff[,c(33:ncol(freq_means_annot_nr_EvS_afdiff))], freq_means_annot_nr_EvPA_afdiff[,c(33:ncol(freq_means_annot_nr_EvPA_afdiff))])


## Reorganize master table for plots: all AF diff stats in one column plus metadata cols
## E vs S differences only
freq_means_annot_nr_E2Sdiff <- data.frame(rbind(cbind("candidate","E2S","spino","raw",freq_means_annot_nr_afdiff$cand.E2S.spino),
	   cbind("candidate","E2S","EvS","raw", freq_means_annot_nr_afdiff$cand.E2S.EvS),
	   cbind("matched","E2S","spino","mean", freq_means_annot_nr_afdiff$bg.E2S.mean.spino),
	   cbind("matched","E2S","EvS","mean", freq_means_annot_nr_afdiff$bg.E2S.mean.EvS),
	   cbind("matched","E2S","spino","perc", freq_means_annot_nr_afdiff$bg.E2S.perc.mean.spino),
	   cbind("matched","E2S","EvS","perc", freq_means_annot_nr_afdiff$bg.E2S.perc.mean.EvS)))
	   
## Rename cols
colnames(freq_means_annot_nr_E2Sdiff) <- c("set","contrast","category","type","af_diff")
## Set AF differences to numeric format
freq_means_annot_nr_E2Sdiff[,5] <- as.numeric(freq_means_annot_nr_E2Sdiff[,5])

##Calculate means of AF differences by set (cand or matched), category (spino or "E vs S" 
## outliers), and type (raw, mean, or perc)
freq_means_annot_nr_E2Sdiff %>%
  group_by(set, category, type) %>%
  summarise(avg = mean(af_diff), med = median(af_diff), stdev = sd(af_diff), .groups="keep")

#T-test for differences in absolute E vs S AF-change in E vs S outliers vs background
t.test(af_diff ~ set, 
	data = freq_means_annot_nr_E2Sdiff[freq_means_annot_nr_E2Sdiff$type!="perc" &
	freq_means_annot_nr_E2Sdiff$category=="EvS",],
	alternative = 'greater')

#T-test for differences in absolute E vs S AF-change in spino candidates vs background
t.test(af_diff ~ set, 
	data = freq_means_annot_nr_E2Sdiff[freq_means_annot_nr_E2Sdiff$type!="perc" &
	freq_means_annot_nr_E2Sdiff$category=="spino",],
	alternative = 'greater')

### Panel F ###
# make plot of percent AF differences by set and category
E2Sdiff <- ggplot(data = freq_means_annot_nr_E2Sdiff[freq_means_annot_nr_E2Sdiff$type!="perc",], aes(x=set, y=af_diff)) + 
	geom_boxplot(aes(colour=category)) + 
	scale_color_manual(values=c("#39568A","#C45C5C")) +
	ggtitle("AF differences: E vs S") + 
	theme_classic()

## Reorganize master table for plots: all AF diff stats in one column plus metadata cols
## E vs PA differences only
freq_means_annot_nr_E2PAdiff <- data.frame(rbind(cbind("candidate","E2PA","spino","raw",freq_means_annot_nr_afdiff$cand.E2PA.spino),
	   cbind("candidate","E2PA","EvS","raw",freq_means_annot_nr_afdiff$cand.E2PA.EvS),
	   cbind("matched","E2PA","spino","mean",freq_means_annot_nr_afdiff$bg.E2PA.mean.spino),
	   cbind("matched","E2PA","EvS","mean",freq_means_annot_nr_afdiff$bg.E2PA.mean.EvS),
	   cbind("matched","E2PA","spino","perc",freq_means_annot_nr_afdiff$bg.E2PA.perc.mean.spino),
	   cbind("matched","E2PA","EvS","perc",freq_means_annot_nr_afdiff$bg.E2PA.perc.mean.EvS)))

## Rename cols
colnames(freq_means_annot_nr_E2PAdiff) <- c("set","contrast","category","type","af_diff")
## Set AF differences to numeric format
freq_means_annot_nr_E2PAdiff[,5] <- as.numeric(freq_means_annot_nr_E2PAdiff[,5])

##Calculate means of AF differences by set (cand or matched), category (spino or "E vs S" 
## outliers), and type (raw, mean, or perc)
freq_means_annot_nr_E2PAdiff %>%
  group_by(set, category, type) %>%
  summarise(avg = mean(af_diff), med = median(af_diff), stdev = sd(af_diff), .groups="keep")

#T-test for differences in absolute E vs PA AF-change in E vs S outliers vs background
t.test(af_diff ~ set, 
	data = freq_means_annot_nr_E2PAdiff[freq_means_annot_nr_E2PAdiff$type!="perc" &
	freq_means_annot_nr_E2PAdiff$category=="EvS",],
	alternative = 'greater')

#T-test for differences in absolute E vs PA AF-change in spino candidates vs background
t.test(af_diff ~ set, 
	data = freq_means_annot_nr_E2PAdiff[freq_means_annot_nr_E2PAdiff$type!="perc" &
	freq_means_annot_nr_E2PAdiff$category=="spino",],
	alternative = 'greater')


### Panel G ###
# make plot of percent AF differences by set and category
E2PAdiff <- ggplot(data = freq_means_annot_nr_E2PAdiff[freq_means_annot_nr_E2PAdiff$type!="perc",], aes(x=set, y=af_diff)) + 
	geom_boxplot(aes(colour=category)) + 
	scale_color_manual(values=c("#39568A","#C45C5C")) +
	ggtitle("AF differences: E vs PA") + 
	theme_classic()


## Save panel F & G plots together
pdf(file = "rudflies_2023_redo.freq_means_annot_nr_E2Sdiff.mean.multiboxplot.pdf", width=3, height=6)
	ggarrange(E2Sdiff, E2PAdiff,
              ncol = 1, nrow = 2)
dev.off()



################
### Figure 3 ###
################

## Panel D
# This manhattan plot shows the contrast between SE and E samples combined in TPT 1
manh.EvSE.T1.rolwin20 <- ggplot(glm.all.rolwin20, aes(POS, EvSE.T1.logp.rolwin20)) +
	geom_line(alpha = 0.5, colour = "gray80") +
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	geom_hline(yintercept=min(na.omit(glm.all.rolwin20[glm.all.rolwin20$EvSE.T1.fdr.rolwin20 < 0.01,]$EvSE.T1.logp.rolwin20)),color="black",linetype="dashed",linewidth=.25) +
	## Highlight points significant in SE vs E GLM contrast
	geom_point(data=na.omit(glm.all.rolwin20[glm.all.rolwin20$EvSE.T1.logp.rolwin20 > min(na.omit(glm.all.rolwin20[glm.all.rolwin20$EvSE.T1.fdr.rolwin20<0.01,]$EvSE.T1.logp.rolwin20)),]),aes(POS, EvSE.T1.logp.rolwin20), color = "#39568A", size = 1) +
	## Highlight spinosad resistance candidate SNPs
	geom_point(data=glm.all.rolwin20.annot.spino,aes(POS, EvSE.T1.logp.rolwin20), colour = "#C45C5C", size = 1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm.all.rolwin20$EvSE.T1.logp.rolwin20)))) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("E vs SE (TPT1): 20-SNP sliding window")

## Unused in figure 3, relegated to Fig. S4 panel E
# This manhattan plot shows the contrast between SP and E samples combined in TPT 1
manh.EvSP.T1.rolwin20 <- ggplot(glm.all.rolwin20, aes(POS, EvSP.T1.logp.rolwin20)) + 
	geom_line(alpha = 0.5, colour = "gray80") +  
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	geom_hline(yintercept=min(na.omit(glm.all.rolwin20[glm.all.rolwin20$EvSP.T1.fdr.rolwin20 < 0.01,]$EvSP.T1.logp.rolwin20)),color="black",linetype="dashed",linewidth=.25) +
	## Highlight points significant in SE vs E GLM contrast
	geom_point(data=na.omit(glm.all.rolwin20[glm.all.rolwin20$EvSE.T1.logp.rolwin20>min(na.omit(glm.all.rolwin20[glm.all.rolwin20$EvSE.T1.fdr.rolwin20<0.01,]$EvSE.T1.logp.rolwin20)),]),aes(POS, EvSP.T1.logp.rolwin20), color = "#39568A", size = 1) +
	## Highlight spinosad resistance candidate SNPs
geom_point(data=glm.all.rolwin20.annot.spino,aes(POS, EvSP.T1.logp.rolwin20), colour = "#C45C5C", size = 1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm.all.rolwin20$EvSP.T1.logp.rolwin20))+0.5)) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("E vs SP (TPT1): 20-SNP sliding window")

## Unused in figure 3
# This manhattan plot shows the contrast between PA and E samples combined in TPT 1
manh.SEvSP.T1.rolwin20 <- ggplot(glm.all.rolwin20, aes(POS, SEvSP.T1.logp.rolwin20)) + 
	geom_line(alpha = 0.5, colour = "gray80") +  
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	geom_hline(yintercept=min(na.omit(glm.all.rolwin20[glm.all.rolwin20$SEvSP.T1.fdr.rolwin20 < 0.01,]$SEvSP.T1.logp.rolwin20)),color="black",linetype="dashed",linewidth=.25) +
	## Highlight points significant in SE vs E GLM contrast
	geom_point(data=na.omit(glm.all.rolwin20[glm.all.rolwin20$EvSE.T1.logp.rolwin20 > min(na.omit(glm.all.rolwin20[glm.all.rolwin20$EvSE.T1.fdr.rolwin20 < 0.01,]$EvSE.T1.logp.rolwin20)),]),aes(POS, SEvSP.T1.logp.rolwin20), color = "#39568A", size = 1) +
	## Highlight spinosad resistance candidate SNPs
	geom_point(data=glm.all.rolwin20.annot.spino,aes(POS, SEvSP.T1.logp.rolwin20), colour = "#C45C5C", size = 1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm.all.rolwin20$SEvSP.T1.logp.rolwin20)))) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("SE vs SP (TPT1): 20-SNP sliding window")

## Panel E
# This manhattan plot shows the contrast between PA and E samples combined in TPT 1
manh.EvPA.T1.rolwin20 <- ggplot(glm.all.rolwin20, aes(POS, EvPA.T1.logp.rolwin20)) + 
	geom_line(alpha = 0.5, colour = "gray80") +  
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	geom_hline(yintercept=min(na.omit(glm.all.rolwin20[glm.all.rolwin20$EvPA.T1.fdr.rolwin20 < 0.01,]$EvPA.T1.logp.rolwin20)),color="black",linetype="dashed",linewidth=.25) +
	## Highlight points significant in SE vs E GLM contrast
	geom_point(data=na.omit(glm.all.rolwin20[glm.all.rolwin20$EvSE.T1.logp.rolwin20 > min(na.omit(glm.all.rolwin20[glm.all.rolwin20$EvSE.T1.fdr.rolwin20<0.01,]$EvSE.T1.logp.rolwin20)),]),aes(POS, EvPA.T1.logp.rolwin20), color = "#39568A", size = 1) +
	## Highlight spinosad resistance candidate SNPs
	geom_point(data=glm.all.rolwin20.annot.spino,aes(POS, EvPA.T1.logp.rolwin20), colour = "#C45C5C", size = 1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm.all.rolwin20$EvPA.T1.logp.rolwin20)))) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("E vs PA (TPT1): 20-SNP sliding window")

## Plot all four pairwise TPT1 treatment comparisons together
pdf(file = "rudflies_2023_redo.EvSE.EvSP.SEvSP.EvPA.T1.rolwin20.glm.manh.pdf", width=7.5, height=14)
	ggarrange(manh.EvSE.T1.rolwin20, manh.EvSP.T1.rolwin20, manh.SEvSP.T1.rolwin20, manh.EvPA.T1.rolwin20,
              ncol = 1, nrow = 4)
dev.off()

## Plot three most relevant pairwise TPT1 treatment comparisons (not SP vs SE) together
pdf(file = "rudflies_2023_redo.EvSE.EvSP.EvPA.T1.rolwin20.glm.manh.pdf", width=7.5, height=9)
	ggarrange(manh.EvSE.T1.rolwin20, manh.EvSP.T1.rolwin20, manh.EvPA.T1.rolwin20,
              ncol = 1, nrow = 3)
dev.off()


### AF change side panels ###

## This time, for spinosad resistance candidate SNPs, we'll be comparing absolute 
## AF differences in TPT1 pairwise treatment contrasts and compare those AF differences to 
## those of matched background SNPs. We'll check whether they we see higher AF differences 
## for spinosad candidates than for matched SNPs. This would tell us whether candidate AFs
## are changing in spino-exposed populations more than expected by chance, suggesting 
## possible adaptative response. We'll do the same analysis for "SE vs E" and "PA vs E" 
## glm contrast top outliers. 

## In TPT1, the S treatment was divided into SE (extinct cages) and SP (persistent cages)

## We'll also check AF means from "SE vs E" contrast significant SNPs in other contrasts 
## to see if we can detect some adaptive convergence. We would have done the same with 
## "SP vs E" outliers, but there were none.

## Start by calculating mean allele frequencies for each treatment at TPT1
freq_means_t1 <- cbind(haf.sites.T1filt,
	E.af=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$treat.fix == "E" & haf.meta.T1filt$tpt=="1"]),
	S.af=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$treat.fix == "S" & haf.meta.T1filt$tpt=="1"]),
	SE.af=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition == "SE" & haf.meta.T1filt$tpt=="1"]),
	SP.af=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition == "SP" & haf.meta.T1filt$tpt=="1"]),
	PA.af=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$treat.fix == "PA"&haf.meta.T1filt$tpt=="1"]))


## Merge AF means only with annotation info
freq_means_t1_annot <- merge(freq_means_t1, vep, by=c("CHROM","POS"))

## Make sure we don't have duplicate loci, generated by overlapping annotations
freq_means_t1_annot_nr <- freq_means_t1_annot[!duplicated(freq_means_t1_annot[,c(1:2)]), ]

## Now filter for spino candidates only
freq_means_t1_annot_nr_spino <- freq_means_t1_annot_nr[freq_means_t1_annot_nr$Gene %in% spino.cand$Gene,]

## Now filter for non-spino candidates as a background SNP list
freq_means_t1_annot_nr_bg <- freq_means_t1_annot_nr[!freq_means_t1_annot_nr$Gene %in% spino.cand$Gene,]

## Since the number of sig "SE vs E" SNPs is different from all N spinosad candidate SNPs, 
## make same-sized list of the top N SNPs from "SE vs E" contrast.
glm.all.annot.nr.topEvSE <- head(glm.all.annot.nr[order(glm.all.annot.nr$EvSE.T1.fdr), ],
	nrow(freq_means_t1_annot_nr_spino))
## Make AF mean frequency table for top "SE vs E" SNP list
freq_means_t1_annot_nr_EvSE <- merge(glm.all.annot.nr.topEvSE[,c(1:2)],
								freq_means_t1_annot_nr,
								by=c("CHROM","POS"))
## Make AF mean frequency table for "SE vs E" background SNPs (not top outliers)
freq_means_t1_annot_nr_EvSE_bg <- freq_means_t1_annot_nr[!freq_means_t1_annot_nr$Uploaded_variation %in% glm.all.annot.nr.topEvSE$Uploaded_variation,]

## Since the number of sig "PA vs E" SNPs is different from all N spinosad candidate SNPs,
## make same-sized list of the top N SNPs from "PA vs E" contrast.
glm.all.annot.nr.topEvPA <- head(glm.all.annot.nr[order(glm.all.annot.nr$PAvE.T1.fdr), ],nrow(freq_means_t1_annot_nr_spino))
## Make AF mean frequency table for top "PA vs E" SNP list
freq_means_t1_annot_nr_EvPA <- merge(glm.all.annot.nr.topEvPA[,c(1:2)],
								freq_means_t1_annot_nr,
								by=c("CHROM","POS"))
## Make AF mean frequency table for "PA vs E" background SNPs (not top outliers)
freq_means_t1_annot_nr_EvPA_bg <- freq_means_t1_annot_nr[!freq_means_t1_annot_nr$Uploaded_variation %in% glm.all.annot.nr.topEvPA$Uploaded_variation,]

## Since there are no sig "SP vs E" SNPs, make a list of the top N SNPs from "SP vs E" 
## contrast to match all N spinosad candidate SNPs.
glm.all.annot.nr.topEvSP <- head(glm.all.annot.nr[order(glm.all.annot.nr$EvSP.T1.fdr), ],
	nrow(freq_means_t1_annot_nr_spino))
## Make AF mean frequency table for top "SP vs E" SNP list
freq_means_t1_annot_nr_EvSP <- merge(glm.all.annot.nr.topEvSP[,c(1:2)],
								freq_means_t1_annot_nr,
								by=c("CHROM","POS"))
## Make AF mean frequency table for "SP vs E" background SNPs (not top outliers)
freq_means_t1_annot_nr_EvSP_bg <- freq_means_t1_annot_nr[!freq_means_t1_annot_nr$Uploaded_variation %in% glm.all.annot.nr.topEvSP$Uploaded_variation,]


## Instead of running full loop, you can just reload the table from a previous run
freq_means_t1_results <- read.table("rudflies_2023_redo.freq_means_t1_results.txt", header=TRUE)

### DON'T RUN CODE BLOCK IF YOU JUST LOADED TABLE FROM PREVIOUS RUN ###
freq_means_t1_results <- c()#create empty object
## Loop through all cols pval values for candidate genes
for(g in c(1:nrow(freq_means_t1_annot_nr_spino))) { #number of sites in spino candidate genes
  tryCatch({
  	temp_chrom <- freq_means_t1_annot_nr_spino[g,1] #find chrom
  	temp_e <- freq_means_t1_annot_nr_spino[g,3] #find E mean frequency
  	temp_se <- freq_means_t1_annot_nr_spino[g,5] #find SE mean frequency
  	temp_sp <- freq_means_t1_annot_nr_spino[g,6] #find SP mean frequency
  	temp_pa <- freq_means_t1_annot_nr_spino[g,7] #find PA mean frequency
  	temp_cons <- freq_means_t1_annot_nr_spino[g,14] #find allele consequence
  	#select bg SNPs by matching chromosome, consequence, and approximate size
  	temp_bg <- freq_means_t1_annot_nr_bg[freq_means_t1_annot_nr_bg$CHROM %in% temp_chrom[1] &
  		freq_means_t1_annot_nr_bg$Consequence %in% temp_cons &
  		freq_means_t1_annot_nr_bg$E.af>(temp_e - (0.25 * min(temp_e,1-temp_e))) &
  		freq_means_t1_annot_nr_bg$E.af<(temp_e + (0.25 * min(temp_e,1-temp_e))),]
  	temp_bg <- temp_bg[sample(nrow(temp_bg), min(c(nrow(temp_bg),100))), ] #randomly sample 100 SNPs
  	#calculate bg allele differences
  	temp_af <- cbind(bg.E2SE=abs(temp_bg$E.af-temp_bg$SE.af),  #raw diff for SE vs E
  		bg.E2SP=abs(temp_bg$E.af-temp_bg$SP.af),  #raw diff for SP vs E
  		bg.E2PA=abs(temp_bg$E.af-temp_bg$PA.af), #raw diff for PA vs E
  		bg.E2SE.perc=abs((temp_bg$E.af-temp_bg$SE.af)/temp_bg$E.af),#perc diff for SE vs E
  		bg.E2SP.perc=abs((temp_bg$E.af-temp_bg$SP.af)/temp_bg$E.af),#perc diff for SP vs E
  		bg.E2PA.perc=abs((temp_bg$E.af-temp_bg$PA.af)/temp_bg$E.af))#perc diff for PA vs E
  	#calculate means for all bg allele differences
  	temp_means <- colMeans(temp_af) #find means for all bg SNP differences
  	names(temp_means) <- paste0(names(temp_means), ".mean")
  	#calculate variances for all bg allele differences
  	temp_vars <- colVars(temp_af)
  	names(temp_vars) <- paste0(names(temp_vars), ".var")
  	#add row to table
  	freq_means_t1_results <- rbind(freq_means_t1_results,append(temp_means,temp_vars))
  }, error=function(e){})
}

#save table so we don't need to run this loop every time
colnames(freq_means_t1_results) <- paste0(colnames(freq_means_t1_results), ".spino")
write.table(freq_means_t1_results, file="rudflies_2023_redo.freq_means_t1_results.txt",sep = "\t", quote = FALSE, row.names = F)

## The following loop will go through list of "SE vs E" top outliers, find all matching
## background background SNPs, sample 100 of them, then calculate absolute frequency
## differences for treatment contrasts of interest.

## Instead of running full loop, you can just reload the table from a previous run
freq_means_t1_results_EvSE <- read.table("rudflies_2023_redo.freq_means_t1_results_EvSE.txt", header=TRUE)

### DON'T RUN CODE BLOCK IF YOU JUST LOADED TABLE FROM PREVIOUS RUN ###
freq_means_t1_results_EvSE <- c()#create empty object
## Loop through all cols pval values for candidate genes
for(g in c(1:nrow(freq_means_t1_annot_nr_EvSE))) { #number of sites in EvSE candidate genes
  tryCatch({
  	temp_chrom <- freq_means_t1_annot_nr_EvSE[g,1] #find chrom
  	temp_e <- freq_means_t1_annot_nr_EvSE[g,3] #find E mean frequency
  	temp_se <- freq_means_t1_annot_nr_EvSE[g,5] #find SE mean frequency
  	temp_sp <- freq_means_t1_annot_nr_EvSE[g,6] #find SP mean frequency
  	temp_pa <- freq_means_t1_annot_nr_EvSE[g,7] #find PA mean frequency
  	temp_cons <- freq_means_t1_annot_nr_EvSE[g,14] #find allele consequence
  	#select bg SNPs by matching chromosome, consequence, and approximate size
  	temp_bg <- freq_means_t1_annot_nr_bg[freq_means_t1_annot_nr_bg$CHROM %in% temp_chrom[1] &
  		freq_means_t1_annot_nr_bg$Consequence %in% temp_cons &
  		freq_means_t1_annot_nr_bg$E.af>(temp_e - (0.25 * min(temp_e,1-temp_e))) &
  		freq_means_t1_annot_nr_bg$E.af<(temp_e + (0.25 * min(temp_e,1-temp_e))),]
  	#randomly sample 100 SNPs, or use all samples if less than 100 matches
  	temp_bg <- temp_bg[sample(nrow(temp_bg), min(c(nrow(temp_bg),100))), ]
  	#calculate bg allele differences
  	temp_af <- cbind(bg.E2SE=abs(temp_bg$E.af-temp_bg$SE.af), #raw diff for SE vs E
  		bg.E2SP=abs(temp_bg$E.af-temp_bg$SP.af), #raw diff for SP vs E
  		bg.E2PA=abs(temp_bg$E.af-temp_bg$PA.af), #raw diff for PA vs E
  		bg.E2SE.perc=abs((temp_bg$E.af-temp_bg$SE.af)/temp_bg$E.af),#perc diff for SE vs E
  		bg.E2SP.perc=abs((temp_bg$E.af-temp_bg$SP.af)/temp_bg$E.af),#perc diff for SP vs E
  		bg.E2PA.perc=abs((temp_bg$E.af-temp_bg$PA.af)/temp_bg$E.af))#perc diff for PA vs E
  	#calculate means for all bg allele differences
  	temp_means <- colMeans(temp_af)
  	names(temp_means) <- paste0(names(temp_means), ".mean")
  	#calculate variances for all bg allele differences
  	temp_vars <- colVars(temp_af)
  	names(temp_vars) <- paste0(names(temp_vars), ".var")
  	#add row to table
  	freq_means_t1_results_EvSE <- rbind(freq_means_t1_results_EvSE,append(temp_means,temp_vars))
  }, error=function(e){})
}

#save table so we don't need to run this loop every time
colnames(freq_means_t1_results_EvSE) <- paste0(colnames(freq_means_t1_results_EvSE), ".EvSE.T1")
write.table(freq_means_t1_results_EvSE, file="rudflies_2023_redo.freq_means_t1_results_EvSE.txt",sep = "\t", quote = FALSE, row.names = F)


## The following loop will go through list of "SP vs E" top outliers, find all matching
## background background SNPs, sample 100 of them, then calculate absolute frequency
## differences for treatment contrasts of interest.

## Instead of running full loop, you can just reload the table from a previous run
freq_means_t1_results_EvSP <- read.table("rudflies_2023_redo.freq_means_t1_results_EvSP.txt", header=TRUE)

### DON'T RUN CODE BLOCK IF YOU JUST LOADED TABLE FROM PREVIOUS RUN ###
freq_means_t1_results_EvSP <- c()#create empty object
## Loop through all cols pval values for candidate genes
for(g in c(1:nrow(freq_means_t1_annot_nr_EvSP))) { #number of sites in EvSP candidate genes
  tryCatch({
  	temp_chrom <- freq_means_t1_annot_nr_EvSP[g,1] #find chrom
  	temp_e <- freq_means_t1_annot_nr_EvSP[g,3] #find E mean frequency
  	temp_se <- freq_means_t1_annot_nr_EvSP[g,5] #find SE mean frequency
  	temp_sp <- freq_means_t1_annot_nr_EvSP[g,6] #find SP mean frequency
  	temp_pa <- freq_means_t1_annot_nr_EvSP[g,7] #find PA mean frequency
  	temp_cons <- freq_means_t1_annot_nr_EvSP[g,14] #find allele consequence
  	#select bg SNPs by matching chromosome, consequence, and approximate size
  	temp_bg <- freq_means_t1_annot_nr_bg[freq_means_t1_annot_nr_bg$CHROM %in% temp_chrom[1] &
  		freq_means_t1_annot_nr_bg$Consequence %in% temp_cons &
  		freq_means_t1_annot_nr_bg$E.af>(temp_e - (0.25 * min(temp_e,1-temp_e))) &
  		freq_means_t1_annot_nr_bg$E.af<(temp_e + (0.25 * min(temp_e,1-temp_e))),]
  	#randomly sample 100 SNPs, or use all samples if less than 100 matches
  	temp_bg <- temp_bg[sample(nrow(temp_bg), min(c(nrow(temp_bg),100))), ]
  	#calculate bg allele differences
  	temp_af <- cbind(bg.E2SE=abs(temp_bg$E.af-temp_bg$SE.af), #raw diff for SE vs E
  		bg.E2SP=abs(temp_bg$E.af-temp_bg$SP.af), #raw diff for SP vs E
  		bg.E2PA=abs(temp_bg$E.af-temp_bg$PA.af), #raw diff for PA vs E
  		bg.E2SE.perc=abs((temp_bg$E.af-temp_bg$SE.af)/temp_bg$E.af),#perc diff for SE vs E
  		bg.E2SP.perc=abs((temp_bg$E.af-temp_bg$SP.af)/temp_bg$E.af),#perc diff for SP vs E
  		bg.E2PA.perc=abs((temp_bg$E.af-temp_bg$PA.af)/temp_bg$E.af))#perc diff for PA vs E
  	#calculate means for all bg allele differences
  	temp_means <- colMeans(temp_af)
  	names(temp_means) <- paste0(names(temp_means), ".mean")
  	#calculate variances for all bg allele differences
  	temp_vars <- colVars(temp_af)
  	names(temp_vars) <- paste0(names(temp_vars), ".var")
  	#add row to table
  	freq_means_t1_results_EvSP <- rbind(freq_means_t1_results_EvSP,append(temp_means,temp_vars))
  }, error=function(e){})
}

#save table so we don't need to run this loop every time
colnames(freq_means_t1_results_EvSP) <- paste0(colnames(freq_means_t1_results_EvSP), ".EvSP.T1")
write.table(freq_means_t1_results_EvSP, file="rudflies_2023_redo.freq_means_t1_results_EvSP.txt",sep = "\t", quote = FALSE, row.names = F)


## The following loop will go through list of "SE vs E" top outliers, find all matching
## background background SNPs, sample 100 of them, then calculate absolute frequency
## differences for treatment contrasts of interest.

## Instead of running full loop, you can just reload the table from a previous run
freq_means_t1_results_EvPA <- read.table("rudflies_2023_redo.freq_means_t1_results_EvPA.txt", header=TRUE)

### DON'T RUN CODE BLOCK IF YOU JUST LOADED TABLE FROM PREVIOUS RUN ###
freq_means_t1_results_EvPA <- c()#create empty object
## Loop through all cols pval values for candidate genes
for(g in c(1:nrow(freq_means_t1_annot_nr_EvPA))) { #number of sites in EvPA candidate genes
  tryCatch({
  	temp_chrom <- freq_means_t1_annot_nr_EvPA[g,1] #find chrom
  	temp_e <- freq_means_t1_annot_nr_EvPA[g,3] #find E mean frequency
  	temp_se <- freq_means_t1_annot_nr_EvPA[g,5] #find SE mean frequency
  	temp_sp <- freq_means_t1_annot_nr_EvPA[g,6] #find SP mean frequency
  	temp_pa <- freq_means_t1_annot_nr_EvPA[g,7] #find PA mean frequency
  	temp_cons <- freq_means_t1_annot_nr_EvPA[g,14] #find allele consequence
  	#select bg SNPs by matching chromosome, consequence, and approximate size
  	temp_bg <- freq_means_t1_annot_nr_bg[freq_means_t1_annot_nr_bg$CHROM %in% temp_chrom[1] &
  		freq_means_t1_annot_nr_bg$Consequence %in% temp_cons &
  		freq_means_t1_annot_nr_bg$E.af>(temp_e - (0.25 * min(temp_e,1-temp_e))) &
  		freq_means_t1_annot_nr_bg$E.af<(temp_e + (0.25 * min(temp_e,1-temp_e))),]
  	#randomly sample 100 SNPs, or use all samples if less than 100 matches
  	temp_bg <- temp_bg[sample(nrow(temp_bg), min(c(nrow(temp_bg),100))), ]
  	#calculate bg allele differences
  	temp_af <- cbind(bg.E2SE=abs(temp_bg$E.af-temp_bg$SE.af), #raw diff for SE vs E
  		bg.E2SP=abs(temp_bg$E.af-temp_bg$SP.af), #raw diff for SP vs E
  		bg.E2PA=abs(temp_bg$E.af-temp_bg$PA.af), #raw diff for PA vs E
  		bg.E2SE.perc=abs((temp_bg$E.af-temp_bg$SE.af)/temp_bg$E.af),#perc diff for SE vs E
  		bg.E2SP.perc=abs((temp_bg$E.af-temp_bg$SP.af)/temp_bg$E.af),#perc diff for SP vs E
  		bg.E2PA.perc=abs((temp_bg$E.af-temp_bg$PA.af)/temp_bg$E.af))#perc diff for PA vs E
  	#calculate means for all bg allele differences
  	temp_means <- colMeans(temp_af)
  	names(temp_means) <- paste0(names(temp_means), ".mean")
  	#calculate variances for all bg allele differences
  	temp_vars <- colVars(temp_af)
  	names(temp_vars) <- paste0(names(temp_vars), ".var")
  	#add row to table
  	freq_means_t1_results_EvPA <- rbind(freq_means_t1_results_EvPA,append(temp_means,temp_vars))
  }, error=function(e){})
}

#save table so we don't need to run this loop every time
colnames(freq_means_t1_results_EvPA) <- paste0(colnames(freq_means_t1_results_EvPA), ".EvPA.T1")
write.table(freq_means_t1_results_EvPA, file="rudflies_2023_redo.freq_means_t1_results_EvPA.txt",sep = "\t", quote = FALSE, row.names = F)


### Now that we've calculated background AF difference, time to calculate the same for 
### actual spino candidates and GLM outlier SNPs, then make a master table for each.

## First, make spino candidate table
freq_means_t1_annot_nr_spino_afdiff <- cbind(freq_means_t1_annot_nr_spino, 
	cand.E2SE.spino=abs(freq_means_t1_annot_nr_spino$E.af - freq_means_t1_annot_nr_spino$SE.af),
	cand.E2SP.spino=abs(freq_means_t1_annot_nr_spino$E.af - freq_means_t1_annot_nr_spino$SP.af), cand.E2PA.spino = abs(freq_means_t1_annot_nr_spino$E.af - freq_means_t1_annot_nr_spino$PA.af),
	cand.E2SE.perc.spino=abs((freq_means_t1_annot_nr_spino$E.af - freq_means_t1_annot_nr_spino$SE.af)/freq_means_t1_annot_nr_spino$E.af),
	cand.E2SP.perc.spino=abs((freq_means_t1_annot_nr_spino$E.af - freq_means_t1_annot_nr_spino$SP.af)/freq_means_t1_annot_nr_spino$E.af),
	cand.E2PA.perc.spino=abs((freq_means_t1_annot_nr_spino$E.af - freq_means_t1_annot_nr_spino$PA.af)/freq_means_t1_annot_nr_spino$E.af),
	freq_means_t1_results)

## Second, make "E vs SE" GLM outlier table
freq_means_t1_annot_nr_EvSE_afdiff <- cbind(freq_means_t1_annot_nr_EvSE,
	cand.E2SE.EvSE=abs(freq_means_t1_annot_nr_EvSE$E.af - freq_means_t1_annot_nr_EvSE$SE.af),
	cand.E2SP.EvSE=abs(freq_means_t1_annot_nr_EvSE$E.af - freq_means_t1_annot_nr_EvSE$SP.af),
	cand.E2PA.EvSE=abs(freq_means_t1_annot_nr_EvSE$E.af - freq_means_t1_annot_nr_EvSE$PA.af),
	cand.E2SE.perc.EvSE=abs((freq_means_t1_annot_nr_EvSE$E.af - freq_means_t1_annot_nr_EvSE$SE.af)/freq_means_t1_annot_nr_EvSE$E.af),
	cand.E2SP.perc.EvSE=abs((freq_means_t1_annot_nr_EvSE$E.af - freq_means_t1_annot_nr_EvSE$SP.af)/freq_means_t1_annot_nr_EvSE$E.af),
	cand.E2PA.perc.EvSE=abs((freq_means_t1_annot_nr_EvSE$E.af - freq_means_t1_annot_nr_EvSE$PA.af)/freq_means_t1_annot_nr_EvSE$E.af), 
	freq_means_t1_results_EvSE)

## Third, make "E vs SP" GLM outlier table
freq_means_t1_annot_nr_EvSP_afdiff <- cbind(freq_means_t1_annot_nr_EvSP,
	cand.E2SE.EvSP=abs(freq_means_t1_annot_nr_EvSP$E.af - freq_means_t1_annot_nr_EvSP$SE.af),
	cand.E2SP.EvSP=abs(freq_means_t1_annot_nr_EvSP$E.af - freq_means_t1_annot_nr_EvSP$SP.af),
	cand.E2PA.EvSP=abs(freq_means_t1_annot_nr_EvSP$E.af - freq_means_t1_annot_nr_EvSP$PA.af),
	cand.E2SE.perc.EvSP=abs((freq_means_t1_annot_nr_EvSP$E.af - freq_means_t1_annot_nr_EvSP$SE.af)/freq_means_t1_annot_nr_EvSP$E.af),
	cand.E2SP.perc.EvSP=abs((freq_means_t1_annot_nr_EvSP$E.af - freq_means_t1_annot_nr_EvSP$SP.af)/freq_means_t1_annot_nr_EvSP$E.af),
	cand.E2PA.perc.EvSP=abs((freq_means_t1_annot_nr_EvSP$E.af - freq_means_t1_annot_nr_EvSP$PA.af)/freq_means_t1_annot_nr_EvSP$E.af),
	freq_means_t1_results_EvSP)

## Fourth, make "E vs PA" GLM outlier table
freq_means_t1_annot_nr_EvPA_afdiff <- cbind(freq_means_t1_annot_nr_EvPA,
	cand.E2SE.EvPA=abs(freq_means_t1_annot_nr_EvPA$E.af - freq_means_t1_annot_nr_EvPA$SE.af),
	cand.E2SP.EvPA=abs(freq_means_t1_annot_nr_EvPA$E.af - freq_means_t1_annot_nr_EvPA$SP.af),
	cand.E2PA.EvPA=abs(freq_means_t1_annot_nr_EvPA$E.af - freq_means_t1_annot_nr_EvPA$PA.af),
	cand.E2SE.perc.EvPA=abs((freq_means_t1_annot_nr_EvPA$E.af - freq_means_t1_annot_nr_EvPA$SE.af)/freq_means_t1_annot_nr_EvPA$E.af),
	cand.E2SP.perc.EvPA=abs((freq_means_t1_annot_nr_EvPA$E.af - freq_means_t1_annot_nr_EvPA$SP.af)/freq_means_t1_annot_nr_EvPA$E.af),
	cand.E2PA.perc.EvPA=abs((freq_means_t1_annot_nr_EvPA$E.af - freq_means_t1_annot_nr_EvPA$PA.af)/freq_means_t1_annot_nr_EvPA$E.af),
	freq_means_t1_results_EvPA)

## Combine all relevant AF difference stats into a master table
freq_means_t1_annot_nr_afdiff <- cbind(freq_means_t1_annot_nr_spino_afdiff[,c(22:ncol(freq_means_t1_annot_nr_spino_afdiff))], freq_means_t1_annot_nr_EvSE_afdiff[,c(22:ncol(freq_means_t1_annot_nr_EvSE_afdiff))], freq_means_t1_annot_nr_EvPA_afdiff[,c(22:ncol(freq_means_t1_annot_nr_EvPA_afdiff))])

## Reorganize master table for plots: all AF diff stats in one column plus metadata cols
## E vs SE differences only
freq_means_t1_annot_nr_E2SEdiff <- data.frame(rbind(cbind("candidate","E2SE","spino","raw",freq_means_t1_annot_nr_afdiff$cand.E2SE.spino),
	   cbind("candidate","E2SE","EvSE","raw",freq_means_t1_annot_nr_afdiff$cand.E2SE.EvSE),
	   cbind("matched","E2SE","spino","mean",freq_means_t1_annot_nr_afdiff$bg.E2SE.mean.spino),
	   cbind("matched","E2SE","EvSE","mean",freq_means_t1_annot_nr_afdiff$bg.E2SE.mean.EvSE),
	   cbind("matched","E2SE","spino","perc",freq_means_t1_annot_nr_afdiff$bg.E2SE.perc.mean.spino),
	   cbind("matched","E2SE","EvSE","perc",freq_means_t1_annot_nr_afdiff$bg.E2SE.perc.mean.EvSE)))

## Rename cols	   
colnames(freq_means_t1_annot_nr_E2SEdiff) <- c("set","contrast","category","type","af_diff")
## Set AF differences to numeric format
freq_means_t1_annot_nr_E2SEdiff[,5] <- as.numeric(freq_means_t1_annot_nr_E2SEdiff[,5])

##Calculate means of AF differences by set (cand or matched), category (spino or "E vs SE" 
## outliers), and type (raw, mean, or perc)
freq_means_t1_annot_nr_E2SEdiff %>%
  group_by(set, category, type) %>%
  summarise(avg = mean(af_diff), med = median(af_diff), stdev = sd(af_diff), .groups="keep")

#T-test for differences in absolute E vs SE AF-change in E vs SE candidates vs background
t.test(af_diff ~ set, 
	data = freq_means_t1_annot_nr_E2SEdiff[freq_means_t1_annot_nr_E2SEdiff$type!="perc" &
	freq_means_t1_annot_nr_E2SEdiff$category=="EvSE",],
	alternative = 'greater')

#T-test for differences in absolute E vs SE AF-change in spino candidates vs background
t.test(af_diff ~ set, 
	data = freq_means_t1_annot_nr_E2SEdiff[freq_means_t1_annot_nr_E2SEdiff$type!="perc" &
	freq_means_t1_annot_nr_E2SEdiff$category=="spino",],
	alternative = 'greater')

### Panel F ###
# make plot of percent AF differences by set and category
E2SEdiff <- ggplot(data = freq_means_t1_annot_nr_E2SEdiff[freq_means_t1_annot_nr_E2SEdiff$type!="perc",], aes(x=set, y=af_diff)) + 
	geom_boxplot(aes(colour=category)) + 
	scale_color_manual(values=c("#39568A","#C45C5C")) +
	ggtitle("AF differences: E vs SE") + 
	theme_classic()

## Reorganize master table for plots: all AF diff stats in one column plus metadata cols
## E vs SP differences only
freq_means_t1_annot_nr_E2SPdiff <- data.frame(rbind(cbind("candidate","E2SP","spino","raw", freq_means_t1_annot_nr_afdiff$cand.E2SP.spino),
	   cbind("candidate","E2SP","EvSE","raw", freq_means_t1_annot_nr_afdiff$cand.E2SP.EvSE),
	   cbind("matched","E2SP","spino","mean", freq_means_t1_annot_nr_afdiff$bg.E2SP.mean.spino),
	   cbind("matched","E2SP","EvSE","mean", freq_means_t1_annot_nr_afdiff$bg.E2SP.mean.EvSE),
	   cbind("matched","E2SP","spino","perc", freq_means_t1_annot_nr_afdiff$bg.E2SP.perc.mean.spino),
	   cbind("matched","E2SP","EvSE","perc", freq_means_t1_annot_nr_afdiff$bg.E2SP.perc.mean.EvSE)))
	   
## Rename cols
colnames(freq_means_t1_annot_nr_E2SPdiff) <- c("set","contrast","category","type","af_diff")
## Set AF differences to numeric format
freq_means_t1_annot_nr_E2SPdiff[,5] <- as.numeric(freq_means_t1_annot_nr_E2SPdiff[,5])

##Calculate means of AF differences by set (cand or matched), category (spino or "E vs SE" 
## outliers), and type (raw, mean, or perc)
freq_means_t1_annot_nr_E2SPdiff %>%
  group_by(set, category, type) %>%
  summarise(avg = mean(af_diff), med = median(af_diff), stdev = sd(af_diff), .groups="keep")

#T-test for differences in absolute E vs SP AF-change in E vs SE candidates vs background
t.test(af_diff ~ set, 
	data = freq_means_t1_annot_nr_E2SPdiff[freq_means_t1_annot_nr_E2SPdiff$type!="perc" &
	freq_means_t1_annot_nr_E2SPdiff$category=="EvSE",],
	 alternative = 'greater')

#T-test for differences in absolute E vs SP AF-change in spino candidates vs background
t.test(af_diff ~ set, 
	data = freq_means_t1_annot_nr_E2SPdiff[freq_means_t1_annot_nr_E2SPdiff$type!="perc" &
	freq_means_t1_annot_nr_E2SPdiff$category=="spino",],
	 alternative = 'greater')
	 
### Panel G ###
# make plot of percent AF differences by set and category
E2SPdiff <- ggplot(data = freq_means_t1_annot_nr_E2SPdiff[freq_means_t1_annot_nr_E2SPdiff$type!="perc",], aes(x=set, y=af_diff)) + 
	geom_boxplot(aes(colour=category)) + 
	scale_color_manual(values=c("#39568A","#C45C5C")) +
	ggtitle("AF differences: E vs SP") + 
	theme_classic()
	
## Reorganize master table for plots: all AF diff stats in one column plus metadata cols
## E vs PA differences only
freq_means_t1_annot_nr_E2PAdiff <- data.frame(rbind(cbind("candidate","E2PA","spino","raw", freq_means_t1_annot_nr_afdiff$cand.E2PA.spino),
	cbind("candidate","E2PA","EvSE","raw", freq_means_t1_annot_nr_afdiff$cand.E2PA.EvSE),
	cbind("matched","E2PA","spino","mean", freq_means_t1_annot_nr_afdiff$bg.E2PA.mean.spino),
	cbind("matched","E2PA","EvSE","mean", freq_means_t1_annot_nr_afdiff$bg.E2PA.mean.EvSE),
	cbind("matched","E2PA","spino","perc", freq_means_t1_annot_nr_afdiff$bg.E2PA.perc.mean.spino),
	cbind("matched","E2PA","EvSE","perc", freq_means_t1_annot_nr_afdiff$bg.E2PA.perc.mean.EvSE)))
	   
## Rename cols
colnames(freq_means_t1_annot_nr_E2PAdiff) <- c("set","contrast","category","type","af_diff")
## Set AF differences to numeric format
freq_means_t1_annot_nr_E2PAdiff[,5] <- as.numeric(freq_means_t1_annot_nr_E2PAdiff[,5])

##Calculate means of AF differences by set (cand or matched), category (spino or "E vs SE" 
## outliers), and type (raw, mean, or perc)
freq_means_t1_annot_nr_E2PAdiff %>%
  group_by(set, category, type) %>%
  summarise(avg = mean(af_diff), med = median(af_diff), stdev = sd(af_diff), .groups="keep")

#T-test for differences in absolute E vs PA AF-change in E vs S candidates vs background
t.test(af_diff ~ set, 
	data = freq_means_t1_annot_nr_E2PAdiff[freq_means_t1_annot_nr_E2PAdiff$type!="perc" &
	freq_means_t1_annot_nr_E2PAdiff$category=="EvSE",],
	 alternative = 'greater')

#T-test for differences in absolute E vs PA AF-change in spino candidates vs background
t.test(af_diff ~ set, 
	data = freq_means_t1_annot_nr_E2PAdiff[freq_means_t1_annot_nr_E2PAdiff$type!="perc" &
	freq_means_t1_annot_nr_E2PAdiff$category=="spino",],
	 alternative = 'greater')

### Unused in main doc, supplementary Fig. S4 panel B ###
# make plot of percent AF differences by set and category
E2PAdiff <- ggplot(data = freq_means_t1_annot_nr_E2PAdiff[freq_means_t1_annot_nr_E2PAdiff$type!="perc",], aes(x=set, y=af_diff)) + 
	geom_boxplot(aes(colour=category)) + 
	scale_color_manual(values=c("#39568A","#C45C5C")) +
	ggtitle("AF differences: E vs PA") + 
	theme_classic()


## Save all plots together
pdf(file = "rudflies_2023_redo.freq_means_t1_annot_nr_E2SEdiff.mean.multiboxplot.pdf", width=3, height=9)
	ggarrange(E2SEdiff, E2SPdiff, E2PAdiff,
              ncol = 1, nrow = 3)
dev.off()


###########################################################
### Non-parametric SNP rank correlations and RRHO plots ###
###########################################################

## This will allow us to visualize convergence between treatments via non-parametic
## rank correlations of different sets of GLM p-values. More specifically, want to test
## for convergence of PA with either SE or SP, relative to control E populations, at TPT1.
## This will be done by taking -log10(p) values for all SNPs the "E vs PA", "E vs SE", 
## and "E vs SP" glm contrasts, making them positive or negative according to the mean 
## change in allele frequency relative to control E populations, then running pairwise
## correlations between contrasts to check for concordance in SNP ranks. 

## These ranks will be closer to 1 if the AF changes were highly significant and changing
## in the positive direction in spinosad-exposed populations relative to controls. 
## Alternatively, ranks will be closer to N (# of SNPs) if AF changes were highly 
## significant and changing in the negative direction. 

## The plots will show areas of density in the quadrants where the SNP ranks in pairwise 
## contrasts are most clustered. If the clusters appear in the lower left or upper right ## quadrants, the contrasts are positively correlated, which suggests convergence in the
## exposed populations. If the clusters appear in the upper left or lower right 
## quadrants, the contrasts are negatively correlated. This may be biologically 
## interesting, but the patterns are discordant. Finally, if there are little or no areas
## of density in the plot, the glm contrasts are uncorrelated, suggesting a lack of 
## parallel adaptation to spinosad, relative to control E populations.

## Overall, the RRHO plots are visualization methods that complement non-parametric
## Spearman Rank correlation analysis, which we will also do in this script.


#Try this color pallette for RRHO plitting
# Define a vector of colors (using names or hex codes)
rrho_cols <- c("#E8C863","#DD9551","#D24D4D","#C03268","#684187","#34215F") #light to dark

# Create a palette function that interpolates between these colors
# The resulting 'color_palette' is a function that takes an integer 'n'
# and returns 'n' interpolated colors.
rrho_palette <- colorRampPalette(rrho_cols)

# Generate 200 colors from the continuous palette
rrho_shades <- rrho_palette(200)


## Let's calculate mean frequencies, then calculate pairwise differences
## We'll use these values to assign positive or negative signs to -log10p values
## The resulting values will reflect both the magnitude and direction of AF change
freq_diff <- cbind(EvSE.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="E"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SE"]),
	EvSP.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="E"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SP"]),
	EvPA.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="E"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="PA"]),
	SEvPA.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SE"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="PA"]),
	SPvPA.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SP"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="PA"]),
	SEvSP.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SE"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SP"]))

## For use in bedtools, we'll create a bed file from AF differences
freq_diff_bed <- cbind(haf.sites.T1filt,STOP=haf.sites.T1filt$POS+1, freq_diff)
freq_diff_bed <- freq_diff_bed[order(freq_diff_bed[,1], freq_diff_bed[,2]), ]
## save
options(scipen=20)
write.table(freq_diff_bed, file="rudflies_2023_redo.freq_diff.T1.bed",sep = "\t", quote = FALSE, row.names = F)
options(scipen=0)

## Join TPT1 GLM results with new freq_diff_bed table
freq_diff_PAvSvSEvE <- merge(contrast.PAvSvSEvE.table, freq_diff_bed[,-3], by=c("CHROM","POS"))

## And sort by locus
freq_diff_PAvSvSEvE <- freq_diff_PAvSvSEvE[order(freq_diff_PAvSvSEvE[,1], freq_diff_PAvSvSEvE[,2]), ]

## Extra filter to ensure we have no rows where p=0 of our contrasts of interest
## These 0 values will result in infinite -log10(p) values, which RRHOs can't handle
freq_diff_PAvSvSEvE <- freq_diff_PAvSvSEvE[rowProds(as.matrix(freq_diff_PAvSvSEvE[,c(3:5,8)]))!=0,]

## Now calculate AF means across treatments for additional filtering of table
haf.freq.T1filt.af_mean <- cbind(haf.sites.T1filt,af_mean=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$treat.fix=="E" | haf.meta.T1filt$treat.fix=="S" | haf.meta.T1filt$treat.fix=="PA"]))
## Merge new means with GLM/freq_diff table
freq_diff_PAvSvSEvE <- merge(freq_diff_PAvSvSEvE, haf.freq.T1filt.af_mean, by=c("CHROM","POS"))

## Filter if AF means are less than 0.15 or greater than 0.85.
## This preserves ranking patterns for biologically relevant loci while removing those 
## that appear to have major treatment differences due to REF or ALT alleles being rare.
freq_diff_PAvSvSEvE <- freq_diff_PAvSvSEvE[freq_diff_PAvSvSEvE$af_mean>0.15 & freq_diff_PAvSvSEvE$af_mean<0.85,]

## Let create locus names to cross reference rankings across contrasts
RRHO.SNPs <- paste(freq_diff_PAvSvSEvE$CHROM, freq_diff_PAvSvSEvE$POS, sep = "_")

### Create a table for each contrast to use in correlations and RRHO plotting
### Here, we finally calculate signed -log10 p-values

## EvSE
EvSE.T1.signed.logp.snp.all <- data.frame(cbind(snps=RRHO.SNPs,logp=-log10(freq_diff_PAvSvSEvE$EvSE.T1)*(abs(freq_diff_PAvSvSEvE$EvSE.T1.diff)/freq_diff_PAvSvSEvE$EvSE.T1.diff)))
# reformat
EvSE.T1.signed.logp.snp.all$logp <- as.numeric(EvSE.T1.signed.logp.snp.all$logp)

## EvSP
EvSP.T1.signed.logp.snp.all <- data.frame(cbind(snps=RRHO.SNPs,logp=-log10(freq_diff_PAvSvSEvE$EvSP.T1)*(abs(freq_diff_PAvSvSEvE$EvSP.T1.diff)/freq_diff_PAvSvSEvE$EvSP.T1.diff)))
# reformat
EvSP.T1.signed.logp.snp.all$logp <- as.numeric(EvSP.T1.signed.logp.snp.all$logp)

## SPvSE
SPvSE.T1.signed.logp.snp.all <- data.frame(cbind(snps=RRHO.SNPs,logp=-log10(freq_diff_PAvSvSEvE$SEvSP.T1)*-(abs(freq_diff_PAvSvSEvE$SEvSP.T1.diff)/freq_diff_PAvSvSEvE$SEvSP.T1.diff)))
# reformat
SPvSE.T1.signed.logp.snp.all$logp <- as.numeric(SPvSE.T1.signed.logp.snp.all$logp)

## EvPA
EvPA.T1.signed.logp.snp.all <- data.frame(cbind(snps=RRHO.SNPs,logp=-log10(freq_diff_PAvSvSEvE$EvPA.T1)*(abs(freq_diff_PAvSvSEvE$EvPA.T1.diff)/freq_diff_PAvSvSEvE$EvPA.T1.diff)))
# reformat
EvPA.T1.signed.logp.snp.all$logp <- as.numeric(EvPA.T1.signed.logp.snp.all$logp)


## Spearman rank correlations
# Non-parametric rank correlation between EvSE.T1 and EvPA.T1
cor.test(as.numeric(EvSE.T1.signed.logp.snp.all[,2]), as.numeric(EvPA.T1.signed.logp.snp.all[,2]), method="spearman")
# Non-parametric rank correlation between EvSP.T1 and EvPA.T1
cor.test(as.numeric(EvSP.T1.signed.logp.snp.all[,2]), as.numeric(EvPA.T1.signed.logp.snp.all[,2]), method="spearman")
# Non-parametric rank correlation between SPvSE.T1 and EvPA.T1
cor.test(as.numeric(SPvSE.T1.signed.logp.snp.all[,2]), as.numeric(EvPA.T1.signed.logp.snp.all[,2]), method="spearman")
# Non-parametric rank correlation between EvSE.T1 and EvSP.T1
cor.test(as.numeric(EvSE.T1.signed.logp.snp.all[,2]), as.numeric(EvSP.T1.signed.logp.snp.all[,2]), method="spearman")

### Results ####
#	Spearman's rank correlation rho
#
#data:  as.numeric(EvSE.T1.signed.logp.snp.all[, 2]) and #as.numeric(EvPA.T1.signed.logp.snp.all[, 2])
#S = 9.0593e+16, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.5392372 
#
#Warning message:
#In cor.test.default(as.numeric(EvSE.T1.signed.logp.snp.all[, 2]),  :
#  Cannot compute exact p-value with ties
#
#	Spearman's rank correlation rho
#
#data:  as.numeric(EvSP.T1.signed.logp.snp.all[, 2]) and #as.numeric(EvPA.T1.signed.logp.snp.all[, 2])
#S = 1.7059e+17, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.1323443 
#
#Warning message:
#In cor.test.default(as.numeric(EvSP.T1.signed.logp.snp.all[, 2]),  :
#  Cannot compute exact p-value with ties
#
#	Spearman's rank correlation rho
#
#data:  as.numeric(SPvSE.T1.signed.logp.snp.all[, 2]) and #as.numeric(EvPA.T1.signed.logp.snp.all[, 2])
#S = 1.294e+17, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.3418404 
#
#Warning message:
#In cor.test.default(as.numeric(SPvSE.T1.signed.logp.snp.all[, 2]),  :
#  Cannot compute exact p-value with ties
#
#	Spearman's rank correlation rho
#
#data:  as.numeric(EvSE.T1.signed.logp.snp.all[, 2]) and #as.numeric(EvSP.T1.signed.logp.snp.all[, 2])
#S = 1.6134e+17, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.1793939 
#
#Warning message:
#In cor.test.default(as.numeric(EvSE.T1.signed.logp.snp.all[, 2]),  :
#  Cannot compute exact p-value with ties



### RRHO plotting ###

## Since the RRHO package is often used for gene expression studies with thousands of 
## genes, processing for millions of SNPs takes excessive computation time. For 
## simplicity, we sampled 25K SNPs multiple times and checked whether each RRHO arrived at 
## a consensus plot. For the record, they did. 

## Randomly sample master table
freq_diff_PAvSvSEvE_samp <- freq_diff_PAvSvSEvE[sample(nrow(freq_diff_PAvSvSEvE),25000), ]

## Let create locus names to cross reference rankings across contrasts
RRHO.SNPs.samp <- paste(freq_diff_PAvSvSEvE_samp$CHROM, freq_diff_PAvSvSEvE_samp$POS, sep = "_")

### Create a table for each contrast to use in correlations and RRHO plotting
### Here, we finally calculate signed -log10 p-values

## EvSE
EvSE.T1.signed.logp.snp <- data.frame(cbind(snps=RRHO.SNPs.samp,logp=-log10(freq_diff_PAvSvSEvE_samp$EvSE.T1)*(abs(freq_diff_PAvSvSEvE_samp$EvSE.T1.diff)/freq_diff_PAvSvSEvE_samp$EvSE.T1.diff)))
# reformat
EvSE.T1.signed.logp.snp$logp <- as.numeric(EvSE.T1.signed.logp.snp$logp)

## EvSP
EvSP.T1.signed.logp.snp <- data.frame(cbind(snps=RRHO.SNPs.samp,logp=-log10(freq_diff_PAvSvSEvE_samp$EvSP.T1)*(abs(freq_diff_PAvSvSEvE_samp$EvSP.T1.diff)/freq_diff_PAvSvSEvE_samp$EvSP.T1.diff)))
# reformat
EvSP.T1.signed.logp.snp$logp <- as.numeric(EvSP.T1.signed.logp.snp$logp)

## SPvSE
SPvSE.T1.signed.logp.snp <- data.frame(cbind(snps=RRHO.SNPs.samp,logp=-log10(freq_diff_PAvSvSEvE_samp$SEvSP.T1)*-(abs(freq_diff_PAvSvSEvE_samp$SEvSP.T1.diff)/freq_diff_PAvSvSEvE_samp$SEvSP.T1.diff)))
# reformat
SPvSE.T1.signed.logp.snp$logp <- as.numeric(SPvSE.T1.signed.logp.snp$logp)

## EvPA
EvPA.T1.signed.logp.snp <- data.frame(cbind(snps=RRHO.SNPs.samp,logp=-log10(freq_diff_PAvSvSEvE_samp$EvPA.T1)*(abs(freq_diff_PAvSvSEvE_samp$EvPA.T1.diff)/freq_diff_PAvSvSEvE_samp$EvPA.T1.diff)))
# reformat
EvPA.T1.signed.logp.snp$logp <- as.numeric(EvPA.T1.signed.logp.snp$logp)


### Finally, we can make RRHOs ###
## On first pass, we found the highest concordance when comparing "EvSE" and "EvPA",
## therefore we create this RRHO object first and normalize the heatmap color scale 
## maximums in all other RRHOs realative to this one. We also set the color gradient 
## to a custom one (rrho_shades).

### Figure 3 panel H
## Create RRHO object for comparison of "EvSE" and "EvPA" GLM contrasts
rrho2.EvSE.T1.EvPA.T1.signed.logp.25Ksnp <-  RRHO2_initialize(EvSE.T1.signed.logp.snp, EvPA.T1.signed.logp.snp, labels = c("EvSE.T1", "EvPA.T1"), log10.ind=TRUE)
#save plot
pdf(file = "rudflies_2023_redo.EvSE.T1.EvPA.T1.signed.logp.25Ksnp.RRHO2.pdf", width=6, height=5.5)
	RRHO2_heatmap(rrho2.EvSE.T1.EvPA.T1.signed.logp.25Ksnp, maximum = max(na.omit(rrho2.EvSE.T1.EvPA.T1.signed.logp.25Ksnp$hypermat[is.finite(rrho2.EvSE.T1.EvPA.T1.signed.logp.25Ksnp$hypermat)])), colorGradient = rrho_shades)
dev.off()

### Figure 3 panel I
## Create RRHO object for comparison of "EvSP" and "EvPA" GLM contrasts
rrho2.EvSP.T1.EvPA.T1.signed.logp.25Ksnp <-  RRHO2_initialize(EvSP.T1.signed.logp.snp, EvPA.T1.signed.logp.snp, labels = c("EvSP.T1", "EvPA.T1"), log10.ind=TRUE)
#save plot
pdf(file = "rudflies_2023_redo.EvSP.T1.EvPA.T1.signed.logp.25Ksnp.RRHO2.pdf", width=6, height=5.5)
	RRHO2_heatmap(rrho2.EvSP.T1.EvPA.T1.signed.logp.25Ksnp, maximum = max(na.omit(rrho2.EvSE.T1.EvPA.T1.signed.logp.25Ksnp$hypermat[is.finite(rrho2.EvSE.T1.EvPA.T1.signed.logp.25Ksnp$hypermat)])), colorGradient = rrho_shades)
dev.off()

## Create RRHO object for comparison of "SPvSE" and "EvPA" GLM contrasts
rrho2.SPvSE.T1.EvPA.T1.signed.logp.25Ksnp <-  RRHO2_initialize(SPvSE.T1.signed.logp.snp, EvPA.T1.signed.logp.snp, labels = c("SPvSE.T1", "EvPA.T1"), log10.ind=TRUE)
#save plot
pdf(file = "rudflies_2023_redo.SPvSE.T1.EvPA.T1.signed.logp.25Ksnp.RRHO2.pdf", width=6, height=5.5)
	RRHO2_heatmap(rrho2.SPvSE.T1.EvPA.T1.signed.logp.25Ksnp, maximum = max(na.omit(rrho2.EvSE.T1.EvPA.T1.signed.logp.25Ksnp$hypermat[is.finite(rrho2.EvSE.T1.EvPA.T1.signed.logp.25Ksnp$hypermat)])), colorGradient = rrho_shades)
dev.off()

### Figure S6 panel D
## Create RRHO object for comparison of "EvSE" and "EvSP" GLM contrasts
rrho2.EvSE.T1.EvSP.T1.signed.logp.25Ksnp <-  RRHO2_initialize(EvSE.T1.signed.logp.snp, EvSP.T1.signed.logp.snp, labels = c("EvSE.T1", "EvSP.T1"), log10.ind=TRUE)
#save plot
pdf(file = "rudflies_2023_redo.EvSE.T1.EvSP.T1.signed.logp.25Ksnp.RRHO2.pdf", width=6, height=5.5)
	RRHO2_heatmap(rrho2.EvSE.T1.EvSP.T1.signed.logp.25Ksnp, maximum = max(na.omit(rrho2.EvSE.T1.EvPA.T1.signed.logp.25Ksnp$hypermat[is.finite(rrho2.EvSE.T1.EvPA.T1.signed.logp.25Ksnp$hypermat)])), colorGradient = rrho_shades)
dev.off()



################
### Figure 4 ###
################

## For improved visualization, we aimed to mask genomic regions in Manhattan plots that
## do not appear to be convergent selective targets across spinosad-exposed treatments 
## (S and PA), relative to E control populations. To do so, we leveraged our experimental
## design to mask -log10(p) values in target GLM contrasts ("S vs E" and "PA vs E") 
## by subracting maximum-normalized values from the "PA vs S" contrast.

## The concept relies on hypothesis that PA and S genomes show major differences 
## since PA populations were split off from the E control populations prior to S and
## were exposed to spinosad selection almost immediately. Clustering analyses
## and Manhattan plots all confirm that PA is more divergent from both E and S pops
## than E and S populations are from each other. 

## Due to strong selection from spinosad, PA and S populations show putative strong 
## sweeps in similar regions relative to control populations, however genomic draft 
## obscures the actual loci subject to selection. PA and S pops should be most similar 
## at adaptive loci though while drift and draft are likely to have more disparate 
## impacts in non-target loci, resulting in more variable signals.  

## Using this information, we subtract the maximum-normalized "PA vs S" -log10(p) values
## from "S vs E" and "PA vs E" values. In effect, convergent adaptive loci should have 
## low -log10(p) values in "PA vs S" and high -log10(p) values in contrasts that include E 
## controls, therefore those regions will not be masked by subtraction and will remain 
## peaks. Conversely, regions where -log10(p) values are moderate to high in 
## "PA vs S" are not likely subject to convergent parallel adaptation, and subtraction 
## of these values will mask these non-target regions in "S vs E" and "PA vs E" GLM 
## contrasts. 



### 3-way GLM comparison
## To look for top empirical candidate loci subject to convergent selection between
## PA and S populations, relative to E, we looked for parallel adaptation in SE and PA
## in TPT1 and SP and PA in TPT4

## Join all relevant GLM results into a single filtering table
glm.all.rolwin20.adaptive.fdr <- na.omit(cbind(glm.all.rolwin20[,c(1:2)], glm.all.rolwin20$PAvE.T1.fdr.rolwin20, glm.all.rolwin20$PAvE.T4.fdr.rolwin20, glm.all.rolwin20$SvE.T4.fdr.rolwin20, glm.all.rolwin20$EvSE.T1.fdr.rolwin20))

## For filtering, calculate row maximums FDR to ensure all contrasts meet threshold
glm.all.rolwin20.adaptive.fdr$fdr_max <- rowMaxs(as.matrix(glm.all.rolwin20.adaptive.fdr[,c(3:6)]))

## Sort by maximum FDR
glm.all.rolwin20.adaptive.fdr <- glm.all.rolwin20.adaptive.fdr[order(glm.all.rolwin20.adaptive.fdr[,7]),]

## Look for loci where all contrasts of interest are significant at FDR < 0.01, control 
## allele frequencies at TPT1 are between 0.1 and 0.9, and AF difference between PA and E
## at TPT1 is greater than 0.1.
dim(merge(glm.all.rolwin20.adaptive.fdr[glm.all.rolwin20.adaptive.fdr$fdr_max < 0.01,c(1:2)], freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by=c("CHROM","POS")))

## Look for loci where all contrasts of interest are significant at FDR < 0.05, control 
## allele frequencies at TPT1 are between 0.1 and 0.9, and AF difference between PA and E
## at TPT1 is greater than 0.1.
dim(merge(glm.all.rolwin20.adaptive.fdr[glm.all.rolwin20.adaptive.fdr$fdr_max < 0.01,c(1:2)], freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by=c("CHROM","POS")))

## Nothing, so what locus has the lowest FDR across all contrasts of interest?
glm.all.rolwin20.adaptive.fdr[glm.all.rolwin20.adaptive.fdr[,7] == min(glm.all.rolwin20.adaptive.fdr[,7]),]

## Since there's no apparent overlap in SE SP and PA at FDR<0.01 or FDR<0.05, we'll find 
## the top genes (FDR < 0.1) that trend towards significance at a relaxed threshold.
glm.all.rolwin20.adaptive.fdr10 <- merge(merge(glm.all.rolwin20.adaptive.fdr[glm.all.rolwin20.adaptive.fdr$fdr_max < 0.1,c(1:2)], freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by=c("CHROM","POS")), vep, by=c("CHROM","POS"))

## Save gene list for external GO enrichment analysis in BiNGO
write.table(unique(glm.all.rolwin20.adaptive.fdr10$Gene), file="glm.all.rolwin20.adaptive.fdr10.txt", sep = "\t", quote = FALSE, row.names = F)

## Save gene list w/ FlyCADD score > 0.6 for external GO enrichment analysis in BiNGO
write.table(unique(glm.all.rolwin20.adaptive.fdr10[glm.all.rolwin20.adaptive.fdr10$FLYCADD > 0.6,]$Gene), file="glm.all.rolwin20.adaptive.fdr10.flycadd60.txt", sep = "\t", quote = FALSE, row.names = F)

### TPT1 SE and PA convergent candidates ###

## Retrieve "TPT1 convergent candidate loci" using the following criteria:
#	TPT1 SE vs E GLM contrast significant at FDR < 0.05 (1st target threshold)
#	TPT1 PA vs E GLM contrast significant at FDR < 0.01 (2nd target threshold)
#	TPT4 PA vs E GLM contrast significant at FDR < 0.01 (ensures PA & E remain divergent)
#	TPT1 E allele frequency < 0.9 (not rare alleles with possible inflated AF diff)
#	TPT1 E allele frequency > 0.1 (not rare alleles with possible inflated AF diff)
#	TPT1 absolute AF difference between PA and E > 0.1 (sizable AF diff between PA & E)
glm.all.rolwin20.TPT1.convergent_candidates <- unique(merge(unique(glm.all.rolwin20.annot[glm.all.rolwin20.annot$PAvE.T1.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$PAvE.T4.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$EvSE.T1.fdr.rolwin20<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by=c("CHROM","POS")))

## Save gene list for external GO enrichment analysis in BiNGO
write.table(unique(glm.all.rolwin20.TPT1.convergent_candidates$Gene), file="glm.all.rolwin20.TPT1.convergent_candidate_genes.txt", sep = "\t", quote = FALSE, row.names = F)

## Retrieve "Functional TPT1 candidate loci": 
# 	Same filters as above plus FlyCADD score > 0.6
#	FlyCADD details available at https://github.com/JuliaBeets/FlyCADD
glm.all.rolwin20.TPT1.functional_convergent_candidates <- unique(merge(merge(unique(glm.all.rolwin20.annot[glm.all.rolwin20.annot$PAvE.T1.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$PAvE.T4.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$EvSE.T1.fdr.rolwin20<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by=c("CHROM","POS")), vep[vep$FLYCADD > 0.6,], by=c("CHROM","POS")))

## Save gene list w/ FlyCADD score > 0.6 for external GO enrichment analysis in BiNGO
write.table(unique(glm.all.rolwin20.TPT1.functional_convergent_candidates$Gene), file="glm.all.rolwin20.TPT1.functional_convergent_candidate_genes.txt", sep = "\t", quote = FALSE, row.names = F)


### TPT4 SP and PA convergent candidates ###

## Retrieve "TPT4 convergent candidate loci" using the following criteria:
#	TPT4 SP vs E GLM contrast significant at FDR < 0.05 (1st target threshold)
#	TPT4 PA vs E GLM contrast significant at FDR < 0.01 (2nd target threshold)
#	TPT1 PA vs E GLM contrast significant at FDR < 0.01 (ensures PA & E start divergent)
#	TPT1 E allele frequency < 0.9 (not rare alleles with possible inflated AF diff)
#	TPT1 E allele frequency > 0.1 (not rare alleles with possible inflated AF diff)
#	TPT1 absolute AF difference between PA and E > 0.1 (sizable AF diff between PA & E)
glm.all.rolwin20.TPT4.convergent_candidates <- unique(merge(unique(glm.all.rolwin20.annot[glm.all.rolwin20.annot$PAvE.T1.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$PAvE.T4.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$SvE.T4.fdr.rolwin20<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by=c("CHROM","POS")))

## Save gene list for external GO enrichment analysis in BiNGO
write.table(unique(glm.all.rolwin20.TPT4.convergent_candidates$Gene), file="glm.all.rolwin20.TPT4.convergent_candidate_genes.txt", sep = "\t", quote = FALSE, row.names = F)

## Retrieve "Functional TPT4 candidate loci": 
# 	Same filters as above plus FlyCADD score > 0.6
#	FlyCADD details available at https://github.com/JuliaBeets/FlyCADD
glm.all.rolwin20.TPT4.functional_convergent_candidates <- unique(merge(merge(unique(glm.all.rolwin20.annot[glm.all.rolwin20.annot$PAvE.T1.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$PAvE.T4.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$SvE.T4.fdr.rolwin20<0.01, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by=c("CHROM","POS")), vep[vep$FLYCADD > 0.6,], by=c("CHROM","POS")))

## Save gene list w/ FlyCADD score > 0.6 for external GO enrichment analysis in BiNGO
write.table(unique(glm.all.rolwin20.TPT4.functional_convergent_candidates$Gene), file="glm.all.rolwin20.TPT4.functional_convergent_candidate_genes.txt", sep = "\t", quote = FALSE, row.names = F)


### Let's build residual Manhattan plots and overlay convergent candidates ###

### TPT1 convergent candidates ###

## Try a subtraction method to calculate residuals, in reference to PA vs E, to look 
## for loci where SE is approaching PA by T1.

## This is a really long R command, but basically we're creating a table with locus info,
## the residual values, the input -log10(p) values, and the corresponding FDR values. 

## To calculate residuals, we subtract the maximum-normalized -log10(p) values for 
## PA vs SE contrasts from PA vs E -log10(p) values. We use the following formula: 
#	residual = An - (Bn / Bmax * Amax)
#		An = TPT1 PA vs E -log10(p) value for locus n
#		Amax = maximum A value across all loci
#		Bn = TPT1 PA vs SE -log10(p) value for locus n
#		Bmax = maximum A value across all loci
glm.all.rolwin20.PAvE.logdiffT1 <- as.data.frame(na.omit(cbind(CHROM=glm.all.rolwin20$CHROM, POS=glm.all.rolwin20$POS, EvPA.T1minusPAvSE.T1=(as.numeric(glm.all.rolwin20$EvPA.T1.logp.rolwin20)-(as.numeric(glm.all.rolwin20$PAvSE.T1.logp.rolwin20)/max(na.omit(as.numeric(glm.all.rolwin20$PAvSE.T1.logp.rolwin20)))*max(na.omit(as.numeric(glm.all.rolwin20$EvPA.T1.logp.rolwin20))))),EvPA.T1.fdr.rolwin20=glm.all.rolwin20$EvPA.T1.fdr.rolwin20,EvPA.T1.logp.rolwin20=glm.all.rolwin20$EvPA.T1.logp.rolwin20,PAvSE.T1.fdr.rolwin20=glm.all.rolwin20$PAvSE.T1.fdr.rolwin20,PAvSE.T1.logp.rolwin20=glm.all.rolwin20$PAvSE.T1.logp.rolwin20)))

## Create new residual and z-score columns for exploratory analyses
glm.all.rolwin20.PAvE.logdiffT1$residual <- as.numeric(glm.all.rolwin20.PAvE.logdiffT1[,3])
glm.all.rolwin20.PAvE.logdiffT1$z <- ave(as.numeric(glm.all.rolwin20.PAvE.logdiffT1[,3]), glm.all.rolwin20.PAvE.logdiffT1[,1], FUN=scale)

## Exploratory plot to visualize spread of z-scores. 
pdf(file = "rudflies_2023_redo.PAvE.T1minusPAvSE.T1.rolwin20.zscore_hist.pdf")
	hist(glm.all.rolwin20.PAvE.logdiffT1$z)
dev.off()


## According to the motivating hypothesis, PA and S populations are both more 
## differentiated from E than PA and S are from each other at convergent adaptive loci. 
## These loci of interest will have high positive residual scores. Some residual scores 
## end up negative though, and by definition these would show higher divergence between 
## S vs P than either have from E control pops. Therefore, these loci are considered 
## non-convergent (divergent) between PA and S, and the values are floored at 0 for 
## plotting purposes. 
glm.all.rolwin20.PAvE.logdiffT1[glm.all.rolwin20.PAvE.logdiffT1[,3] < 0,3] <- 0

## Format data types for plotting
glm.all.rolwin20.PAvE.logdiffT1[,3] <- as.numeric(glm.all.rolwin20.PAvE.logdiffT1[,3])
glm.all.rolwin20.PAvE.logdiffT1[,2] <- as.integer(glm.all.rolwin20.PAvE.logdiffT1[,2])
glm.all.rolwin20.PAvE.logdiffT1[,8] <- as.numeric(glm.all.rolwin20.PAvE.logdiffT1[,8])

## Now build plot
manh.EvPA.T1minusPAvSE.T1.rolwin20 <- ggplot(glm.all.rolwin20.PAvE.logdiffT1, aes(POS, EvPA.T1minusPAvSE.T1)) + geom_line(alpha = 0.5, colour = "grey") +  
	## Add threshold that all sig. SNPs (FDR<0.05) fall above, not all SNPs over line significant though
	geom_hline(yintercept = min(na.omit(glm.all.rolwin20[glm.all.rolwin20$EvPA.T1.fdr.rolwin20<0.05,]$EvPA.T1.logp.rolwin20)), color = "black",linetype = "dashed",linewidth = .25) +
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	geom_hline(yintercept = min(na.omit(glm.all.rolwin20[glm.all.rolwin20$EvPA.T1.fdr.rolwin20<0.01,]$EvPA.T1.logp.rolwin20)), color = "black",linetype = "dashed",linewidth = .25) +
	## Highlight TPT4 convergent candidate loci, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin20.annot[glm.all.rolwin20.annot$PAvE.T1.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$PAvE.T4.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$SvE.T4.fdr.rolwin20<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")),glm.all.rolwin20.PAvE.logdiffT1, by = c("CHROM","POS")),aes(POS, as.numeric(EvPA.T1minusPAvSE.T1)), color = "#97AA71", size = 1) +
	## Highlight TPT4 convergent candidate loci w/ FlyCADD > 0.6, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin20.annot[glm.all.rolwin20.annot$PAvE.T1.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$PAvE.T4.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$SvE.T4.fdr.rolwin20<0.05 & glm.all.rolwin20.annot$FLYCADD > 0.6, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")),glm.all.rolwin20.PAvE.logdiffT1, by = c("CHROM","POS")),aes(POS, as.numeric(EvPA.T1minusPAvSE.T1)), color = "#C95B5B", size = 1) +
	## Highlight TPT1 convergent candidate loci, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin20.annot[glm.all.rolwin20.annot$PAvE.T1.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$PAvE.T4.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$EvSE.T1.fdr.rolwin20<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin20.PAvE.logdiffT1, by = c("CHROM","POS")),aes(POS, as.numeric(EvPA.T1minusPAvSE.T1)), color = "#78B0BF", size = 1) +
	## Highlight TPT4 convergent candidate loci w/ FlyCADD > 0.6, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin20.annot[glm.all.rolwin20.annot$PAvE.T1.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$PAvE.T4.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$EvSE.T1.fdr.rolwin20<0.05 & glm.all.rolwin20.annot$FLYCADD > 0.6, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin20.PAvE.logdiffT1, by = c("CHROM","POS")),aes(POS, as.numeric(EvPA.T1minusPAvSE.T1)), color = "#E8B53E", size = 1) +
	## Highlight top convergent candidate loci, filters described above
	geom_point(data = merge(merge(glm.all.rolwin20.adaptive.fdr[glm.all.rolwin20.adaptive.fdr$fdr_max < 0.1,c(1:2)], freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin20.PAvE.logdiffT1, by = c("CHROM","POS")),aes(POS, as.numeric(EvPA.T1minusPAvSE.T1)), color = "#4C4C4C", size = 1) +
  	#visual options
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, round(max(as.numeric(na.omit(glm.all.rolwin20.PAvE.logdiffT1[,3])))))) +
  	scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col = "candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("E vs PA minus PA vs SE (TPT1): 20-SNP sliding window")


## Try a subtraction method to calculate residuals, in reference to SE vs E, to look 
## for loci where SE is approaching PA by T1.

## This is a really long R command, but basically we're creating a table with locus info,
## the residual values, the input -log10(p) values, and the corresponding FDR values. 

## To calculate residuals, we subtract the maximum-normalized -log10(p) values for 
## PA vs SE contrasts from SE vs E -log10(p) values. We use the following formula: 
#	residual = An - (Bn / Bmax * Amax)
#		An = TPT1 SE vs E -log10(p) value for locus n
#		Amax = maximum A value across all loci
#		Bn = TPT1 PA vs SE -log10(p) value for locus n
#		Bmax = maximum A value across all loci
glm.all.rolwin20.SEvE.logdiffT1 <- as.data.frame(na.omit(cbind(CHROM = glm.all.rolwin20$CHROM, POS = glm.all.rolwin20$POS, EvSE.T1minusPAvSE.T1 = (as.numeric(glm.all.rolwin20$EvSE.T1.logp.rolwin20) - (as.numeric(glm.all.rolwin20$PAvSE.T1.logp.rolwin20)/max(na.omit(as.numeric(glm.all.rolwin20$PAvSE.T1.logp.rolwin20)))*max(na.omit(as.numeric(glm.all.rolwin20$EvSE.T1.logp.rolwin20))))),EvSE.T1.fdr.rolwin20 = glm.all.rolwin20$EvSE.T1.fdr.rolwin20,EvSE.T1.logp.rolwin20 = glm.all.rolwin20$EvSE.T1.logp.rolwin20,PAvSE.T1.fdr.rolwin20 = glm.all.rolwin20$PAvSE.T1.fdr.rolwin20,PAvSE.T1.logp.rolwin20 = glm.all.rolwin20$PAvSE.T1.logp.rolwin20)))

## Create new residual and z-score columns for exploratory analyses
glm.all.rolwin20.SEvE.logdiffT1$residual <- as.numeric(glm.all.rolwin20.SEvE.logdiffT1[,3])
glm.all.rolwin20.SEvE.logdiffT1$z <- ave(as.numeric(glm.all.rolwin20.SEvE.logdiffT1[,3]), glm.all.rolwin20.SEvE.logdiffT1[,1], FUN = scale)

## Exploratory plot to visualize spread of z-scores. 
pdf(file = "rudflies_2023_redo.SEvE.T1minusPAvSE.T1.rolwin20.zscore_hist.pdf")
	hist(glm.all.rolwin20.SEvE.logdiffT1$z)
dev.off()

## According to the motivating hypothesis, PA and S populations are both more 
## differentiated from E than PA and S are from each other at convergent adaptive loci. 
## These loci of interest will have high positive residual scores. Some residual scores 
## end up negative though, and by definition these would show higher divergence between 
## S vs P than either have from E control pops. Therefore, these loci are considered 
## non-convergent (divergent) between PA and S, and the values are floored at 0 for 
## plotting purposes. 
glm.all.rolwin20.SEvE.logdiffT1[glm.all.rolwin20.SEvE.logdiffT1[,3] < 0,3] <- 0

## Format data types for plotting
glm.all.rolwin20.SEvE.logdiffT1[,3] <- as.numeric(glm.all.rolwin20.SEvE.logdiffT1[,3])
glm.all.rolwin20.SEvE.logdiffT1[,2] <- as.integer(glm.all.rolwin20.SEvE.logdiffT1[,2])
glm.all.rolwin20.SEvE.logdiffT1[,8] <- as.numeric(glm.all.rolwin20.SEvE.logdiffT1[,8])

## Now build plot
manh.EvSE.T1minusPAvSE.T1.rolwin20 <- ggplot(glm.all.rolwin20.SEvE.logdiffT1, aes(POS, EvSE.T1minusPAvSE.T1)) + geom_line(alpha = 0.5, colour = "grey") +  
	## Add threshold that all sig. SNPs (FDR<0.05) fall above, not all SNPs over line significant though
	geom_hline(yintercept = min(na.omit(glm.all.rolwin20[glm.all.rolwin20$EvSE.T1.fdr.rolwin20<0.05,]$EvSE.T1.logp.rolwin20)), color = "black",linetype = "dashed",linewidth = .25) +
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	geom_hline(yintercept = min(na.omit(glm.all.rolwin20[glm.all.rolwin20$EvSE.T1.fdr.rolwin20<0.01,]$EvSE.T1.logp.rolwin20)), color = "black",linetype = "dashed",linewidth = .25) +
	## Highlight TPT4 convergent candidate loci, filters described above
geom_point(data = merge(merge(unique(glm.all.rolwin20.annot[glm.all.rolwin20.annot$PAvE.T1.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$PAvE.T4.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$SEvE.T4.fdr.rolwin20<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")),glm.all.rolwin20.SEvE.logdiffT1, by = c("CHROM","POS")),aes(POS, as.numeric(EvSE.T1minusPAvSE.T1)), color = "#97AA71", size = 1) +
	## Highlight TPT4 convergent candidate loci w/ FlyCADD > 0.6, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin20.annot[glm.all.rolwin20.annot$PAvE.T1.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$PAvE.T4.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$SEvE.T4.fdr.rolwin20<0.05 & glm.all.rolwin20.annot$FLYCADD > 0.6, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")),glm.all.rolwin20.SEvE.logdiffT1, by = c("CHROM","POS")),aes(POS, as.numeric(EvSE.T1minusPAvSE.T1)), color = "#C95B5B", size = 1) +
	## Highlight TPT1 convergent candidate loci, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin20.annot[glm.all.rolwin20.annot$PAvE.T1.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$PAvE.T4.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$EvSE.T1.fdr.rolwin20<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin20.SEvE.logdiffT1, by = c("CHROM","POS")),aes(POS, as.numeric(EvSE.T1minusPAvSE.T1)), color = "#78B0BF", size = 1) +
	## Highlight TPT4 convergent candidate loci w/ FlyCADD > 0.6, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin20.annot[glm.all.rolwin20.annot$PAvE.T1.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$PAvE.T4.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$EvSE.T1.fdr.rolwin20<0.05 & glm.all.rolwin20.annot$FLYCADD > 0.6, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin20.SEvE.logdiffT1, by = c("CHROM","POS")),aes(POS, as.numeric(EvSE.T1minusPAvSE.T1)), color = "#E8B53E", size = 1) +
	## Highlight top convergent candidate loci, filters described above
	geom_point(data = merge(merge(glm.all.rolwin20.adaptive.fdr[glm.all.rolwin20.adaptive.fdr$fdr_max < 0.1,c(1:2)], freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin20.SEvE.logdiffT1, by = c("CHROM","POS")),aes(POS, as.numeric(EvSE.T1minusPAvSE.T1)), color = "#4C4C4C", size = 1) +
	#visual options
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, round(max(as.numeric(na.omit(glm.all.rolwin20.SEvE.logdiffT1[,3])))))) +
  	scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col = "candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("E vs SE minus PA vs SE (TPT1): 20-SNP sliding window")


### TPT4 convergent candidates ###

## Try a subtraction method to calculate residuals, in reference to PA vs E, to look 
## for loci where SP is approaching PA by T4.

## This is a really long R command, but basically we're creating a table with locus info,
## the residual values, the input -log10(p) values, and the corresponding FDR values. 

## To calculate residuals, we subtract the maximum-normalized -log10(p) values for 
## PA vs SP contrasts from PA vs E -log10(p) values. We use the following formula: 
#	residual = An - (Bn / Bmax * Amax)
#		An = TPT4 PA vs E -log10(p) value for locus n
#		Amax = maximum A value across all loci
#		Bn = TPT4 PA vs SP -log10(p) value for locus n
#		Bmax = maximum A value across all loci
glm.all.rolwin20.PAvE.logdiffT4 <- as.data.frame(na.omit(cbind(CHROM = glm.all.rolwin20$CHROM, POS = glm.all.rolwin20$POS, PAvE.T4minusPAvS.T4 = (as.numeric(glm.all.rolwin20$PAvE.T4.logp.rolwin20) - (as.numeric(glm.all.rolwin20$PAvS.T4.logp.rolwin20)/max(na.omit(as.numeric(glm.all.rolwin20$PAvS.T4.logp.rolwin20)))*max(na.omit(as.numeric(glm.all.rolwin20$PAvE.T4.logp.rolwin20))))),PAvE.T4.fdr.rolwin20 = glm.all.rolwin20$PAvE.T4.fdr.rolwin20,PAvE.T4.logp.rolwin20 = glm.all.rolwin20$PAvE.T4.logp.rolwin20,PAvS.T4.fdr.rolwin20 = glm.all.rolwin20$PAvS.T4.fdr.rolwin20,PAvS.T4.logp.rolwin20 = glm.all.rolwin20$PAvS.T4.logp.rolwin20)))

## Create new residual and z-score columns for exploratory analyses
glm.all.rolwin20.PAvE.logdiffT4$residual <- as.numeric(glm.all.rolwin20.PAvE.logdiffT4[,3])
glm.all.rolwin20.PAvE.logdiffT4$z <- ave(as.numeric(glm.all.rolwin20.PAvE.logdiffT4[,3]), glm.all.rolwin20.PAvE.logdiffT4[,1], FUN = scale)

## Exploratory plot to visualize spread of z-scores. 
pdf(file = "rudflies_2023_redo.PAvE.T4minusPAvS.T4.rolwin20.zscore_hist.pdf")
	hist(glm.all.rolwin20.PAvE.logdiffT4$z)
dev.off()

## According to the motivating hypothesis, PA and S populations are both more 
## differentiated from E than PA and S are from each other at convergent adaptive loci. 
## These loci of interest will have high positive residual scores. Some residual scores 
## end up negative though, and by definition these would show higher divergence between 
## S vs P than either have from E control pops. Therefore, these loci are considered 
## non-convergent (divergent) between PA and S, and the values are floored at 0 for 
## plotting purposes. 
glm.all.rolwin20.PAvE.logdiffT4[glm.all.rolwin20.PAvE.logdiffT4[,3] < 0,3] <- 0

## Format data types for plotting
glm.all.rolwin20.PAvE.logdiffT4[,3] <- as.numeric(glm.all.rolwin20.PAvE.logdiffT4[,3])
glm.all.rolwin20.PAvE.logdiffT4[,2] <- as.integer(glm.all.rolwin20.PAvE.logdiffT4[,2])
glm.all.rolwin20.PAvE.logdiffT4[,8] <- as.numeric(glm.all.rolwin20.PAvE.logdiffT4[,8])

## Now build plot
manh.PAvE.T4minusPAvS.T4.rolwin20 <- ggplot(glm.all.rolwin20.PAvE.logdiffT4, aes(POS, PAvE.T4minusPAvS.T4)) + geom_line(alpha = 0.5, colour = "grey") +  
	## Add threshold that all sig. SNPs (FDR<0.05) fall above, not all SNPs over line significant though
	geom_hline(yintercept = min(na.omit(glm.all.rolwin20[glm.all.rolwin20$PAvE.T4.fdr.rolwin20<0.05,]$PAvE.T4.logp.rolwin20)), color = "black",linetype = "dashed",linewidth = .25) +
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	geom_hline(yintercept = min(na.omit(glm.all.rolwin20[glm.all.rolwin20$PAvE.T4.fdr.rolwin20<0.01,]$PAvE.T4.logp.rolwin20)), color = "black",linetype = "dashed",linewidth = .25) +
	## Highlight TPT4 convergent candidate loci, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin20.annot[glm.all.rolwin20.annot$PAvE.T1.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$PAvE.T4.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$SvE.T4.fdr.rolwin20<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")),glm.all.rolwin20.PAvE.logdiffT4, by = c("CHROM","POS")),aes(POS, as.numeric(PAvE.T4minusPAvS.T4)), color = "#97AA71", size = 1) +
	## Highlight TPT4 convergent candidate loci w/ FlyCADD > 0.6, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin20.annot[glm.all.rolwin20.annot$PAvE.T1.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$PAvE.T4.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$SvE.T4.fdr.rolwin20<0.05 & glm.all.rolwin20.annot$FLYCADD > 0.6, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")),glm.all.rolwin20.PAvE.logdiffT4, by = c("CHROM","POS")),aes(POS, as.numeric(PAvE.T4minusPAvS.T4)), color = "#C95B5B", size = 1) +
	## Highlight TPT1 convergent candidate loci, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin20.annot[glm.all.rolwin20.annot$PAvE.T1.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$PAvE.T4.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$EvSE.T1.fdr.rolwin20<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin20.PAvE.logdiffT4, by = c("CHROM","POS")),aes(POS, as.numeric(PAvE.T4minusPAvS.T4)), color = "#78B0BF", size = 1) +
	## Highlight TPT4 convergent candidate loci w/ FlyCADD > 0.6, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin20.annot[glm.all.rolwin20.annot$PAvE.T1.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$PAvE.T4.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$EvSE.T1.fdr.rolwin20<0.05 & glm.all.rolwin20.annot$FLYCADD > 0.6, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin20.PAvE.logdiffT4, by = c("CHROM","POS")),aes(POS, as.numeric(PAvE.T4minusPAvS.T4)), color = "#E8B53E", size = 1) +
	## Highlight top convergent candidate loci, filters described above
	geom_point(data = merge(merge(glm.all.rolwin20.adaptive.fdr[glm.all.rolwin20.adaptive.fdr$fdr_max < 0.1,c(1:2)], freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin20.PAvE.logdiffT4, by = c("CHROM","POS")),aes(POS, as.numeric(PAvE.T4minusPAvS.T4)), color = "#4C4C4C", size = 1) +
	#visual options
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, round(max(as.numeric(na.omit(glm.all.rolwin20.PAvE.logdiffT4[,3])))))) +
  	scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col = "candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("PA vs E minus PA vs S (TPT4): 20-SNP sliding window")


## Try a subtraction method to calculate residuals, in reference to SP vs E, to look 
## for loci where SP is approaching PA by T4.

## This is a really long R command, but basically we're creating a table with locus info,
## the residual values, the input -log10(p) values, and the corresponding FDR values. 

## To calculate residuals, we subtract the maximum-normalized -log10(p) values for 
## PA vs SP contrasts from SP vs E -log10(p) values. We use the following formula: 
#	residual = An - (Bn / Bmax * Amax)
#		An = TPT4 SP vs E -log10(p) value for locus n
#		Amax = maximum A value across all loci
#		Bn = TPT4 PA vs SP -log10(p) value for locus n
#		Bmax = maximum A value across all loci
glm.all.rolwin20.SPvE.logdiffT4 <- as.data.frame(na.omit(cbind(CHROM = glm.all.rolwin20$CHROM, POS = glm.all.rolwin20$POS, SvE.T4minusPAvS.T4 = (as.numeric(glm.all.rolwin20$SvE.T4.logp.rolwin20) - (as.numeric(glm.all.rolwin20$PAvS.T4.logp.rolwin20)/max(na.omit(as.numeric(glm.all.rolwin20$PAvS.T4.logp.rolwin20)))*max(na.omit(as.numeric(glm.all.rolwin20$SvE.T4.logp.rolwin20))))),SvE.T4.fdr.rolwin20 = glm.all.rolwin20$SvE.T4.fdr.rolwin20,SvE.T4.logp.rolwin20 = glm.all.rolwin20$SvE.T4.logp.rolwin20,PAvS.T4.fdr.rolwin20 = glm.all.rolwin20$PAvS.T4.fdr.rolwin20,PAvS.T4.logp.rolwin20 = glm.all.rolwin20$PAvS.T4.logp.rolwin20)))

## Create new residual and z-score columns for exploratory analyses
glm.all.rolwin20.SPvE.logdiffT4$residual <- as.numeric(glm.all.rolwin20.SPvE.logdiffT4[,3])
glm.all.rolwin20.SPvE.logdiffT4$z <- ave(as.numeric(glm.all.rolwin20.SPvE.logdiffT4[,3]), glm.all.rolwin20.SPvE.logdiffT4[,1], FUN = scale)

## Exploratory plot to visualize spread of z-scores. 
pdf(file = "rudflies_2023_redo.SPvE.T4minusPAvS.T4.rolwin20.zscore_hist.pdf")
	hist(glm.all.rolwin20.SPvE.logdiffT4$z)
dev.off()

## According to the motivating hypothesis, PA and S populations are both more 
## differentiated from E than PA and S are from each other at convergent adaptive loci. 
## These loci of interest will have high positive residual scores. Some residual scores 
## end up negative though, and by definition these would show higher divergence between 
## S vs P than either have from E control pops. Therefore, these loci are considered 
## non-convergent (divergent) between PA and S, and the values are floored at 0 for 
## plotting purposes. 
glm.all.rolwin20.SPvE.logdiffT4[glm.all.rolwin20.SPvE.logdiffT4[,3] < 0,3] <- 0

## Format data types for plotting
glm.all.rolwin20.SPvE.logdiffT4[,3] <- as.numeric(glm.all.rolwin20.SPvE.logdiffT4[,3])
glm.all.rolwin20.SPvE.logdiffT4[,2] <- as.integer(glm.all.rolwin20.SPvE.logdiffT4[,2])
glm.all.rolwin20.SPvE.logdiffT4[,8] <- as.numeric(glm.all.rolwin20.SPvE.logdiffT4[,8])

## Now build plot
manh.SPvE.T4minusPAvS.T4.rolwin20 <- ggplot(glm.all.rolwin20.SPvE.logdiffT4, aes(POS, SvE.T4minusPAvS.T4)) + geom_line(alpha = 0.5, colour = "grey") +  
	## Add threshold that all sig. SNPs (FDR<0.05) fall above, not all SNPs over line significant though
	geom_hline(yintercept = min(na.omit(glm.all.rolwin20[glm.all.rolwin20$SvE.T4.fdr.rolwin20<0.05,]$SvE.T4.logp.rolwin20)), color = "black",linetype = "dashed",linewidth = .25) +
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	geom_hline(yintercept = min(na.omit(glm.all.rolwin20[glm.all.rolwin20$SvE.T4.fdr.rolwin20<0.01,]$SvE.T4.logp.rolwin20)), color = "black",linetype = "dashed",linewidth = .25) +
	## Highlight TPT4 convergent candidate loci, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin20.annot[glm.all.rolwin20.annot$PAvE.T1.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$PAvE.T4.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$SvE.T4.fdr.rolwin20<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")),glm.all.rolwin20.SPvE.logdiffT4, by = c("CHROM","POS")),aes(POS, as.numeric(SvE.T4minusPAvS.T4)), color = "#97AA71", size = 1) +
	## Highlight TPT4 convergent candidate loci w/ FlyCADD > 0.6, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin20.annot[glm.all.rolwin20.annot$PAvE.T1.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$PAvE.T4.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$SvE.T4.fdr.rolwin20<0.05 & glm.all.rolwin20.annot$FLYCADD > 0.6, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")),glm.all.rolwin20.SPvE.logdiffT4, by = c("CHROM","POS")),aes(POS, as.numeric(SvE.T4minusPAvS.T4)), color = "#C95B5B", size = 1) +
	## Highlight TPT1 convergent candidate loci, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin20.annot[glm.all.rolwin20.annot$PAvE.T1.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$PAvE.T4.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$EvSE.T1.fdr.rolwin20<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin20.SPvE.logdiffT4, by = c("CHROM","POS")),aes(POS, as.numeric(SvE.T4minusPAvS.T4)), color = "#78B0BF", size = 1) +
	## Highlight TPT4 convergent candidate loci w/ FlyCADD > 0.6, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin20.annot[glm.all.rolwin20.annot$PAvE.T1.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$PAvE.T4.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$EvSE.T1.fdr.rolwin20<0.05 & glm.all.rolwin20.annot$FLYCADD > 0.6, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin20.SPvE.logdiffT4, by = c("CHROM","POS")),aes(POS, as.numeric(SvE.T4minusPAvS.T4)), color = "#E8B53E", size = 1) +
	## Highlight top convergent candidate loci, filters described above
	geom_point(data = merge(merge(glm.all.rolwin20.adaptive.fdr[glm.all.rolwin20.adaptive.fdr$fdr_max < 0.1,c(1:2)], freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin20.SPvE.logdiffT4, by = c("CHROM","POS")),aes(POS, as.numeric(SvE.T4minusPAvS.T4)), color = "#4C4C4C", size = 1) +
	#visual options
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, round(max(as.numeric(na.omit(glm.all.rolwin20.SPvE.logdiffT4[,3])))))) +
  	scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col = "candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("S vs E minus PA vs S (TPT4): 20-SNP sliding window")
  
  
### All TPT convergent candidates ###

## Try a subtraction method to calculate residuals, in reference to PA vs E, to look 
## for loci where S is approaching PA over all timepoints.

## This is a really long R command, but basically we're creating a table with locus info,
## the residual values, the input -log10(p) values, and the corresponding FDR values. 

## To calculate residuals, we subtract the maximum-normalized -log10(p) values for 
## PA vs SE contrasts from PA vs E -log10(p) values. We use the following formula: 
#	residual = An - (Bn / Bmax * Amax)
#		An = All TPT PA vs E -log10(p) value for locus n
#		Amax = maximum A value across all loci
#		Bn = All TPT PA vs S -log10(p) value for locus n
#		Bmax = maximum A value across all loci
glm.all.rolwin20.PAvE.logdiff <- as.data.frame(na.omit(cbind(CHROM = glm.all.rolwin20$CHROM, POS = glm.all.rolwin20$POS, PAvEminusPAvS = (as.numeric(glm.all.rolwin20$PAvE.logp.rolwin20) - (as.numeric(glm.all.rolwin20$PAvS.logp.rolwin20)/max(na.omit(as.numeric(glm.all.rolwin20$PAvS.logp.rolwin20)))*max(na.omit(as.numeric(glm.all.rolwin20$PAvE.logp.rolwin20))))),PAvE.fdr.rolwin20 = glm.all.rolwin20$PAvE.fdr.rolwin20,PAvE.logp.rolwin20 = glm.all.rolwin20$PAvE.logp.rolwin20,PAvS.fdr.rolwin20 = glm.all.rolwin20$PAvS.fdr.rolwin20,PAvS.logp.rolwin20 = glm.all.rolwin20$PAvS.logp.rolwin20)))

## Create new residual and z-score columns for exploratory analyses
glm.all.rolwin20.PAvE.logdiff$residual <- as.numeric(glm.all.rolwin20.PAvE.logdiff[,3])
glm.all.rolwin20.PAvE.logdiff$z <- ave(as.numeric(glm.all.rolwin20.PAvE.logdiff[,3]), glm.all.rolwin20.PAvE.logdiff[,1], FUN = scale)

## Exploratory plot to visualize spread of z-scores. 
pdf(file = "rudflies_2023_redo.PAvEminusPAvS.rolwin20.zscore_hist.pdf")
	hist(glm.all.rolwin20.PAvE.logdiff$z)
dev.off()

## According to the motivating hypothesis, PA and S populations are both more 
## differentiated from E than PA and S are from each other at convergent adaptive loci. 
## These loci of interest will have high positive residual scores. Some residual scores 
## end up negative though, and by definition these would show higher divergence between 
## S vs P than either have from E control pops. Therefore, these loci are considered 
## non-convergent (divergent) between PA and S, and the values are floored at 0 for 
## plotting purposes. 
glm.all.rolwin20.PAvE.logdiff[glm.all.rolwin20.PAvE.logdiff[,3] < 0,3] <- 0

## Format data types for plotting
glm.all.rolwin20.PAvE.logdiff[,3] <- as.numeric(glm.all.rolwin20.PAvE.logdiff[,3])
glm.all.rolwin20.PAvE.logdiff[,2] <- as.integer(glm.all.rolwin20.PAvE.logdiff[,2])
glm.all.rolwin20.PAvE.logdiff[,8] <- as.numeric(glm.all.rolwin20.PAvE.logdiff[,8])

## Now build plot
manh.PAvEminusPAvS.rolwin20 <- ggplot(glm.all.rolwin20.PAvE.logdiff, aes(POS, PAvEminusPAvS)) + geom_line(alpha = 0.5, colour = "grey") +  
	## Add threshold that all sig. SNPs (FDR<0.05) fall above, not all SNPs over line significant though
	geom_hline(yintercept = min(na.omit(glm.all.rolwin20[glm.all.rolwin20$PAvE.fdr.rolwin20<0.05,]$PAvE.logp.rolwin20)), color = "black",linetype = "dashed",linewidth = .25) +
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	geom_hline(yintercept = min(na.omit(glm.all.rolwin20[glm.all.rolwin20$PAvE.fdr.rolwin20<0.01,]$PAvE.logp.rolwin20)), color = "black",linetype = "dashed",linewidth = .25) +
	## Highlight TPT4 convergent candidate loci, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin20.annot[glm.all.rolwin20.annot$PAvE.T1.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$PAvE.T4.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$SvE.T4.fdr.rolwin20<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")),glm.all.rolwin20.PAvE.logdiff, by = c("CHROM","POS")),aes(POS, as.numeric(PAvEminusPAvS)), color = "#97AA71", size = 1) +
	## Highlight TPT4 convergent candidate loci w/ FlyCADD > 0.6, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin20.annot[glm.all.rolwin20.annot$PAvE.T1.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$PAvE.T4.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$SvE.T4.fdr.rolwin20<0.05 & glm.all.rolwin20.annot$FLYCADD > 0.6, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")),glm.all.rolwin20.PAvE.logdiff, by = c("CHROM","POS")),aes(POS, as.numeric(PAvEminusPAvS)), color = "#C95B5B", size = 1) +
	## Highlight TPT1 convergent candidate loci, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin20.annot[glm.all.rolwin20.annot$PAvE.T1.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$PAvE.T4.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$EvSE.T1.fdr.rolwin20<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin20.PAvE.logdiff, by = c("CHROM","POS")),aes(POS, as.numeric(PAvEminusPAvS)), color = "#78B0BF", size = 1) +
	## Highlight TPT4 convergent candidate loci w/ FlyCADD > 0.6, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin20.annot[glm.all.rolwin20.annot$PAvE.T1.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$PAvE.T4.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$EvSE.T1.fdr.rolwin20<0.05 & glm.all.rolwin20.annot$FLYCADD > 0.6, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin20.PAvE.logdiff, by = c("CHROM","POS")),aes(POS, as.numeric(PAvEminusPAvS)), color = "#E8B53E", size = 1) +
	## Highlight top convergent candidate loci, filters described above
	geom_point(data = merge(merge(glm.all.rolwin20.adaptive.fdr[glm.all.rolwin20.adaptive.fdr$fdr_max < 0.1,c(1:2)], freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin20.PAvE.logdiff, by = c("CHROM","POS")),aes(POS, as.numeric(PAvEminusPAvS)), color = "#4C4C4C", size = 1) +
  	#visual options
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, round(max(as.numeric(na.omit(glm.all.rolwin20.PAvE.logdiff[,3])))))) +
  	scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col = "candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("PA vs E minus PA vs S: 20-SNP sliding window")


## Try a subtraction method to calculate residuals, in reference to S vs E, to look 
## for loci where S is approaching PA over all timepoints.

## This is a really long R command, but basically we're creating a table with locus info,
## the residual values, the input -log10(p) values, and the corresponding FDR values. 

## To calculate residuals, we subtract the maximum-normalized -log10(p) values for 
## PA vs SE contrasts from PA vs E -log10(p) values. We use the following formula: 
#	residual = An - (Bn / Bmax * Amax)
#		An = All TPT S vs E -log10(p) value for locus n
#		Amax = maximum A value across all loci
#		Bn = All TPT PA vs S -log10(p) value for locus n
#		Bmax = maximum A value across all loci
glm.all.rolwin20.SvE.logdiff <- as.data.frame(na.omit(cbind(CHROM = glm.all.rolwin20$CHROM, POS = glm.all.rolwin20$POS, SvEminusPAvS = (as.numeric(glm.all.rolwin20$SvE.logp.rolwin20) - (as.numeric(glm.all.rolwin20$PAvS.logp.rolwin20)/max(na.omit(as.numeric(glm.all.rolwin20$PAvS.logp.rolwin20)))*max(na.omit(as.numeric(glm.all.rolwin20$SvE.logp.rolwin20))))),SvE.fdr.rolwin20 = glm.all.rolwin20$SvE.fdr.rolwin20,SvE.logp.rolwin20 = glm.all.rolwin20$SvE.logp.rolwin20,PAvS.fdr.rolwin20 = glm.all.rolwin20$PAvS.fdr.rolwin20,PAvS.logp.rolwin20 = glm.all.rolwin20$PAvS.logp.rolwin20)))

## Create new residual and z-score columns for exploratory analyses
glm.all.rolwin20.SvE.logdiff$residual <- as.numeric(glm.all.rolwin20.SvE.logdiff[,3])
glm.all.rolwin20.SvE.logdiff$z <- ave(as.numeric(glm.all.rolwin20.SvE.logdiff[,3]), glm.all.rolwin20.SvE.logdiff[,1], FUN = scale)

## Exploratory plot to visualize spread of z-scores. 
pdf(file = "rudflies_2023_redo.SvEminusPAvS.rolwin20.zscore_hist.pdf")
	hist(glm.all.rolwin20.SvE.logdiff$z)
dev.off()

## According to the motivating hypothesis, PA and S populations are both more 
## differentiated from E than PA and S are from each other at convergent adaptive loci. 
## These loci of interest will have high positive residual scores. Some residual scores 
## end up negative though, and by definition these would show higher divergence between 
## S vs P than either have from E control pops. Therefore, these loci are considered 
## non-convergent (divergent) between PA and S, and the values are floored at 0 for 
## plotting purposes. 
glm.all.rolwin20.SvE.logdiff[glm.all.rolwin20.SvE.logdiff[,3] < 0,3] <- 0

## Format data types for plotting
glm.all.rolwin20.SvE.logdiff[,3] <- as.numeric(glm.all.rolwin20.SvE.logdiff[,3])
glm.all.rolwin20.SvE.logdiff[,2] <- as.integer(glm.all.rolwin20.SvE.logdiff[,2])
glm.all.rolwin20.SvE.logdiff[,8] <- as.numeric(glm.all.rolwin20.SvE.logdiff[,8])

## Now build plot
manh.SvEminusPAvS.rolwin20 <- ggplot(glm.all.rolwin20.SvE.logdiff, aes(POS, SvEminusPAvS)) + geom_line(alpha = 0.5, colour = "grey") +  
	## Add threshold that all sig. SNPs (FDR<0.05) fall above, not all SNPs over line significant though
	geom_hline(yintercept = min(na.omit(glm.all.rolwin20[glm.all.rolwin20$SvE.fdr.rolwin20<0.05,]$SvE.logp.rolwin20)), color = "black",linetype = "dashed",linewidth = .25) +
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	geom_hline(yintercept = min(na.omit(glm.all.rolwin20[glm.all.rolwin20$SvE.fdr.rolwin20<0.01,]$SvE.logp.rolwin20)), color = "black",linetype = "dashed",linewidth = .25) +
	## Highlight TPT4 convergent candidate loci, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin20.annot[glm.all.rolwin20.annot$PAvE.T1.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$PAvE.T4.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$SvE.T4.fdr.rolwin20<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")),glm.all.rolwin20.SvE.logdiff, by = c("CHROM","POS")),aes(POS, as.numeric(SvEminusPAvS)), color = "#97AA71", size = 1) +
	## Highlight TPT4 convergent candidate loci w/ FlyCADD > 0.6, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin20.annot[glm.all.rolwin20.annot$PAvE.T1.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$PAvE.T4.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$SvE.T4.fdr.rolwin20<0.05 & glm.all.rolwin20.annot$FLYCADD > 0.6, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")),glm.all.rolwin20.SvE.logdiff, by = c("CHROM","POS")),aes(POS, as.numeric(SvEminusPAvS)), color = "#C95B5B", size = 1) +
	## Highlight TPT1 convergent candidate loci, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin20.annot[glm.all.rolwin20.annot$PAvE.T1.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$PAvE.T4.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$EvSE.T1.fdr.rolwin20<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin20.SvE.logdiff, by = c("CHROM","POS")),aes(POS, as.numeric(SvEminusPAvS)), color = "#78B0BF", size = 1) +
	## Highlight TPT4 convergent candidate loci w/ FlyCADD > 0.6, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin20.annot[glm.all.rolwin20.annot$PAvE.T1.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$PAvE.T4.fdr.rolwin20<0.01 & glm.all.rolwin20.annot$EvSE.T1.fdr.rolwin20<0.05 & glm.all.rolwin20.annot$FLYCADD > 0.6, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin20.SvE.logdiff, by = c("CHROM","POS")),aes(POS, as.numeric(SvEminusPAvS)), color = "#E8B53E", size = 1) +
	## Highlight top convergent candidate loci, filters described above
	geom_point(data = merge(merge(glm.all.rolwin20.adaptive.fdr[glm.all.rolwin20.adaptive.fdr$fdr_max < 0.1,c(1:2)], freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin20.SvE.logdiff, by = c("CHROM","POS")),aes(POS, as.numeric(SvEminusPAvS)), color = "#4C4C4C", size = 1) +
	#visual options
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, round(max(as.numeric(na.omit(glm.all.rolwin20.SvE.logdiff[,3])))))) +
  	scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col = "candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("S vs E minus PA vs S: 20-SNP sliding window")

### Plot all reciprocal residual plots, modified for Figure 4 and Figure S8 ###

## 	In the plots below the highlighted points are as follows:
##	blue - TPT1 SE and PA convergent adaptive loci
##	yellow - TPT1 SE and PA convergent adaptive loci w/ FlyCADD score > 0.06
##	green - TPT4 SP and PA convergent adaptive loci
##	brown - TPT4 SP and PA convergent adaptive loci w/ FlyCADD score > 0.06
##	black - TPT1 SE, TPT4 SP, and PA convergent adaptive loci
pdf(file = "rudflies_2023_redo.convergent_evo_residuals.rolwin20.glm.manh.pdf", width = 11.25, height = 9)
	ggarrange(manh.EvPA.T1minusPAvSE.T1.rolwin20, 
			  manh.EvSE.T1minusPAvSE.T1.rolwin20, 
			  manh.PAvE.T4minusPAvS.T4.rolwin20, 
			  manh.SPvE.T4minusPAvS.T4.rolwin20, 
			  manh.PAvEminusPAvS.rolwin20, 
			  manh.SvEminusPAvS.rolwin20, 
              ncol = 2, nrow = 3)
dev.off()





### Allele frequency plots over time for top candidates ###

## TPT1 SE & PA convergent adaptive candidates
## Merge functional candidate list with residual score tables
TPT1_top_candidates <- merge(merge(glm.all.rolwin20.TPT1.functional_convergent_candidates, glm.all.rolwin20.SEvE.logdiffT1, by=c("CHROM","POS")), glm.all.rolwin20.PAvE.logdiffT1, by=c("CHROM","POS"))

TPT1_top_candidates$flycadd_rank <- NA
TPT1_top_candidates$flycadd_rank <- rank(TPT1_top_candidates$FLYCADD)
TPT1_top_candidates$flycadd_rank <- TPT1_top_candidates$flycadd_rank/max(TPT1_top_candidates$flycadd_rank)*100
TPT1_top_candidates$SE_rank <- NA
TPT1_top_candidates$SE_rank <- rank(TPT1_top_candidates$EvSE.T1minusPAvSE.T1)
TPT1_top_candidates$SE_rank <- TPT1_top_candidates$SE_rank/max(TPT1_top_candidates$SE_rank)*100
TPT1_top_candidates$PA_rank <- NA
TPT1_top_candidates$PA_rank <- rank(TPT1_top_candidates$EvPA.T1minusPAvSE.T1)
TPT1_top_candidates$PA_rank <- TPT1_top_candidates$PA_rank/max(TPT1_top_candidates$PA_rank)*100

TPT1_top_candidates$mean_rank1 <- NA
TPT1_top_candidates$mean_rank1 <- (TPT1_top_candidates$flycadd_rank + TPT1_top_candidates$SE_rank + TPT1_top_candidates$PA_rank)/3
TPT1_top_candidates[TPT1_top_candidates$mean_rank1==max(TPT1_top_candidates$mean_rank1),]

## 3R_15541740 - top TPT1 candidate in SE, PA, and ranked flycadd scores
tempplot <- cbind(haf.meta.filt,t(na.omit(haf.freq.filt[haf.sites.filt$CHROM=="3R" & haf.sites.filt$POS=="15541740",])))
colnames(tempplot)[13] <- "af"
tempplot$treat.fix <- as.character(tempplot$treat.fix)
tempplot[tempplot$tpt=="1" & (tempplot$cage=="11" | tempplot$cage=="15" | tempplot$cage=="21" | tempplot$cage=="27" | tempplot$cage=="41" | tempplot$cage=="45"),8] <- "SE"
tempplot[tempplot$treat.fix=="S",8] <- "SP"
tempplot$treat.fix <- as.factor(tempplot$treat.fix)

tempplot.mean <- tempplot %>%
  group_by(treat.fix, tpt) %>% 
  summarise_at(vars("af"), mean)
  
tempplot.mean <- tempplot.mean[tempplot.mean$tpt != 0,]

##stip chart
pdf(file = "rudflies_2023_redo.haf.3R_15541740.geneplot.PAvSvE.strip.pdf", width=4, height=4)
	ggplot(tempplot) + geom_point(position=position_dodge(width=0.5), aes(x=tpt, y=af, color=treat.fix),
              width=0.2, height=0.1, alpha=0.4, size=1, shape=19) + theme_classic() +  ggtitle("3R:15541740") + scale_color_manual(values = sample_cols) + geom_line(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix)) + geom_point(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix), size=2.5)
dev.off()

TPT1_top_candidates$mean_rank2 <- (TPT1_top_candidates$SE_rank + TPT1_top_candidates$PA_rank)/2
TPT1_top_candidates[TPT1_top_candidates$mean_rank2==max(TPT1_top_candidates$mean_rank2),]

## 3R_15319344 - top TPT1 candidate in SE and PA (plus minimum flycadd score of 0.6)
tempplot <- cbind(haf.meta.filt,t(na.omit(haf.freq.filt[haf.sites.filt$CHROM=="3R" & haf.sites.filt$POS=="15319344",])))
colnames(tempplot)[13] <- "af"
tempplot$treat.fix <- as.character(tempplot$treat.fix)
tempplot[tempplot$tpt=="1" & (tempplot$cage=="11" | tempplot$cage=="15" | tempplot$cage=="21" | tempplot$cage=="27" | tempplot$cage=="41" | tempplot$cage=="45"),8] <- "SE"
tempplot[tempplot$treat.fix=="S",8] <- "SP"
tempplot$treat.fix <- as.factor(tempplot$treat.fix)

tempplot.mean <- tempplot %>%
  group_by(treat.fix, tpt) %>% 
  summarise_at(vars("af"), mean)
  
tempplot.mean <- tempplot.mean[tempplot.mean$tpt != 0,]

##stip chart
pdf(file = "rudflies_2023_redo.haf.3R_15319344.geneplot.PAvSvE.strip.pdf", width=4, height=4)
	ggplot(tempplot) + geom_point(position=position_dodge(width=0.5), aes(x=tpt, y=af, color=treat.fix),
              width=0.2, height=0.1, alpha=0.4, size=1, shape=19) + theme_classic() +  ggtitle("3R:15319344") + scale_color_manual(values = sample_cols) + geom_line(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix)) + geom_point(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix), size=2.5)
dev.off()

## TPT4 SP & PA convergent adaptive candidates
TPT4_top_candidates <- merge(merge(glm.all.rolwin20.TPT4.functional_convergent_candidates, glm.all.rolwin20.SPvE.logdiffT4, by=c("CHROM","POS")), glm.all.rolwin20.PAvE.logdiffT4, by=c("CHROM","POS"))

TPT4_top_candidates$flycadd_rank <- NA
TPT4_top_candidates$flycadd_rank <- rank(TPT4_top_candidates$FLYCADD)
TPT4_top_candidates$flycadd_rank <- TPT4_top_candidates$flycadd_rank/max(TPT4_top_candidates$flycadd_rank)*100

TPT4_top_candidates$SP_rank <- NA
TPT4_top_candidates$SP_rank <- rank(TPT4_top_candidates$SvE.T4minusPAvS.T4)
TPT4_top_candidates$SP_rank <- TPT4_top_candidates$SP_rank/max(TPT4_top_candidates$SP_rank)*100

TPT4_top_candidates$PA_rank <- NA
TPT4_top_candidates$PA_rank <- rank(TPT4_top_candidates$PAvE.T4minusPAvS.T4)
TPT4_top_candidates$PA_rank <- TPT4_top_candidates$PA_rank/max(TPT4_top_candidates$PA_rank)*100

TPT4_top_candidates$mean_rank1 <- NA
TPT4_top_candidates$mean_rank1 <- (TPT4_top_candidates$flycadd_rank + TPT4_top_candidates$SP_rank + TPT4_top_candidates$PA_rank)/3
TPT4_top_candidates[TPT4_top_candidates$mean_rank1==max(TPT4_top_candidates$mean_rank1),]

## 3R_7856392 - top TPT4 candidate in SP, PA, and ranked flycadd scores
tempplot <- cbind(haf.meta.filt,t(na.omit(haf.freq.filt[haf.sites.filt$CHROM=="3R" & haf.sites.filt$POS=="7856392",])))
colnames(tempplot)[13] <- "af"
tempplot$treat.fix <- as.character(tempplot$treat.fix)
tempplot[tempplot$tpt=="1" & (tempplot$cage=="11" | tempplot$cage=="15" | tempplot$cage=="21" | tempplot$cage=="27" | tempplot$cage=="41" | tempplot$cage=="45"),8] <- "SE"
tempplot[tempplot$treat.fix=="S",8] <- "SP"
tempplot$treat.fix <- as.factor(tempplot$treat.fix)

tempplot.mean <- tempplot %>%
  group_by(treat.fix, tpt) %>% 
  summarise_at(vars("af"), mean)
  
tempplot.mean <- tempplot.mean[tempplot.mean$tpt != 0,]

##stip chart
pdf(file = "rudflies_2023_redo.haf.3R_7856392.geneplot.PAvSvE.strip.pdf", width=4, height=4)
	ggplot(tempplot) + geom_point(position=position_dodge(width=0.5), aes(x=tpt, y=af, color=treat.fix),
              width=0.2, height=0.1, alpha=0.4, size=1, shape=19) + theme_classic() +  ggtitle("3R:7856392") + scale_color_manual(values = sample_cols) + geom_line(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix)) + geom_point(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix), size=2.5)
dev.off()

TPT4_top_candidates$mean_rank2 <- (TPT4_top_candidates$SP_rank + TPT4_top_candidates$PA_rank)/2
TPT4_top_candidates[TPT4_top_candidates$mean_rank2==max(TPT4_top_candidates$mean_rank2),]

## 3R_7856044 - top TPT4 candidate in SP and PA (plus minimum flycadd score of 0.6)
tempplot <- cbind(haf.meta.filt,t(na.omit(haf.freq.filt[haf.sites.filt$CHROM=="3R" & haf.sites.filt$POS=="7856044",])))
colnames(tempplot)[13] <- "af"
tempplot$treat.fix <- as.character(tempplot$treat.fix)
tempplot[tempplot$tpt=="1" & (tempplot$cage=="11" | tempplot$cage=="15" | tempplot$cage=="21" | tempplot$cage=="27" | tempplot$cage=="41" | tempplot$cage=="45"),8] <- "SE"
tempplot[tempplot$treat.fix=="S",8] <- "SP"
tempplot$treat.fix <- as.factor(tempplot$treat.fix)

tempplot.mean <- tempplot %>%
  group_by(treat.fix, tpt) %>% 
  summarise_at(vars("af"), mean)
  
tempplot.mean <- tempplot.mean[tempplot.mean$tpt != 0,]

##stip chart
pdf(file = "rudflies_2023_redo.haf.3R_7856044.geneplot.PAvSvE.strip.pdf", width=4, height=4)
	ggplot(tempplot) + geom_point(position=position_dodge(width=0.5), aes(x=tpt, y=af, color=treat.fix),
              width=0.2, height=0.1, alpha=0.4, size=1, shape=19) + theme_classic() +  ggtitle("3R:7856044") + scale_color_manual(values = sample_cols) + geom_line(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix)) + geom_point(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix), size=2.5)
dev.off()


## All TPT S & PA convergent adaptive candidates
AllTPT_top_candidates <- merge(merge(glm.all.rolwin20.adaptive.fdr10, glm.all.rolwin20.SvE.logdiff, by=c("CHROM","POS")), glm.all.rolwin20.PAvE.logdiff, by=c("CHROM","POS"))

AllTPT_top_candidates$flycadd_rank <- NA
AllTPT_top_candidates$flycadd_rank <- rank(AllTPT_top_candidates$FLYCADD)
AllTPT_top_candidates$flycadd_rank <- AllTPT_top_candidates$flycadd_rank/max(AllTPT_top_candidates$flycadd_rank)*100

AllTPT_top_candidates$S_rank <- NA
AllTPT_top_candidates$S_rank <- rank(AllTPT_top_candidates$SvEminusPAvS)
AllTPT_top_candidates$S_rank <- AllTPT_top_candidates$S_rank/max(AllTPT_top_candidates$S_rank)*100

AllTPT_top_candidates$PA_rank <- NA
AllTPT_top_candidates$PA_rank <- rank(AllTPT_top_candidates$PAvEminusPAvS)
AllTPT_top_candidates$PA_rank <- AllTPT_top_candidates$PA_rank/max(AllTPT_top_candidates$PA_rank)*100

AllTPT_top_candidates$mean_rank1 <- NA
AllTPT_top_candidates$mean_rank1 <- (AllTPT_top_candidates$flycadd_rank + AllTPT_top_candidates$S_rank + AllTPT_top_candidates$PA_rank)/3
AllTPT_top_candidates[AllTPT_top_candidates$mean_rank1==max(AllTPT_top_candidates$mean_rank1),]

## 3L_18437252 - top candidate in all TPT in S, PA, and ranked flycadd scores
tempplot <- cbind(haf.meta.filt,t(na.omit(haf.freq.filt[haf.sites.filt$CHROM=="3L" & haf.sites.filt$POS=="18437252",])))
colnames(tempplot)[13] <- "af"
tempplot$treat.fix <- as.character(tempplot$treat.fix)
tempplot[tempplot$tpt=="1" & (tempplot$cage=="11" | tempplot$cage=="15" | tempplot$cage=="21" | tempplot$cage=="27" | tempplot$cage=="41" | tempplot$cage=="45"),8] <- "SE"
tempplot[tempplot$treat.fix=="S",8] <- "SP"
tempplot$treat.fix <- as.factor(tempplot$treat.fix)

tempplot.mean <- tempplot %>%
  group_by(treat.fix, tpt) %>% 
  summarise_at(vars("af"), mean)
  
tempplot.mean <- tempplot.mean[tempplot.mean$tpt != 0,]

##stip chart
pdf(file = "rudflies_2023_redo.haf.3L_18437252.geneplot.PAvSvE.strip.pdf", width=4, height=4)
	ggplot(tempplot) + geom_point(position=position_dodge(width=0.5), aes(x=tpt, y=af, color=treat.fix),
              width=0.2, height=0.1, alpha=0.4, size=1, shape=19) + theme_classic() +  ggtitle("3L:18437252") + scale_color_manual(values = sample_cols) + geom_line(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix)) + geom_point(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix), size=2.5)
dev.off()

AllTPT_top_candidates$mean_rank2 <- (AllTPT_top_candidates$S_rank + AllTPT_top_candidates$PA_rank)/2
AllTPT_top_candidates[AllTPT_top_candidates$mean_rank2==max(AllTPT_top_candidates$mean_rank2),]

## 3R_13398258 - top candidate in all TPT in S and PA (but no flycadd score weighting)
tempplot <- cbind(haf.meta.filt,t(na.omit(haf.freq.filt[haf.sites.filt$CHROM=="3R" & haf.sites.filt$POS=="13398258",])))
colnames(tempplot)[13] <- "af"
tempplot$treat.fix <- as.character(tempplot$treat.fix)
tempplot[tempplot$tpt=="1" & (tempplot$cage=="11" | tempplot$cage=="15" | tempplot$cage=="21" | tempplot$cage=="27" | tempplot$cage=="41" | tempplot$cage=="45"),8] <- "SE"
tempplot[tempplot$treat.fix=="S",8] <- "SP"
tempplot$treat.fix <- as.factor(tempplot$treat.fix)

tempplot.mean <- tempplot %>%
  group_by(treat.fix, tpt) %>% 
  summarise_at(vars("af"), mean)
  
tempplot.mean <- tempplot.mean[tempplot.mean$tpt != 0,]

##stip chart
pdf(file = "rudflies_2023_redo.haf.3R_13398258.geneplot.PAvSvE.strip.pdf", width=4, height=4)
	ggplot(tempplot) + geom_point(position=position_dodge(width=0.5), aes(x=tpt, y=af, color=treat.fix),
              width=0.2, height=0.1, alpha=0.4, size=1, shape=19) + theme_classic() +  ggtitle("3R:13398258") + scale_color_manual(values = sample_cols) + geom_line(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix)) + geom_point(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix), size=2.5)
dev.off()



### Final note about Z-scores ###
## Z-scores were helpful in determining local outliers for exploratory purposes (not 
## reported here), but the first pass Z-scores were normalized across the genome and 
## differentiation was not uniform across different chromosomes or LD blocks. To address 
## this, we calculated Z-scores again normalizing within local LD windows. This was 
## performed outside of R using bedools implemented in the script "bedtools_LD_zscore.sh".

## Save all residual/zscore tables for external use 
options(scipen=20)
#PA TPT1 convergent candidates
write.table(cbind(glm.all.rolwin20.PAvE.logdiffT1[,c(1,2)], STOP=glm.all.rolwin20.PAvE.logdiffT1[,2]+1, glm.all.rolwin20.PAvE.logdiffT1[,c(3:ncol(glm.all.rolwin20.PAvE.logdiffT1))]), file="glm.all.rolwin20.PAvE.logdiffT1.bed", sep = "\t", quote = FALSE, row.names = F)

#SE TPT1 convergent candidates
write.table(cbind(glm.all.rolwin20.SEvE.logdiffT1[,c(1,2)], STOP=glm.all.rolwin20.SEvE.logdiffT1[,2]+1, glm.all.rolwin20.SEvE.logdiffT1[,c(3:ncol(glm.all.rolwin20.SEvE.logdiffT1))]), file="glm.all.rolwin20.SEvE.logdiffT1.bed", sep = "\t", quote = FALSE, row.names = F)

#PA TPT4 convergent candidates
write.table(cbind(glm.all.rolwin20.PAvE.logdiffT4[,c(1,2)], STOP=glm.all.rolwin20.PAvE.logdiffT4[,2]+1, glm.all.rolwin20.PAvE.logdiffT4[,c(3:ncol(glm.all.rolwin20.PAvE.logdiffT4))]), file="glm.all.rolwin20.PAvE.logdiffT4.bed", sep = "\t", quote = FALSE, row.names = F)

#SP TPT4 convergent candidates
write.table(cbind(glm.all.rolwin20.SPvE.logdiffT4[,c(1,2)], STOP=glm.all.rolwin20.SPvE.logdiffT4[,2]+1, glm.all.rolwin20.SPvE.logdiffT4[,c(3:ncol(glm.all.rolwin20.SPvE.logdiffT4))]), file="glm.all.rolwin20.SPvE.logdiffT4.bed", sep = "\t", quote = FALSE, row.names = F)

#PA All TPT convergent candidates
write.table(cbind(glm.all.rolwin20.PAvE.logdiff[,c(1,2)], STOP=glm.all.rolwin20.PAvE.logdiff[,2]+1, glm.all.rolwin20.PAvE.logdiff[,c(3:ncol(glm.all.rolwin20.PAvE.logdiff))]), file="glm.all.rolwin20.PAvE.logdiff.bed", sep = "\t", quote = FALSE, row.names = F)

#S All TPT convergent candidates
write.table(cbind(glm.all.rolwin20.SvE.logdiff[,c(1,2)], STOP=glm.all.rolwin20.SvE.logdiff[,2]+1, glm.all.rolwin20.SvE.logdiff[,c(3:ncol(glm.all.rolwin20.SvE.logdiff))]), file="glm.all.rolwin20.SvE.logdiff.bed", sep = "\t", quote = FALSE, row.names = F)
options(scipen=0)



########################################
### Supplementary figures and tables ###
########################################


### Figure S3 ###

## This is a multi-panel manhattan plot that shows all leave-one-out GLM results. These
## GLM analyses were performed while iteratively dropping each S population and running
## TPT1 filtered pairwise treatment contrasts. Due to the low N for both SE (extinct) 
## and SP (persistent) populations, we can assess whether any indivual S sample had an 
## outsized impact on results, compared to the GLM results utilizing all samples.

## First load all GLM results generated via the separate R script:
## "glm.rudflies2023.PAvSvSEvE.leave1out.r"

# Dropped cage 3 (SP)
glm_no3 <- read.table("rudflies_2023_PAvSvSEvE.wLociGLMcontrast.LOO.no3.txt", header=FALSE)
names(glm_no3) <- c("CHROM","POS","EvPA.T1_no3","EvSE.T1_no3","EvSP.T1_no3","PAvSE.T1_no3","PAvSP.T1_no3","SEvSP.T1_no3")

# Dropped cage 7 (SP)
glm_no7 <- read.table("rudflies_2023_PAvSvSEvE.wLociGLMcontrast.LOO.no7.txt", header=FALSE)
names(glm_no7) <- c("CHROM","POS","EvPA.T1_no7","EvSE.T1_no7","EvSP.T1_no7","PAvSE.T1_no7","PAvSP.T1_no7","SEvSP.T1_no7")

# Dropped cage 11 (SE)
glm_no11 <- read.table("rudflies_2023_PAvSvSEvE.wLociGLMcontrast.LOO.no11.txt", header=FALSE)
names(glm_no11) <- c("CHROM","POS","EvPA.T1_no11","EvSE.T1_no11","EvSP.T1_no11","PAvSE.T1_no11","PAvSP.T1_no11","SEvSP.T1_no11")

# Dropped cage 15 (SE)
glm_no15 <- read.table("rudflies_2023_PAvSvSEvE.wLociGLMcontrast.LOO.no15.txt", header=FALSE)
names(glm_no15) <- c("CHROM","POS","EvPA.T1_no15","EvSE.T1_no15","EvSP.T1_no15","PAvSE.T1_no15","PAvSP.T1_no15","SEvSP.T1_no15")

# Dropped cage 21 (SE)
glm_no21 <- read.table("rudflies_2023_PAvSvSEvE.wLociGLMcontrast.LOO.no21.txt", header=FALSE)
names(glm_no21) <- c("CHROM","POS","EvPA.T1_no21","EvSE.T1_no21","EvSP.T1_no21","PAvSE.T1_no21","PAvSP.T1_no21","SEvSP.T1_no21")

# Dropped cage 27 (SE)
glm_no27 <- read.table("rudflies_2023_PAvSvSEvE.wLociGLMcontrast.LOO.no27.txt", header=FALSE)
names(glm_no27) <- c("CHROM","POS","EvPA.T1_no27","EvSE.T1_no27","EvSP.T1_no27","PAvSE.T1_no27","PAvSP.T1_no27","SEvSP.T1_no27")

# Dropped cage 33 (SP)
glm_no33 <- read.table("rudflies_2023_PAvSvSEvE.wLociGLMcontrast.LOO.no33.txt", header=FALSE)
names(glm_no33) <- c("CHROM","POS","EvPA.T1_no33","EvSE.T1_no33","EvSP.T1_no33","PAvSE.T1_no33","PAvSP.T1_no33","SEvSP.T1_no33")

# Dropped cage 37 (SP)
glm_no37 <- read.table("rudflies_2023_PAvSvSEvE.wLociGLMcontrast.LOO.no37.txt", header=FALSE)
names(glm_no37) <- c("CHROM","POS","EvPA.T1_no37","EvSE.T1_no37","EvSP.T1_no37","PAvSE.T1_no37","PAvSP.T1_no37","SEvSP.T1_no37")

# Dropped cage 41 (SE)
glm_no41 <- read.table("rudflies_2023_PAvSvSEvE.wLociGLMcontrast.LOO.no41.txt", header=FALSE)
names(glm_no41) <- c("CHROM","POS","EvPA.T1_no41","EvSE.T1_no41","EvSP.T1_no41","PAvSE.T1_no41","PAvSP.T1_no41","SEvSP.T1_no41")

# Dropped cage 45 (SE)
glm_no45 <- read.table("rudflies_2023_PAvSvSEvE.wLociGLMcontrast.LOO.no45.txt", header=FALSE)
names(glm_no45) <- c("CHROM","POS","EvPA.T1_no45","EvSE.T1_no45","EvSP.T1_no45","PAvSE.T1_no45","PAvSP.T1_no45","SEvSP.T1_no45")

# No dropped cages for reference
glm_PAvSvSEvE_all <- read.table("rudflies_2023_PAvSvSEvE.wLoci.GLMcontrast.txt", header=TRUE)
glm_PAvSvSEvE_all <- na.omit(glm_PAvSvSEvE_all)
names(glm_PAvSvSEvE_all) <- c("CHROM","POS","EvPA.T1","EvSE.T1","EvSP.T1","PAvSE.T1","PAvSP.T1","SEvSP.T1")

## Make a master file of all results by merging each table by locus
glm_PAvSvSEvE_all <- merge(glm_PAvSvSEvE_all, glm_no3, by=c("CHROM","POS"))
glm_PAvSvSEvE_all <- merge(glm_PAvSvSEvE_all, glm_no7, by=c("CHROM","POS"))
glm_PAvSvSEvE_all <- merge(glm_PAvSvSEvE_all, glm_no11, by=c("CHROM","POS"))
glm_PAvSvSEvE_all <- merge(glm_PAvSvSEvE_all, glm_no15, by=c("CHROM","POS"))
glm_PAvSvSEvE_all <- merge(glm_PAvSvSEvE_all, glm_no21, by=c("CHROM","POS"))
glm_PAvSvSEvE_all <- merge(glm_PAvSvSEvE_all, glm_no27, by=c("CHROM","POS"))
glm_PAvSvSEvE_all <- merge(glm_PAvSvSEvE_all, glm_no33, by=c("CHROM","POS"))
glm_PAvSvSEvE_all <- merge(glm_PAvSvSEvE_all, glm_no37, by=c("CHROM","POS"))
glm_PAvSvSEvE_all <- merge(glm_PAvSvSEvE_all, glm_no41, by=c("CHROM","POS"))
glm_PAvSvSEvE_all <- merge(glm_PAvSvSEvE_all, glm_no45, by=c("CHROM","POS"))

# filter out chromosome 4
glm_PAvSvSEvE_all <- glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$CHROM != "4",]

# p-value correction
for(i in c(3:68)) {
  glm_PAvSvSEvE_all <- cbind(glm_PAvSvSEvE_all, p.adjust(glm_PAvSvSEvE_all[,i], method = "fdr"))
  colnames(glm_PAvSvSEvE_all)[i+66] <- paste(names(glm_PAvSvSEvE_all)[i], ".fdr", sep="")
}

# calculate -log10 for p-values
for(i in 3:68) {
  glm_PAvSvSEvE_all <- cbind(glm_PAvSvSEvE_all, as.numeric(-log10(glm_PAvSvSEvE_all[,i])))
  colnames(glm_PAvSvSEvE_all)[i+132] <- paste(names(glm_PAvSvSEvE_all)[i], ".logp", sep="")
}

#################################################################
### Check correlations of full glm results vs. leave-one-out  ###
### results for the "SEvSP" contrast.						  ###
#################################################################

## All samples vs. no cage 3 (SP)
cor.test(as.numeric(glm_PAvSvSEvE_all$SEvSP.T1),
		 as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no3), method="pearson")
#	Pearson's product-moment correlation
#
#data:  as.numeric(glm_PAvSvSEvE_all$SEvSP.T1) and 
#as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no3)
#t = 1843.7, df = 1587113, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.8251579 0.8261483
#sample estimates:
#      cor 
#0.8256537 

## All samples vs. no cage 7 (SP)
cor.test(as.numeric(glm_PAvSvSEvE_all$SEvSP.T1), as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no7), method="pearson")
#	Pearson's product-moment correlation
#
#data:  as.numeric(glm_PAvSvSEvE_all$SEvSP.T1) and #as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no7)
#t = 2338.1, df = 1587113, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.8799931 0.8806932
#sample estimates:
#      cor 
#0.8803437 

## All samples vs. no cage 11 (SE)
cor.test(as.numeric(glm_PAvSvSEvE_all$SEvSP.T1), as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no11), method="pearson")
#	Pearson's product-moment correlation
#
#data:  as.numeric(glm_PAvSvSEvE_all$SEvSP.T1) and #as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no11)
#t = 4458.6, df = 1587113, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.9622079 0.9624379
#sample estimates:
#      cor 
#0.9623231 

## All samples vs. no cage 15 (SE)
cor.test(as.numeric(glm_PAvSvSEvE_all$SEvSP.T1), as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no15), method="pearson")
#	Pearson's product-moment correlation
#
#data:  as.numeric(glm_PAvSvSEvE_all$SEvSP.T1) and #as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no15)
#t = 2418.9, df = 1587113, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.8865894 0.8872533
#sample estimates:
#      cor 
#0.8869218 

## All samples vs. no cage 21 (SE)
cor.test(as.numeric(glm_PAvSvSEvE_all$SEvSP.T1), as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no21), method="pearson")
#	Pearson's product-moment correlation
#
#data:  as.numeric(glm_PAvSvSEvE_all$SEvSP.T1) and #as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no21)
#t = 2871, df = 1587113, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.9154661 0.9159685
#sample estimates:
#      cor 
#0.9157176 

## All samples vs. no cage 27 (SE)
cor.test(as.numeric(glm_PAvSvSEvE_all$SEvSP.T1), as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no27), method="pearson")
#	Pearson's product-moment correlation
#
#data:  as.numeric(glm_PAvSvSEvE_all$SEvSP.T1) and #as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no27)
#t = 3654.8, df = 1587113, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.9452438 0.9455743
#sample estimates:
#      cor 
#0.9454093 

## All samples vs. no cage 33 (SP)
cor.test(as.numeric(glm_PAvSvSEvE_all$SEvSP.T1), as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no33), method="pearson")
#	Pearson's product-moment correlation
#
#data:  as.numeric(glm_PAvSvSEvE_all$SEvSP.T1) and #as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no33)
#t = 2244.2, df = 1587113, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.8716256 0.8723712
#sample estimates:
#      cor 
#30.8719989 

## All samples vs. no cage 37 (SP)
cor.test(as.numeric(glm_PAvSvSEvE_all$SEvSP.T1), as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no37), method="pearson")
#	Pearson's product-moment correlation
#
#data:  as.numeric(glm_PAvSvSEvE_all$SEvSP.T1) and #as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no37)
#t = 2394.9, df = 1587113, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.8846802 0.8853546
#sample estimates:
#      cor 
#0.8850178 

## All samples vs. no cage 41 (SE)
cor.test(as.numeric(glm_PAvSvSEvE_all$SEvSP.T1), as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no41), method="pearson")
#	Pearson's product-moment correlation
#
#data:  as.numeric(glm_PAvSvSEvE_all$SEvSP.T1) and #as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no41)
#t = 2685.1, df = 1587113, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.9050262 0.9055875
#sample estimates:
#      cor 
#0.9053072 

## All samples vs. no cage 45 (SE)
cor.test(as.numeric(glm_PAvSvSEvE_all$SEvSP.T1), as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no45), method="pearson")
#	Pearson's product-moment correlation
#
#data:  as.numeric(glm_PAvSvSEvE_all$SEvSP.T1) and #as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no45)
#t = 3588.9, df = 1587113, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.9433834 0.9437248
#sample estimates:
#      cor 
#0.9435543 



### Build Manhattan Plots ###

# No dropped cages for reference
manh_SEvSP_all <- ggplot(glm_PAvSvSEvE_all, aes(POS, SEvSP.T1.logp)) + 
	geom_line(alpha = 0.5, colour = "grey") + 
	 # highlight significant loci (FDR < 0.05) from GLM contrast using all samples
	geom_point(data=na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1.fdr < 0.05,]), aes(POS, SEvSP.T1.logp, color = "red"), size = 1) +
	# draw a dotted line at the minimum value of all significant loci (FDR < 0.05) 
	geom_hline(yintercept = min(na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1.fdr < 0.05,]$SEvSP.T1.logp)), color="black",linetype="dashed",linewidth=.25) +
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$SEvSP.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("SE vs SP: all samples")

# Dropped cage 3 (SP)
manh_SEvSP_no3 <- ggplot(glm_PAvSvSEvE_all, aes(POS, SEvSP.T1_no3.logp)) +
	geom_line(alpha = 0.5, colour = "grey") +  
	# highlight significant loci (FDR < 0.05) from GLM contrast using all samples
	geom_point(data=na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1.fdr < 0.05,]), aes(POS, SEvSP.T1_no3.logp, color = "red"), size = 1) +
	# draw a dotted line at the minimum value of all significant loci (FDR < 0.05) 
	geom_hline(yintercept = min(na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1_no3.fdr < 0.05,]$SEvSP.T1_no3.logp)), color="black",linetype="dashed",linewidth=.25) +
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$SEvSP.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("SE vs SP: no 3 (SP)")

# Dropped cage 7 (SP)
manh_SEvSP_no7 <- ggplot(glm_PAvSvSEvE_all, aes(POS, SEvSP.T1_no7.logp)) +
	geom_line(alpha = 0.5, colour = "grey") +  
	# highlight significant loci (FDR < 0.05) from GLM contrast using all samples
	geom_point(data=na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1.fdr < 0.05,]), aes(POS, SEvSP.T1_no7.logp, color = "red"), size = 1) +
	# draw a dotted line at the minimum value of all significant loci (FDR < 0.05) 
	geom_hline(yintercept = min(na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1_no7.fdr < 0.05,]$SEvSP.T1_no7.logp)) ,color="black",linetype="dashed",linewidth=.25) +
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$SEvSP.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("SE vs SP: no 7 (SP)")

# Dropped cage 11 (SE)
manh_SEvSP_no11 <- ggplot(glm_PAvSvSEvE_all, aes(POS, SEvSP.T1_no11.logp)) +
	geom_line(alpha = 0.5, colour = "grey") +  
	# highlight significant loci (FDR < 0.05) from GLM contrast using all samples
	geom_point(data=na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1.fdr < 0.05,]), aes(POS, SEvSP.T1_no11.logp, color = "red"), size = 1) +
	# draw a dotted line at the minimum value of all significant loci (FDR < 0.05) 
	geom_hline(yintercept = min(na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1_no11.fdr < 0.05,]$SEvSP.T1_no11.logp)), color="black",linetype="dashed",linewidth=.25) +
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$SEvSP.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("SE vs SP: no 11 (SE)")

# Dropped cage 15 (SE)
manh_SEvSP_no15 <- ggplot(glm_PAvSvSEvE_all, aes(POS, SEvSP.T1_no15.logp)) +
	geom_line(alpha = 0.5, colour = "grey") +  
	# highlight significant loci (FDR < 0.05) from GLM contrast using all samples
	geom_point(data=na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1.fdr < 0.05,]), aes(POS, SEvSP.T1_no15.logp, color = "red"), size = 1) +
	# draw a dotted line at the minimum value of all significant loci (FDR < 0.05) 
	geom_hline(yintercept = min(na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1_no15.fdr < 0.05,]$SEvSP.T1_no15.logp)), color="black",linetype="dashed",linewidth=.25) +
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$SEvSP.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("SE vs SP: no 15 (SE)")
  	
# Dropped cage 21 (SE)
manh_SEvSP_no21 <- ggplot(glm_PAvSvSEvE_all, aes(POS, SEvSP.T1_no21.logp)) +
	geom_line(alpha = 0.5, colour = "grey") + 
	# highlight significant loci (FDR < 0.05) from GLM contrast using all samples
	geom_point(data=na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1.fdr < 0.05,]), aes(POS, SEvSP.T1_no21.logp, color = "red"), size = 1) +
	# draw a dotted line at the minimum value of all significant loci (FDR < 0.05) 
	geom_hline(yintercept = min(na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1_no21.fdr < 0.05,]$SEvSP.T1_no21.logp)), color="black",linetype="dashed",linewidth=.25) +
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$SEvSP.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("SE vs SP: no 21 (SE)")

# Dropped cage 27 (SE)
manh_SEvSP_no27 <- ggplot(glm_PAvSvSEvE_all, aes(POS, SEvSP.T1_no27.logp)) +
	geom_line(alpha = 0.5, colour = "grey") +  
	# highlight significant loci (FDR < 0.05) from GLM contrast using all samples
	geom_point(data=na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1.fdr < 0.05,]), aes(POS, SEvSP.T1_no27.logp, color = "red"), size = 1) +
	# draw a dotted line at the minimum value of all significant loci (FDR < 0.05) 
	geom_hline(yintercept = min(na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1_no27.fdr < 0.05,]$SEvSP.T1_no27.logp)), color="black",linetype="dashed",linewidth=.25) +
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$SEvSP.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("SE vs SP: no 27 (SE)")

# Dropped cage 33 (SP)
manh_SEvSP_no33 <- ggplot(glm_PAvSvSEvE_all, aes(POS, SEvSP.T1_no33.logp)) +
	geom_line(alpha = 0.5, colour = "grey") +  
	# highlight significant loci (FDR < 0.05) from GLM contrast using all samples
	geom_point(data=na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1.fdr < 0.05,]), aes(POS, SEvSP.T1_no33.logp, color = "red"), size = 1) +
	# draw a dotted line at the minimum value of all significant loci (FDR < 0.05) 
	geom_hline(yintercept = min(na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1_no33.fdr < 0.05,]$SEvSP.T1_no33.logp)), color="black",linetype="dashed",linewidth=.25) +
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$SEvSP.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("SE vs SP: no 33 (SP)")
  
# Dropped cage 37 (SP)
manh_SEvSP_no37 <- ggplot(glm_PAvSvSEvE_all, aes(POS, SEvSP.T1_no37.logp)) + 
	geom_line(alpha = 0.5, colour = "grey") +  
	# highlight significant loci (FDR < 0.05) from GLM contrast using all samples
	geom_point(data=na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1.fdr < 0.05,]), aes(POS, SEvSP.T1_no37.logp, color = "red"), size = 1) +
	# draw a dotted line at the minimum value of all significant loci (FDR < 0.05) 
	geom_hline(yintercept = min(na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1_no37.fdr < 0.05,]$SEvSP.T1_no37.logp)), color="black",linetype="dashed",linewidth=.25) +
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$SEvSP.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("SE vs SP: no 37 (SP)")
  	
# Dropped cage 41 (SE)
manh_SEvSP_no41 <- ggplot(glm_PAvSvSEvE_all, aes(POS, SEvSP.T1_no41.logp)) + 
	geom_line(alpha = 0.5, colour = "grey") + 
	# highlight significant loci (FDR < 0.05) from GLM contrast using all samples
	geom_point(data=na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1.fdr < 0.05,]), aes(POS, SEvSP.T1_no41.logp, color = "red"), size = 1) +
	# draw a dotted line at the minimum value of all significant loci (FDR < 0.05) 
	geom_hline(yintercept = min(na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1_no41.fdr < 0.05,]$SEvSP.T1_no41.logp)), color="black",linetype="dashed",linewidth=.25) +
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$SEvSP.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("SE vs SP: no 41 (SE)")

# Dropped cage 45 (SE)
manh_SEvSP_no45 <- ggplot(glm_PAvSvSEvE_all, aes(POS, SEvSP.T1_no45.logp)) +
	geom_line(alpha = 0.5, colour = "grey") +  
	# highlight significant loci (FDR < 0.05) from GLM contrast using all samples
	geom_point(data=na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1.fdr < 0.05,]), aes(POS, SEvSP.T1_no45.logp, color = "red"), size = 1) +
	# draw a dotted line at the minimum value of all significant loci (FDR < 0.05) 
	geom_hline(yintercept = min(na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1_no45.fdr < 0.05,]$SEvSP.T1_no45.logp)), color="black",linetype="dashed",linewidth=.25) +
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$SEvSP.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("SE vs SP: no 45 (SE)")

 
### Figure S3 - Plot all leave-one-out manhattans together 
pdf(file = "rudflies_2023_redo.SEvSP.T1.LOO.manh.pdf", width=10, height=12)
	ggarrange(manh_SEvSP_all, manh_SEvSP_no3, manh_SEvSP_no7, manh_SEvSP_no33, manh_SEvSP_no37, manh_SEvSP_no11, manh_SEvSP_no15, manh_SEvSP_no21, manh_SEvSP_no27, manh_SEvSP_no41, manh_SEvSP_no45,
              ncol = 2, nrow = 6)
dev.off()


############################################################
### Build bonus Manhattan plots for the "EvSP" contrast. ###
############################################################

# No dropped cages for reference
manh_EvSP_all <- ggplot(glm_PAvSvSEvE_all, aes(POS, EvSP.T1.logp)) + 
	geom_line(alpha = 0.5, colour = "grey") +  
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$EvSP.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("E vs SP: all samples")

# Dropped cage 3 (SP)
manh_EvSP_no3 <- ggplot(glm_PAvSvSEvE_all, aes(POS, EvSP.T1_no3.logp)) + 
	geom_line(alpha = 0.5, colour = "grey") +  	
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$EvSP.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("E vs SP: no 3")

# Dropped cage 7 (SP)
manh_EvSP_no7 <- ggplot(glm_PAvSvSEvE_all, aes(POS, EvSP.T1_no7.logp)) + 
	geom_line(alpha = 0.5, colour = "grey") +  
	# highlight significant loci (FDR < 0.05) from GLM contrast using all samples
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$EvSP.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("E vs SP: no 7")

# Dropped cage 33 (SP)
manh_EvSP_no33 <- ggplot(glm_PAvSvSEvE_all, aes(POS, EvSP.T1_no33.logp)) + 
	geom_line(alpha = 0.5, colour = "grey") +  
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$EvSP.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("E vs SP: no 33")

# Dropped cage 37 (SP)
manh_EvSP_no37 <- ggplot(glm_PAvSvSEvE_all, aes(POS, EvSP.T1_no37.logp)) + 
	geom_line(alpha = 0.5, colour = "grey") +  
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$EvSP.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("E vs SP: no 37")
  
### Unused in paper - Plot all leave-one-out manhattans together 
pdf(file = "rudflies_2023_redo.EvSP.T1.LOO.manh.pdf", width=10, height=14)
	ggarrange(manh_EvSP_all, manh_EvSP_no3, manh_EvSP_no7, manh_EvSP_no33, manh_EvSP_no37,
              ncol = 1, nrow = 5)
dev.off()


############################################################
### Build bonus Manhattan plots for the "EvSE" contrast. ###
############################################################

# No dropped cages for reference
manh_EvSE_all <- ggplot(glm_PAvSvSEvE_all, aes(POS, EvSE.T1.logp)) + 
	geom_line(alpha = 0.5, colour = "grey") +  
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$EvSE.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("E vs SE: all samples")

# Dropped cage 11 (SE)
manh_EvSE_no11 <- ggplot(glm_PAvSvSEvE_all, aes(POS, EvSE.T1_no11.logp)) + 
	geom_line(alpha = 0.5, colour = "grey") +  
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$EvSE.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("E vs SE: no 11")

# Dropped cage 15 (SE)
manh_EvSE_no15 <- ggplot(glm_PAvSvSEvE_all, aes(POS, EvSE.T1_no15.logp)) + 
	geom_line(alpha = 0.5, colour = "grey") + 
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$EvSE.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("E vs SE: no 15")

# Dropped cage 21 (SE)
manh_EvSE_no21 <- ggplot(glm_PAvSvSEvE_all, aes(POS, EvSE.T1_no21.logp)) + 
	geom_line(alpha = 0.5, colour = "grey") +  
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$EvSE.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("E vs SE: no 21")
  
# Dropped cage 27 (SE)
manh_EvSE_no27 <- ggplot(glm_PAvSvSEvE_all, aes(POS, EvSE.T1_no27.logp)) + 
	geom_line(alpha = 0.5, colour = "grey") +  
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$EvSE.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("E vs SE: no 27")
  	
# Dropped cage 41 (SE)
manh_EvSE_no41 <- ggplot(glm_PAvSvSEvE_all, aes(POS, EvSE.T1_no41.logp)) + 
	geom_line(alpha = 0.5, colour = "grey") +  
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$EvSE.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("E vs SE: no 41")
  	
# Dropped cage 45 (SE)
manh_EvSE_no45 <- ggplot(glm_PAvSvSEvE_all, aes(POS, EvSE.T1_no45.logp)) + 
	geom_line(alpha = 0.5, colour = "grey") +  
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$EvSE.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("E vs SE: no 45")

### Unused in paper - Plot all leave-one-out manhattans together 
pdf(file = "rudflies_2023_redo.EvSE.T1.LOO.manh.pdf", width=10, height=18)
	ggarrange(manh_EvSE_all, manh_EvSE_no11, manh_EvSE_no15, manh_EvSE_no21, manh_EvSE_no27, manh_EvSE_no41, manh_EvSE_no45,
              ncol = 2, nrow = 7)
dev.off()



#################
### Figure S8 ###
#################

## Goal here is to calculate pairwise distances of TPT4 SP vs E outlier SNPs among T1 E 
## populations, then and compare with distances between E and other TPT1 treatments.
## Mainly, we want to see if T4 adaptive convergent candidates are within the range 
## negative control population variation to confirm they did not start divergent through 
## drift or founder effect.

## First, let's grab the allele frequencies for SNPs in the filtered list of GLM results
haf.freq.glm.filt <- merge(glm.all.rolwin20[,c(1:2)], haf.freq, by=c("CHROM","POS"))

## Select TPT4 SP vs E outliers (FDR < 0.05) from frequency table
haf.freq.glm.SvE.T4.filt <- na.omit(haf.freq.glm.filt[glm.all.rolwin20$SvE.T4.fdr.rolwin20 < 0.05,])

## Alternatively, just pick the previously-discovered TPT4 convergent candidates
haf.freq.glm.SvE.T4.filt <- merge(unique(glm.all.rolwin20.TPT4.convergent_candidates[,c(1:2)]), haf.freq, by=c("CHROM","POS"))

## How many are there?
dim(haf.freq.glm.SvE.T4.filt)
#[1] 1988  191

## Now we want to grab some key metadata columns for TPT1 treatment comparisons 
## Specifically, we'll grab sample ID and condition
haf.meta.T1filt.slim <- haf.meta.T1filt[,c(3,13)]

## Pull only TPT1 samples from the TPT4 convergent candidate AF table
haf.freq.glm.SvE.T4.cand.T1filt <- unique(subset(haf.freq.glm.SvE.T4.filt, select = as.character(haf.meta.T1filt$samp)))

## Make a Euclidean distance matrix from the filtered frequency table
haf.freq.glm.SvE.T4.cand.T1.dist <- as.matrix(dist(t(haf.freq.glm.SvE.T4.cand.T1filt), method="euclidean", diag=TRUE, upper=FALSE))

## Convert 0's to NA (since they are self-self comparisons and not useful)
haf.freq.glm.SvE.T4.cand.T1.dist[haf.freq.glm.SvE.T4.cand.T1.dist == 0] <- NA

## Melt the distance matrix into a dataframe that we can use for plotting
haf.freq.glm.SvE.T4.cand.T1.dist.df <- na.omit(reshape2::melt(as.matrix(haf.freq.glm.SvE.T4.cand.T1.dist), varnames = c("row", "col")))

## Now we want to grab some key metadata columns for TPT1 treatment comparisons 
## Specifically, we'll grab sample ID and condition
haf.meta.T1filt.slim <- haf.meta.T1filt[,c(3,13)]

## Merge them with the above distance dataframe by sample A
haf.freq.glm.SvE.T4.cand.T1.dist.df <- merge(haf.freq.glm.SvE.T4.cand.T1.dist.df, haf.meta.T1filt.slim, by.x="row", by.y="samp")

## Merge them with the above distance dataframe by sample B
haf.freq.glm.SvE.T4.cand.T1.dist.df <- merge(haf.freq.glm.SvE.T4.cand.T1.dist.df, haf.meta.T1filt.slim, by.x="col", by.y="samp")

## Now, using the metadata for samples A and B, we'll label treatment contrasts 
## First, we'll make a black treatment contrast
haf.freq.glm.SvE.T4.cand.T1.dist.df$contrast <- "NA"

## Label SEvSE
haf.freq.glm.SvE.T4.cand.T1.dist.df[haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.x == "SE" & haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.y == "SE",]$contrast <- "SEvSE"

## Label SPvSP
haf.freq.glm.SvE.T4.cand.T1.dist.df[haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.x == "SP" & haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.y == "SP",]$contrast <- "SPvSP"

## Label EvE
haf.freq.glm.SvE.T4.cand.T1.dist.df[haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.x == "E" & haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.y == "E",]$contrast <- "EvE"

## Label PAvPA
haf.freq.glm.SvE.T4.cand.T1.dist.df[haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.x == "PA" & haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.y == "PA",]$contrast <- "PAvPA"

## Label SEvSP
haf.freq.glm.SvE.T4.cand.T1.dist.df[haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.x == "SP" & haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.y == "SE",]$contrast <- "SEvSP"

haf.freq.glm.SvE.T4.cand.T1.dist.df[haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.x == "SE" & haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.y == "SP",]$contrast <- "SEvSP"

## Label EvSE
haf.freq.glm.SvE.T4.cand.T1.dist.df[haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.x == "E" & haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.y == "SE",]$contrast <- "EvSE"

haf.freq.glm.SvE.T4.cand.T1.dist.df[haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.x == "SE" & haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.y == "E",]$contrast <- "EvSE"

## Label EvSP
haf.freq.glm.SvE.T4.cand.T1.dist.df[haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.x == "E" & haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.y == "SP",]$contrast <- "EvSP"

haf.freq.glm.SvE.T4.cand.T1.dist.df[haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.x == "SP" & haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.y == "E",]$contrast <- "EvSP"

## Label EvPA
haf.freq.glm.SvE.T4.cand.T1.dist.df[haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.x == "E" & haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.y == "PA",]$contrast <- "EvPA"

haf.freq.glm.SvE.T4.cand.T1.dist.df[haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.x == "PA" & haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.y == "E",]$contrast <- "EvPA"

## Label PAvSE
haf.freq.glm.SvE.T4.cand.T1.dist.df[haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.x == "SE" & haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.y == "PA",]$contrast <- "PAvSE"

haf.freq.glm.SvE.T4.cand.T1.dist.df[haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.x == "PA" & haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.y == "SE",]$contrast <- "PAvSE"

## Label PAvSP
haf.freq.glm.SvE.T4.cand.T1.dist.df[haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.x == "SP" & haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.y == "PA",]$contrast <- "PAvSP"

haf.freq.glm.SvE.T4.cand.T1.dist.df[haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.x == "PA" & haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.y == "SP",]$contrast <- "PAvSP"

## Format from plotting
haf.freq.glm.SvE.T4.dist.df <- as.data.frame(unique(haf.freq.glm.SvE.T4.cand.T1.dist.df[,c(6,3)]))
haf.freq.glm.SvE.T4.dist.df$value <- as.numeric(haf.freq.glm.SvE.T4.dist.df$value)


### Plot all contrasts
## Boxplot
pdf(file = "haf.freq.glm.SvE.T4.T1dist.df.boxplot.pdf")
	ggplot(data = haf.freq.glm.SvE.T4.dist.df, aes(x=contrast, y=value)) + 
	geom_boxplot(aes(colour=contrast)) + 
	ggtitle("Pairwise Euclidean distances: T4 candidate loci (AF) among T1 populations") + 
	theme_classic()
dev.off()

## Jitter plot
pdf(file = "haf.freq.glm.SvE.T4.T1dist.df.jitter.pdf")
	ggplot(data = haf.freq.glm.SvE.T4.dist.df, aes(x=contrast, y=value)) + 
	geom_jitter(aes(colour=contrast)) + 
	ggtitle("Pairwise Euclidean distances: T4 candidate loci (AF) among T1 populations") + 
	theme_classic()
dev.off()

### Plot without PA
## Take out PA-related contrasts
haf.freq.glm.SvE.T4.dist.df.noPA <- haf.freq.glm.SvE.T4.dist.df[haf.freq.glm.SvE.T4.dist.df$contrast != "EvPA" & haf.freq.glm.SvE.T4.dist.df$contrast != "PAvPA" & haf.freq.glm.SvE.T4.dist.df$contrast != "PAvSE" & haf.freq.glm.SvE.T4.dist.df$contrast != "PAvSP",]

## Boxplot
pdf(file = "haf.freq.glm.SvE.T4.T1dist.df.noPA.boxplot.pdf")
	ggplot(data = haf.freq.glm.SvE.T4.dist.df.noPA, aes(x=contrast, y=value)) + 
	geom_boxplot(aes(colour=contrast)) + 
	ggtitle("Pairwise Euclidean distances: T4 candidate loci (AF) among T1 populations") + 
	theme_classic()
dev.off()

## Jitter plot
pdf(file = "haf.freq.glm.SvE.T4.T1dist.df.noPA.jitter.pdf")
	ggplot(data = haf.freq.glm.SvE.T4.dist.df.noPA, aes(x=contrast, y=value)) + 
	geom_jitter(aes(colour=contrast)) + 
	ggtitle("Pairwise Euclidean distances: T4 candidate loci (AF) among T1 populations") + 
	theme_classic()
dev.off()



##################
### Figure S11 ###
##################

### Checking whether the hafpipe frequencies match the dgrp2 founder file frequencies ###
## Load mean frequencies of the relevant DGRP2 samples from the DGRP Freeze.2 VCF file
dgrp2.freq <- read.delim("../../dgrp2_founder_list.frq", header=TRUE, sep = "\t")
## Load mean frequencies of experimental samples from VCF generated by BCFtools
vcf.freq <- read.delim("../calling/filtered-all.frq", header=TRUE, sep = "\t")
## Load mean frequencies of experimental sample from frequency table generated by Hafpipe
haf.freq.means <- cbind(haf.freq[,c(1:2)],rowMeans(haf.freq[,c(3:ncol(haf.freq))]))
## Turns out Hafpipe frequencies are generated for the ALT alleles frequencies instead of 
## REF alleles, so we'll calculate REF frequencies (1 - ALT) to positively correlate with 
## other datasets.
colnames(haf.freq.means) <- c("CHROM","POS","hafpipe.means.A")
haf.freq.means$hafpipe.means.R <- 1-haf.freq.means$hafpipe.means.A
## Merge all three tables by locus
all.freq <- merge(dgrp2.freq, vcf.freq, by=c("CHROM","POS"))
all.freq <- merge(all.freq, haf.freq.means, by=c("CHROM","POS"))

## Create a density plot comparing DGRP2 founder frequencies and BCFtools VCF frequencies
pdf(file = "rudflies_2023_redo.dgrp2_vs_vcf_freqs.density.pdf", width=6.5, height=5)
	ggplot(all.freq, 
	  aes(x=R.FREQ.x, y=R.FREQ.y)) + 
	  geom_bin2d(bins = 50) +
	  scale_fill_continuous(type = "viridis") +
	  xlim(0.05, 0.95) +
	  ylim(0.05, 0.95) +
	  xlab("dgrp2 AF") + 
	  ylab("bcftools AF") +
	  theme_classic() +
	  theme_update(axis.ticks.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.y = element_blank(),
               axis.text.y = element_blank())
dev.off()

## Test correlation
cor.test(as.numeric(all.freq$R.FREQ.x), as.numeric(all.freq$R.FREQ.y), method="pearson")

## Create a density plot comparing DGRP2 founder frequencies and Hafpipe frequencies
pdf(file = "rudflies_2023_redo.dgrp2_vs_hafpipe_freqs.density.pdf", width=6.5, height=5)
	ggplot(all.freq, 
	  aes(x=R.FREQ.x, y=hafpipe.means.R)) + 
	  geom_bin2d(bins = 50) +
	  scale_fill_continuous(type = "viridis") +
	  xlim(0.05, 0.95) +
	  ylim(0.05, 0.95) +
	  xlab("dgrp2 AF") + 
	  ylab("hafpipe AF") +
	  theme_classic() +
	  theme_update(axis.ticks.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.y = element_blank(),
               axis.text.y = element_blank())
dev.off()

## Test correlation
cor.test(as.numeric(all.freq$R.FREQ.x), as.numeric(all.freq$hafpipe.means.R), method="pearson")

## Create a density plot comparing BCFtools VCF frequencies and Hafpipe frequencies
pdf(file = "rudflies_2023_redo.vcf_vs_hafpipe_freqs.density.pdf", width=6.5, height=5)
	ggplot(all.freq, 
	  aes(x=R.FREQ.y, y=hafpipe.means)) + 
	  geom_bin2d(bins = 50) +
	  scale_fill_continuous(type = "viridis") +
	  xlim(0.05, 0.95) +
	  ylim(0.05, 0.95) +
	  xlab("bcftools AF") + 
	  ylab("hafpipe AF") +
	  theme_classic() +
	  theme_update(axis.ticks.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.y = element_blank(),
               axis.text.y = element_blank())
dev.off()

## Test correlation 
cor.test(as.numeric(all.freq$R.FREQ.y), as.numeric(all.freq$hafpipe.means.R), method="pearson")

## Overall takeaway, hafpipe output is in reference to the alt-allele, not the ref
## Plus, hafpipe frequencies are much better correlated with the founder file than the bcftools output.


##################
### Figure S12 ###
##################

## Use this color scheme for "condition" side colors
sidecols <- haf.meta.T1filt$condition 
sidecols <- gsub("PA","#495184",sidecols)
sidecols <- gsub("SP","#D9B851",sidecols)
sidecols <- gsub("SE","#848556",sidecols)
sidecols <- gsub("E","#D26183",sidecols)

### Heatmaps for TPT1 SE outliers SNPs among TPT1 AF frequencies ###

## Select SNPs from AF table significant (FDR<0.05) in GLM contrasts between E vs SE:TPT1
haf.freq.SE.fdr05 <- as.matrix(merge(na.omit(glm.all.rolwin20[glm.all.rolwin20$EvSE.T1.fdr.rolwin20 < 0.05,c(1:2)]),cbind(haf.sites.T1filt,haf.freq.T1filt), by=c("CHROM","POS"))[,-c(1:2)])

## Make the heatmap
pdf(file = "rudflies_2023_redo.haf.freq.SE.fdr05.TPT1.heatmap.pdf", width=10, height=10)
	heatmap.2(unique(haf.freq.SE.fdr05), Rowv=FALSE, dendrogram="col", ColSideColors=sidecols, labRow=FALSE, trace = "none")
dev.off()

## Select SNPs from AF table significant (FDR<0.01) in GLM contrasts between E vs SE:TPT1
haf.freq.SE.fdr01 <- as.matrix(merge(na.omit(glm.all.rolwin20[glm.all.rolwin20$EvSE.T1.fdr.rolwin20 < 0.01,c(1:2)]),cbind(haf.sites.T1filt,haf.freq.T1filt), by=c("CHROM","POS"))[,-c(1:2)])

## Make the heatmap
pdf(file = "rudflies_2023_redo.haf.freq.SE.fdr01.TPT1.heatmap.pdf", width=10, height=10)
	heatmap.2(haf.freq.SE.fdr01, Rowv=FALSE, dendrogram="col", ColSideColors=sidecols, labRow=FALSE, trace = "none")
dev.off()


### Heatmaps for TPT4 SP outliers SNPs among TPT1 AF frequencies ###

## Select SNPs from AF table significant (FDR<0.05) in GLM contrasts between E vs SP:TPT4
haf.freq.SP.fdr05 <- as.matrix(merge(na.omit(glm.all.rolwin20[glm.all.rolwin20$SvE.T4.fdr.rolwin20 < 0.05,c(1:2)]),cbind(haf.sites.T1filt,haf.freq.T1filt), by=c("CHROM","POS"))[,-c(1:2)])

## Make the heatmap
pdf(file = "rudflies_2023_redo.haf.freq.SP.fdr05.TPT1.heatmap.pdf", width=10, height=10)
	heatmap.2(haf.freq.SP.fdr05, Rowv=FALSE, dendrogram="col", ColSideColors=sidecols, labRow=FALSE, trace = "none")
dev.off()

## Select SNPs from AF table significant (FDR<0.01) in GLM contrasts between E vs SP:TPT4
haf.freq.SP.fdr01 <- as.matrix(merge(na.omit(glm.all.rolwin20[glm.all.rolwin20$SvE.T4.fdr.rolwin20 < 0.01,c(1:2)]),cbind(haf.sites.T1filt,haf.freq.T1filt), by=c("CHROM","POS"))[,-c(1:2)])

## Make the heatmap
pdf(file = "rudflies_2023_redo.haf.freq.SP.fdr01.TPT1.heatmap.pdf", width=10, height=10)
	heatmap.2(haf.freq.SP.fdr01, Rowv=FALSE, dendrogram="col", ColSideColors=sidecols, labRow=FALSE, trace = "none")
dev.off()


### Now we're switching to outlier heatmaps among TPT4-sample frequencies

## Make TPT4 AF table subset
## Filter metadata for TPT4 samples for PA, S, and E
haf.meta.T4filt <- haf.meta[haf.meta$tpt == "4" & haf.meta$batch == "a" & (haf.meta$treat.fix == "PA" | haf.meta$treat.fix == "S" | haf.meta$treat.fix == "E"),]

## Filter AF table for TPT4 samples for PA, S, and E
haf.freq.T4filt <- haf.freq[, which((names(haf.freq) %in% haf.meta.T4filt$samp)==TRUE)]
haf.meta.T4filt <- haf.meta.T4filt[which((haf.meta.T4filt$samp %in% names(haf.freq.T4filt))==TRUE),]

#additional filter to remove extreme low variance loci
haf.sites.T4filt <- haf.freq[rowVars(as.matrix(haf.freq.T4filt))>0.001,c(1:2)]
haf.freq.T4filt <- haf.freq.T4filt[rowVars(as.matrix(haf.freq.T4filt))>0.001,]

#tables need a bit of reformatting for plotting
haf.meta.T4filt$treat.fix <- as.factor(haf.meta.T4filt$treat.fix)
haf.meta.T4filt$tpt <- as.factor(haf.meta.T4filt$tpt)

## Use this color scheme for "treatment" side colors
sidecols2 <- haf.meta.T4filt$treat.fix 
sidecols2 <- gsub("PA","#495184",sidecols2)
sidecols2 <- gsub("S","#D9B851",sidecols2)
sidecols2 <- gsub("E","#D26183",sidecols2)

### Heatmaps for TPT4 SP outliers SNPs among TPT1 AF frequencies ###

## Select SNPs from AF table significant (FDR<0.05) in GLM contrasts between E vs SE:TPT1
haf.freq.SE.t4.fdr05 <- as.matrix(merge(na.omit(glm.all.rolwin20[glm.all.rolwin20$EvSE.T1.fdr.rolwin20 < 0.05,c(1:2)]),cbind(haf.sites.T4filt,haf.freq.T4filt), by=c("CHROM","POS"))[,-c(1:2)])

## Make the heatmap
pdf(file = "rudflies_2023_redo.haf.freq.SE.fdr05.TPT4.heatmap.pdf", width=10, height=10)
	heatmap.2(haf.freq.SE.t4.fdr05, Rowv=FALSE, dendrogram="col", ColSideColors=sidecols2, labRow=FALSE, trace = "none")
dev.off()

## Select SNPs from AF table significant (FDR<0.01) in GLM contrasts between E vs SE:TPT1
haf.freq.SE.t4.fdr01 <- as.matrix(merge(na.omit(glm.all.rolwin20[glm.all.rolwin20$EvSE.T1.fdr.rolwin20 < 0.01,c(1:2)]),cbind(haf.sites.T4filt,haf.freq.T4filt), by=c("CHROM","POS"))[,-c(1:2)])

## Make the heatmap
pdf(file = "rudflies_2023_redo.haf.freq.SE.fdr01.TPT4.heatmap.pdf", width=10, height=10)
	heatmap.2(haf.freq.SE.t4.fdr01, Rowv=FALSE, dendrogram="col", ColSideColors=sidecols2, labRow=FALSE, trace = "none")
dev.off()

## Select SNPs from AF table significant (FDR<0.05) in GLM contrasts between E vs SP:TPT4
haf.freq.SP.t4.fdr05 <- as.matrix(merge(na.omit(glm.all.rolwin20[glm.all.rolwin20$SvE.T4.fdr.rolwin20 < 0.05,c(1:2)]),cbind(haf.sites.T4filt,haf.freq.T4filt), by=c("CHROM","POS"))[,-c(1:2)])

## Make the heatmap
pdf(file = "rudflies_2023_redo.haf.freq.SP.fdr05.TPT4.heatmap.pdf", width=10, height=10)
	heatmap.2(haf.freq.SP.t4.fdr05, Rowv=FALSE, dendrogram="col", ColSideColors=sidecols2, labRow=FALSE, trace = "none")
dev.off()

## Select SNPs from AF table significant (FDR<0.01) in GLM contrasts between E vs SP:TPT4
haf.freq.SP.t4.fdr01 <- as.matrix(merge(na.omit(glm.all.rolwin20[glm.all.rolwin20$SvE.T4.fdr.rolwin20 < 0.01,c(1:2)]),cbind(haf.sites.T4filt,haf.freq.T4filt), by=c("CHROM","POS"))[,-c(1:2)])

## Make the heatmap
pdf(file = "rudflies_2023_redo.haf.freq.SP.fdr01.TPT4.heatmap.pdf", width=10, height=10)
	heatmap.2(haf.freq.SP.t4.fdr01, Rowv=FALSE, dendrogram="col", ColSideColors=sidecols2, labRow=FALSE, trace = "none")
dev.off()

