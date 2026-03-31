### In R ###
#configure r environment
setwd("/scratch/user/jfaberha/20260325_052109/admera/gp_analysis/rudflies_2023_redo/r")
library(emmeans)
library(matrixStats)
#install.packages("~/Downloads/ACER-master", repos=NULL, type="source")
library(ACER)
library(poolSeq)
library(ggplot2)
library(ggpubr)
library(stats)
library(ggfortify)
library(dplyr)
library(zoo)
library(rstatix)
#install.packages("slider")
library(slider)
library(RColorBrewer)


#################################
### Load required input files ###
#################################

## Variant Effect Predictor (VEP) annotation file with appended FLYCADD scores
vep <- read.table("filtered-all.annot.vcf.FLYCADD.tsv", header=TRUE) 
## More detailed gene information from gff
snp.gff.overlap <- read.delim("rudflies_2023_hafpipe_loci_genic_overlap_info.bed", sep="\t", header=TRUE)
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

## Merge annotations for filtering and highlighting
glm.all.rolwin20.annot <- merge(glm.all.rolwin20, vep, by=c("CHROM","POS"))

## Load spinosad-resistance candidate gene list
spino.cand <- read.table("spino.cand.list.txt", header=FALSE)
names(spino.cand) <- "Gene"

## Now merge to find All SNPs in and around candidate genes
glm.all.rolwin20.annot.spino <- merge(spino.cand, glm.all.rolwin20.annot, by="Gene")


####################################
### Hypergeometric overlap tests ###
####################################

## Calculate minimum p-values and FDR values per annotated spino candidate genes. We will 
## assign the most significant FDR value overlapping the gene, plus its immediate upstream 
## and downstream regions, to represent the gene in comparison with background gene sets.

## Pull FDR columns only for spino candidates
glm.all.fdr.rolwin20.annot.spino <- unique(cbind(glm.all.rolwin20.annot.spino[,c(1:3)],glm.all.rolwin20.annot.spino$Consequence,select(glm.all.rolwin20.annot.spino,contains("fdr"))))

## Loop through each GLM-results FDR column to calculate spino candidate minimum values
glm.all.fdr.rolwin20.annot.spino.min <- c() # create an object
#loop through all cols fdr values for candidate genes
for(g in 5:ncol(glm.all.fdr.rolwin20.annot.spino)) {
  tryCatch({
  	#create 2-column table with gene name and one FDR column at a time
    temp.min <- as.data.frame(cbind(glm.all.fdr.rolwin20.annot.spino$Gene,glm.all.fdr.rolwin20.annot.spino[,g]))
    #name columns
    colnames(temp.min)[1] <- "Gene" 
    colnames(temp.min)[2] <- "fdr"
    #assign minimum FDR values
    glm.all.fdr.rolwin20.annot.spino.min <- cbind(glm.all.fdr.rolwin20.annot.spino.min,
    	temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(fdr),
        	list(min = min)) %>%
        	pull(min))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.fdr.rolwin20.annot.spino.min) <- names(glm.all.fdr.rolwin20.annot.spino)[c(5:ncol(glm.all.fdr.rolwin20.annot.spino))]
## Add gene names to the table
glm.all.fdr.rolwin20.annot.spino.min <- cbind(Gene=temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(fdr),
        	list(min = min)) %>%
        	pull(Gene),glm.all.fdr.rolwin20.annot.spino.min)
        	
## Save minimum FDR results for downstream use
write.table(glm.all.fdr.rolwin20.annot.spino.min, file="rudflies_2023_redo.glm.all.rolwin20.minfdr.spino.txt", quote = FALSE, sep = "\t", row.names = F)
        	
## Pull FDR columns only for all genes/features
glm.all.fdr.rolwin20.annot <- unique(cbind(Gene=glm.all.rolwin20.annot$Gene,glm.all.rolwin20.annot[,c(1,2)],Consequence=glm.all.rolwin20.annot$Consequence,select(glm.all.rolwin20.annot,contains("fdr"))))

## Loop through each GLM-results FDR column to calculate minimum values for all genes
glm.all.fdr.rolwin20.annot.min <- c() # create an object
#loop through all cols fdr values for candidate genes
for(g in 5:ncol(glm.all.fdr.rolwin20.annot)) { #nrow(haf.freq.cand)
  tryCatch({
  	#create 2-column table with gene name and one FDR column at a time
    temp.min <- as.data.frame(cbind(glm.all.fdr.rolwin20.annot$Gene,glm.all.fdr.rolwin20.annot[,g]))
    #name columns
    colnames(temp.min)[1] <- "Gene" 
    colnames(temp.min)[2] <- "fdr"    
    #assign minimum FDR values
    glm.all.fdr.rolwin20.annot.min <- cbind(glm.all.fdr.rolwin20.annot.min,
    	temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(fdr),
        	list(min = min)) %>%
        	pull(min))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.fdr.rolwin20.annot.min) <- names(glm.all.fdr.rolwin20.annot)[c(5:ncol(glm.all.fdr.rolwin20.annot))]
## Add gene names to the table
glm.all.fdr.rolwin20.annot.min <- cbind(Gene=temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(fdr),
        	list(min = min)) %>%
        	pull(Gene),glm.all.fdr.rolwin20.annot.min)


### Redo above with missense SNPs in candidates only ###

## Merge to find missense SNPs in and around candidate genes
glm.all.rolwin20.annot.spino.missense <- merge(spino.cand, glm.all.rolwin20.annot[glm.all.rolwin20.annot$Consequence=="missense_variant",], by="Gene")

## Pull FDR columns only for spino candidate missense SNPs
glm.all.fdr.rolwin20.annot.spino.missense <- unique(cbind(glm.all.rolwin20.annot.spino.missense[,c(1:3)],glm.all.rolwin20.annot.spino.missense$Consequence,select(glm.all.rolwin20.annot.spino.missense,contains("fdr"))))

## Loop through each GLM-results FDR column to calculate spino candidate minimum values
glm.all.fdr.rolwin20.annot.spino.missense.min <- c() # create an object
#loop through all cols fdr values for candidate genes
for(g in 5:ncol(glm.all.fdr.rolwin20.annot.spino.missense)) { #nrow(haf.freq.cand)
  tryCatch({
  	#create 2-column table with gene name and one FDR column at a time
    temp.min <- as.data.frame(cbind(glm.all.fdr.rolwin20.annot.spino.missense$Gene,glm.all.fdr.rolwin20.annot.spino.missense[,g]))
    #name columns
    colnames(temp.min)[1] <- "Gene" 
    colnames(temp.min)[2] <- "fdr"    
    #assign minimum FDR values
    glm.all.fdr.rolwin20.annot.spino.missense.min <- cbind(glm.all.fdr.rolwin20.annot.spino.missense.min,
    	temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(fdr),
        	list(min = min)) %>%
        	pull(min))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.fdr.rolwin20.annot.spino.missense.min) <- names(glm.all.fdr.rolwin20.annot.spino.missense)[c(5:ncol(glm.all.fdr.rolwin20.annot.spino.missense))]
## Add gene names to the table
glm.all.fdr.rolwin20.annot.spino.missense.min <- cbind(Gene=temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(fdr),
        	list(min = min)) %>%
        	pull(Gene),glm.all.fdr.rolwin20.annot.spino.missense.min)

## Save minimum FDR results from missense SNPs for downstream use
write.table(glm.all.fdr.rolwin20.annot.spino.missense.min, file="rudflies_2023_redo.glm.all.rolwin20.minfdr.spino.missense.txt", quote = FALSE, sep = "\t", row.names = F)


## Pull FDR columns only for all genes/features, filtered for missense SNPs only
glm.all.fdr.rolwin20.annot.missense <- glm.all.fdr.rolwin20.annot[glm.all.fdr.rolwin20.annot$Consequence == "missense_variant",]

## Loop through each GLM-results FDR column to calculate minimum values for all genes
glm.all.fdr.rolwin20.annot.missense.min <- c() # create an object
#loop through all cols fdr values for candidate genes
for(g in 5:ncol(glm.all.fdr.rolwin20.annot.missense)) {
  tryCatch({
  	#create 2-column table with gene name and one FDR column at a time
    temp.min <- as.data.frame(cbind(glm.all.fdr.rolwin20.annot.missense$Gene,glm.all.fdr.rolwin20.annot.missense[,g]))
    #name columns
    colnames(temp.min)[1] <- "Gene" 
    colnames(temp.min)[2] <- "fdr"    
    #assign minimum FDR values
    glm.all.fdr.rolwin20.annot.missense.min <- cbind(glm.all.fdr.rolwin20.annot.missense.min,
    	temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(fdr),
        	list(min = min)) %>%
        	pull(min))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.fdr.rolwin20.annot.missense.min) <- names(glm.all.fdr.rolwin20.annot.missense)[c(5:ncol(glm.all.fdr.rolwin20.annot.missense))]
## Add gene names to the table
glm.all.fdr.rolwin20.annot.missense.min <- cbind(Gene=temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(fdr),
        	list(min = min)) %>%
        	pull(Gene),glm.all.fdr.rolwin20.annot.missense.min)


## Now we've found minimum FDR values for candidate gene and background gene lists, 
## let's run hypergeometric tests on the overlap of candidates and GLM outliers.

## run hypergeometric overlap test per GLM contrast, 20-SNP rolling windows
## First with all gene-related/adjacent SNPS
## Significance threshold FDR<0.05
phyper_all_contrasts05 <- c() #create empty object
q_all <- c() #candidate - outlier overlap gene count
m_all <- c() #candidate gene count
n_all <- c() #non-candidate gene count
k_all <- c() #outlier gene count
## Loop through all GLM columns
for(contrast in 2:ncol(glm.all.fdr.rolwin20.annot.spino.min)) {
  q <- length(unique(glm.all.fdr.rolwin20.annot.spino.min[glm.all.fdr.rolwin20.annot.spino.min[,contrast] < 0.05,1]))
  m <- length(unique(glm.all.fdr.rolwin20.annot.spino.min[,1]))
  n <- length(unique(glm.all.fdr.rolwin20.annot.min[,1])) - m
  k <- length(unique(glm.all.fdr.rolwin20.annot.min[glm.all.fdr.rolwin20.annot.min[,contrast] < 0.05,1]))
  ## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  phyper_all_contrasts05 <- append(phyper_all_contrasts05,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE))
  ## save results as we go
  q_all <- append(q_all,q)
  m_all <- append(m_all,m)
  n_all <- append(n_all,n)
  k_all <- append(k_all,k)    
}

## Create results summary table
phyper_spinoCand_FDR05 <- cbind(contrast=colnames(glm.all.fdr.rolwin20.annot.spino.min[,2:ncol(glm.all.fdr.rolwin20.annot.spino.min)]),phyper_pval=phyper_all_contrasts05,phyper_sig=phyper_all_contrasts05 < 0.05,q_all,m_all,k_all,n_all)

## Significance threshold FDR<0.01
phyper_all_contrasts01 <- c() #create empty object
q_all <- c() #candidate - outlier overlap gene count
m_all <- c() #candidate gene count
n_all <- c() #non-candidate gene count
k_all <- c() #outlier gene count
## Loop through all GLM columns
for(contrast in 2:ncol(glm.all.fdr.rolwin20.annot.spino.min)) {
  q <- length(unique(glm.all.fdr.rolwin20.annot.spino.min[glm.all.fdr.rolwin20.annot.spino.min[,contrast] < 0.01,1]))
  m <- length(unique(glm.all.fdr.rolwin20.annot.spino.min[,1]))
  n <- length(unique(glm.all.fdr.rolwin20.annot.min[,1])) - m
  k <- length(unique(glm.all.fdr.rolwin20.annot.min[glm.all.fdr.rolwin20.annot.min[,contrast] < 0.01,1]))
  ## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  phyper_all_contrasts01 <- append(phyper_all_contrasts01,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE))
  ## save results as we go
  q_all <- append(q_all,q)
  m_all <- append(m_all,m)
  n_all <- append(n_all,n)
  k_all <- append(k_all,k)  
}

## Create results summary table
phyper_spinoCand_FDR01 <- cbind(contrast=colnames(glm.all.fdr.rolwin20.annot.spino.min[,2:ncol(glm.all.fdr.rolwin20.annot.spino.min)]),phyper_pval=phyper_all_contrasts01,phyper_sig=phyper_all_contrasts01 < 0.05,q_all,m_all,k_all,n_all)


### Next with only missense SNPS ###

## Significance threshold FDR<0.05
phyper_missense_contrasts05 <- c() #create empty object
q_all <- c() #candidate - outlier overlap gene count
m_all <- c() #candidate gene count
n_all <- c() #non-candidate gene count
k_all <- c() #outlier gene count
## Loop through all GLM columns
for(contrast in 2:ncol(glm.all.fdr.rolwin20.annot.spino.missense.min)) {
  q <- length(unique(glm.all.fdr.rolwin20.annot.spino.missense.min[glm.all.fdr.rolwin20.annot.spino.missense.min[,contrast] < 0.05,1]))
  m <- length(unique(glm.all.fdr.rolwin20.annot.spino.missense.min[,1]))
  n <- length(unique(glm.all.fdr.rolwin20.annot.missense.min[,1])) - m
  k <- length(unique(glm.all.fdr.rolwin20.annot.missense.min[glm.all.fdr.rolwin20.annot.missense.min[,contrast] < 0.05,1]))
  ## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  phyper_missense_contrasts05 <- append(phyper_missense_contrasts05,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE))
  ## save results as we go
  q_all <- append(q_all,q)
  m_all <- append(m_all,m)
  n_all <- append(n_all,n)
  k_all <- append(k_all,k)  
}

## Create results summary table
phyper_missense_spinoCand_FDR05 <- cbind(contrast=colnames(glm.all.fdr.rolwin20.annot.spino.missense.min[,2:ncol(glm.all.fdr.rolwin20.annot.spino.missense.min)]),phyper_pval=phyper_missense_contrasts05,phyper_sig=phyper_missense_contrasts05 < 0.05,q_all,m_all,k_all,n_all)

## Significance threshold FDR<0.01
phyper_missense_contrasts01 <- c() #create empty object
q_all <- c() #candidate - outlier overlap gene count
m_all <- c() #candidate gene count
n_all <- c() #non-candidate gene count
k_all <- c() #outlier gene count
## Loop through all GLM columns
for(contrast in 2:ncol(glm.all.fdr.rolwin20.annot.spino.missense.min)) {
  q <- length(unique(glm.all.fdr.rolwin20.annot.spino.missense.min[glm.all.fdr.rolwin20.annot.spino.missense.min[,contrast] < 0.01,1]))
  m <- length(unique(glm.all.fdr.rolwin20.annot.spino.missense.min[,1]))
  n <- length(unique(glm.all.fdr.rolwin20.annot.missense.min[,1])) - m
  k <- length(unique(glm.all.fdr.rolwin20.annot.missense.min[glm.all.fdr.rolwin20.annot.missense.min[,contrast] < 0.01,1]))
  ## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  phyper_missense_contrasts01 <- append(phyper_missense_contrasts01,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE))
  ## save results as we go
  q_all <- append(q_all,q)
  m_all <- append(m_all,m)
  n_all <- append(n_all,n)
  k_all <- append(k_all,k)  
}

## Create results summary table
phyper_missense_spinoCand_FDR01 <- cbind(contrast=colnames(glm.all.fdr.rolwin20.annot.spino.missense.min[,2:ncol(glm.all.fdr.rolwin20.annot.spino.missense.min)]),phyper_pval=phyper_missense_contrasts01,phyper_sig=phyper_missense_contrasts01 < 0.05,q_all,m_all,k_all,n_all)

## Save results table for spinosyn candidates vs GLM outliers (FDR < 0.05)
write.table(phyper_spinoCand_FDR05, file="rudflies_2023.phyper_allSNPs_spinoCand_FDR05.txt", quote = FALSE, sep = "\t", row.names = F)

## Save results table for spinosyn candidates vs GLM outliers (FDR < 0.01)
write.table(phyper_spinoCand_FDR01, file="rudflies_2023.phyper_allSNPs_spinoCand_FDR01.txt", quote = FALSE, sep = "\t", row.names = F)

## Save results table for spinosyn candidates vs GLM outliers (FDR < 0.05), missense only
write.table(phyper_missense_spinoCand_FDR05, file="rudflies_2023.phyper_missenseSNPs_spinoCand_FDR05.txt", quote = FALSE, sep = "\t", row.names = F)

## Save results table for spinosyn candidates vs GLM outliers (FDR < 0.01), missense only
write.table(phyper_missense_spinoCand_FDR01, file="rudflies_2023.phyper_missenseSNPs_spinoCand_FDR01.txt", quote = FALSE, sep = "\t", row.names = F)

### Which contrasts show significant overlap, prior to p-value correction

## spino candidates vs GLM outliers (FDR < 0.05), all SNPs
phyper_spinoCand_FDR05[phyper_spinoCand_FDR05[,2]<0.05,]
## spino candidates vs GLM outliers (FDR < 0.01), all SNPs
phyper_spinoCand_FDR01[phyper_spinoCand_FDR01[,2]<0.05,]
## spino candidates vs GLM outliers (FDR < 0.05), missense SNPs
phyper_missense_spinoCand_FDR05[phyper_missense_spinoCand_FDR05[,2]<0.05,]
## spino candidates vs GLM outliers (FDR < 0.01), missense SNPs
phyper_missense_spinoCand_FDR01[phyper_missense_spinoCand_FDR01[,2]<0.05,]


#############################################################
### Examine candidate gene p-values via Wilcoxon rank sum ###
#############################################################

## Hypergeometric overlap analysis relies on GLM outlier cutoff decisions, which can 
## impact results. Wilcoxon Rank Sum tests are non-parametric and rely on rank order of
## all genes in the dataset. Therefore, this test is more sensitive to more modest 
## AF changes across a broader gene network or pathway.

### Remake annotation stat summary lists for pvalues instead of fdr values ###

## Pull raw p-value columns only for spino candidates
glm.all.pval.rolwin20.annot.spino <- unique(cbind(glm.all.rolwin20.annot.spino[,c(1:3)],Consequence=glm.all.rolwin20.annot.spino$Consequence,select(glm.all.rolwin20.annot.spino,contains("rolwin20"),-contains("fdr"),-contains("logp"))))

## Loop through each GLM-results pval column to calculate spino candidate minimum values
glm.all.pval.rolwin20.annot.spino.min <- c() # create an object
#loop through all cols p-values for candidate genes
for(g in 5:ncol(glm.all.pval.rolwin20.annot.spino)) {
  tryCatch({
  	#create 2-column table with gene name and one pval column at a time
    temp.min <- as.data.frame(cbind(glm.all.pval.rolwin20.annot.spino$Gene,glm.all.pval.rolwin20.annot.spino[,g]))
    #name columns
    colnames(temp.min)[1] <- "Gene" 
    colnames(temp.min)[2] <- "pval"
    #assign minimum p-values    
    glm.all.pval.rolwin20.annot.spino.min <- cbind(glm.all.pval.rolwin20.annot.spino.min,
    	temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(pval),
        	list(min = min)) %>%
        	pull(min))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.pval.rolwin20.annot.spino.min) <- names(glm.all.pval.rolwin20.annot.spino)[c(5:ncol(glm.all.pval.rolwin20.annot.spino))]
## Add gene names to the table
glm.all.pval.rolwin20.annot.spino.min <- cbind(Gene=temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(pval),
        	list(min = min)) %>%
        	pull(Gene),glm.all.pval.rolwin20.annot.spino.min)
        	
## Save minimum pval results for downstream use
write.table(glm.all.pval.rolwin20.annot.spino.min, file="rudflies_2023_redo.glm.all.rolwin20.minpval.spino.txt", quote = FALSE, sep = "\t", row.names = F)
        	
        	
## Pull pval columns only for all genes/features
glm.all.pval.rolwin20.annot <- unique(cbind(Gene=glm.all.rolwin20.annot$Gene,glm.all.rolwin20.annot[,c(1,2)],Consequence=glm.all.rolwin20.annot$Consequence,select(glm.all.rolwin20.annot,contains("rolwin20"),-contains("fdr"),-contains("logp"))))
glm.all.pval.rolwin20.annot <- na.omit(glm.all.pval.rolwin20.annot)

## Loop through each GLM-results pval column to calculate minimum values for all genes
glm.all.pval.rolwin20.annot.min <- c() # create an object
#loop through all cols p-values for candidate genes
for(g in 5:ncol(glm.all.pval.rolwin20.annot)) {
  tryCatch({
    #create 2-column table with gene name and one pval column at a time
    temp.min <- as.data.frame(cbind(glm.all.pval.rolwin20.annot$Gene,glm.all.pval.rolwin20.annot[,g]))
    #name columns
    colnames(temp.min)[1] <- "Gene" 
    colnames(temp.min)[2] <- "pval"  
    #assign minimum p-values  
    glm.all.pval.rolwin20.annot.min <- cbind(glm.all.pval.rolwin20.annot.min,
    	temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(pval),
        	list(min = min)) %>%
        	pull(min))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.pval.rolwin20.annot.min) <- names(glm.all.pval.rolwin20.annot)[c(5:ncol(glm.all.pval.rolwin20.annot))]
## Add gene names to the table
glm.all.pval.rolwin20.annot.min <- cbind(Gene=temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(pval),
        	list(min = min)) %>%
        	pull(Gene),glm.all.pval.rolwin20.annot.min)

## Calculate medians instead for plotting
glm.all.pval.rolwin20.annot.med <- unique(glm.all.pval.rolwin20.annot[,1])
glm.all.pval.rolwin20.annot.med <- glm.all.pval.rolwin20.annot.med[order(glm.all.pval.rolwin20.annot.med) ]
#loop through all cols p-values for candidate genes
for(g in c(3,5:ncol(glm.all.pval.rolwin20.annot))) {
  tryCatch({
    temp.med <- as.data.frame(cbind(glm.all.pval.rolwin20.annot$Gene,glm.all.pval.rolwin20.annot[,g]))
    colnames(temp.med)[1] <- "Gene" 
    colnames(temp.med)[2] <- "pval" 
    #assign median p-values
    glm.all.pval.rolwin20.annot.med <- cbind(glm.all.pval.rolwin20.annot.med,
    	temp.med %>%
			group_by(Gene) %>% 
			summarise(median = median(as.numeric(pval), na.rm = TRUE)) %>%
        	pull(median))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.pval.rolwin20.annot.med)[c(2:ncol(glm.all.pval.rolwin20.annot.med))] <- names(glm.all.pval.rolwin20.annot)[c(3,5:ncol(glm.all.pval.rolwin20.annot))]
colnames(glm.all.pval.rolwin20.annot.med)[1]="Gene"

## Add chromosome info to table for plotting
chrom.list <- unique(glm.all.pval.rolwin20.annot[,c(1:2)])
chrom.list <- chrom.list[order(chrom.list$Gene),]
chrom.list <- chrom.list[-c(1:5),]
glm.all.pval.rolwin20.annot.med <- merge(chrom.list, glm.all.pval.rolwin20.annot.med, by="Gene")
        	
        	
## Pull pval columns only for all genes/features, filtered for missense SNPs only
glm.all.pval.rolwin20.annot.spino.missense <- unique(cbind(glm.all.rolwin20.annot.spino.missense[,c(1:3)],glm.all.rolwin20.annot.spino.missense$Consequence,select(glm.all.rolwin20.annot.spino.missense,contains("rolwin20"),-contains("fdr"),-contains("logp"))))

## Loop through each GLM-results pval column to calculate minimum values for all genes
glm.all.pval.rolwin20.annot.spino.missense.min <- c() # create an object
#loop through all cols p-values for candidate genes
for(g in 5:ncol(glm.all.pval.rolwin20.annot.spino.missense)) { #nrow(haf.freq.cand)
  tryCatch({
	#create 2-column table with gene name and one pval column at a time
    temp.min <- as.data.frame(cbind(glm.all.pval.rolwin20.annot.spino.missense$Gene,glm.all.pval.rolwin20.annot.spino.missense[,g]))
    #name columns
    colnames(temp.min)[1] <- "Gene" 
    colnames(temp.min)[2] <- "pval"
    #assign minimum FDR values    
    glm.all.pval.rolwin20.annot.spino.missense.min <- cbind(glm.all.pval.rolwin20.annot.spino.missense.min,
    	temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(pval),
        	list(min = min)) %>%
        	pull(min))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.pval.rolwin20.annot.spino.missense.min) <- names(glm.all.pval.rolwin20.annot.spino.missense)[c(5:ncol(glm.all.pval.rolwin20.annot.spino.missense))]
## Add gene names to the table
glm.all.pval.rolwin20.annot.spino.missense.min <- cbind(Gene=temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(pval),
        	list(min = min)) %>%
        	pull(Gene),glm.all.pval.rolwin20.annot.spino.missense.min)


## Save minimum pval results from missense SNPs for downstream use
write.table(glm.all.pval.rolwin20.annot.spino.missense.min, file="rudflies_2023_redo.glm.all.rolwin20.minpval.spino.missense.txt", quote = FALSE, sep = "\t", row.names = F)


## Pull pval columns only for all genes/features, filtered for missense SNPs only
glm.all.pval.rolwin20.annot.missense <- glm.all.pval.rolwin20.annot[glm.all.pval.rolwin20.annot$Consequence == "missense_variant",]

## Loop through each GLM-results pval column to calculate minimum values for all genes
glm.all.pval.rolwin20.annot.missense.min <- c() # create an object
#loop through all cols p-values for candidate genes
for(g in 5:ncol(glm.all.pval.rolwin20.annot.missense)) { #nrow(haf.freq.cand)
  tryCatch({
  	#create 2-column table with gene name and one pval column at a time
    temp.min <- as.data.frame(cbind(glm.all.pval.rolwin20.annot.missense$Gene,glm.all.pval.rolwin20.annot.missense[,g]))
    #name columns
    colnames(temp.min)[1] <- "Gene" 
    colnames(temp.min)[2] <- "pval"   
    #assign minimum p-values 
    glm.all.pval.rolwin20.annot.missense.min <- cbind(glm.all.pval.rolwin20.annot.missense.min,
    	temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(pval),
        	list(min = min)) %>%
        	pull(min))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.pval.rolwin20.annot.missense.min) <- names(glm.all.pval.rolwin20.annot.missense)[c(5:ncol(glm.all.pval.rolwin20.annot.missense))]
## Add gene names to the table
glm.all.pval.rolwin20.annot.missense.min <- cbind(Gene=temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(pval),
        	list(min = min)) %>%
        	pull(Gene),glm.all.pval.rolwin20.annot.missense.min)
        	
        	


## Time to run the Wicoxon Rank Sum analyses. To do so, we will select one matched 
## background gene for every spinosad candidate gene. This matched set selection is 
## based on the following criteria:
## 	 -Gene must be similar length from start to stop (+/-25%)
## 	 -Gene must have similar number of SNPs (+/-25%)
## 	 -Gene must be the same type (protein-coding, ncRNA, etc)
## 	 -Gene must be on the same chromosome

## After all matching genes are discovered, we randomly select one match for each 
## candidate gene, then perform the Wilcoxon rank sum test between the candidate and 
## background genes to see if candidates show significantly higher p-value ranks (lower 
## p-values), or greater parallel allele frequency differences in the respective GLM 
## contrast. This process is repeated for 1K iterations of matched set selection to 
## mitigate any potential random sampling bias. 


##Find gene characteristics of candidates from GFF for matching
glm.all.pval.rolwin20.annot.spino.min.gff <- merge(snp.gff.overlap, glm.all.pval.rolwin20.annot.spino.min, by.x="ID", by.y="Gene")

##Find gene characteristics of background genes from GFF for matching
glm.all.pval.rolwin20.annot.min.gff <- merge(snp.gff.overlap, glm.all.pval.rolwin20.annot.min, by.x="ID", by.y="Gene")

## Pull spino candidate gene IDs
spino_IDs <- glm.all.pval.rolwin20.annot.spino.min.gff$ID

## Remove spino candidates from potential background gene list
glm.all.pval.rolwin20.annot.bg.min.gff <- glm.all.pval.rolwin20.annot.min.gff[!(glm.all.pval.rolwin20.annot.min.gff$ID %in% spino_IDs),]

## Create empty results table with 1K rows (for each iteration) and one column for 
## each GLM contrast.
wilcox.rolwin20.p <- data.frame(matrix(NA, nrow = 1000, ncol = ncol(glm.all.pval.rolwin20.annot.spino.min.gff)-ncol(snp.gff.overlap)))
## Modify column names to reflect GLM contrasts
names(wilcox.rolwin20.p) <- names(glm.all.pval.rolwin20.annot.spino.min.gff)[c(9:ncol(glm.all.pval.rolwin20.annot.spino.min.gff))]

#generate candidate gene filtering tables
for(k in 1:1000) { #set number of iterations
  #Establish background non-candidate set
  bg.samp.list <- c()
  #This loop finds a random set of matched background genes
  for(i in c(1:nrow(glm.all.pval.rolwin20.annot.spino.min.gff))) {
  	#Select matched lists based on following criteria
	length <- glm.all.pval.rolwin20.annot.spino.min.gff[i,]$LENGTH
	count <- glm.all.pval.rolwin20.annot.spino.min.gff[i,]$SNP_COUNT
	type <- glm.all.pval.rolwin20.annot.spino.min.gff[i,]$TYPE
	chrom <- glm.all.pval.rolwin20.annot.spino.min.gff[i,]$CHROM
	#Pull one match per candidate gene
	bg.samp.temp <- sample_n(glm.all.pval.rolwin20.annot.bg.min.gff %>%
	 filter( 
	    TYPE==type,
		CHROM==chrom,
		LENGTH < length*1.25,
		LENGTH > length*0.75,
		SNP_COUNT < count*1.25,
		SNP_COUNT > count*0.75
	  ),1)
	#append to background list
	bg.samp.list <- rbind(bg.samp.list,bg.samp.temp)
  }
  #This loop will run Wilcoxon Rank Sum tests for each of "j" GLM contrasts
  for(j in 1:(ncol(glm.all.pval.rolwin20.annot.spino.min.gff)-ncol(snp.gff.overlap))) { #set number of iterations to the number of contrasts
    #Establish Rank Sum test sets for p values
    #candidate set
    candidates <- glm.all.pval.rolwin20.annot.spino.min.gff[,8+j]  
    #background non-candidate set
    bg.samp <- as.numeric(bg.samp.list[,8+j])
    candidates.samp <- as.numeric(candidates)
    #this will subsample the SNPs contained in one list to the number of SNPs in the shorter list
    if(length(candidates.samp) > length(bg.samp)){
      candidates.samp <- sample(candidates.samp, length(bg.samp))
    } else {
      bg.samp <- sample(bg.samp, length(candidates.samp))
    }
    #finally, run the tests and append results
    wilcox.rolwin20.p[k,j] <- wilcox.test(candidates.samp, bg.samp, alternative = "less")$p.value
}
}

## Save all results
# Wilcoxon rank sum p-values for all iterations
write.table(wilcox.rolwin20.p, file="rudflies_2023_redo.glm.all.wilcox.rolwin20.p.txt", quote = FALSE, sep = "\t", row.names = F)

# Median Wilcoxon rank sum p-values across iterations
write.table(colMedians(as.matrix(wilcox.rolwin20.p)), file="rudflies_2023_redo.glm.all.wilcox.rolwin20.p.medians.txt", quote = FALSE, sep = "\t", row.names = F)

# Mean Wilcoxon rank sum p-values across iterations
write.table(colMeans(as.matrix(wilcox.rolwin20.p)), file="rudflies_2023_redo.glm.all.wilcox.rolwin20.p.means.txt", quote = FALSE, sep = "\t", row.names = F)

## Calculate the percent of p-values significant at a raw p<0.05
colMedians(as.matrix(wilcox.rolwin20.p))
colMeans(as.matrix(wilcox.rolwin20.p))
wilcox.rolwin20.contrast <- c()
wilcox.rolwin20.contrast.p05 <- wilcox.rolwin20.p<0.05
for(l in 1:ncol(wilcox.rolwin20.p)) {
   wilcox.rolwin20.contrast <- append(wilcox.rolwin20.contrast,sum(wilcox.rolwin20.contrast.p05[,l], na.rm=TRUE)*0.1)
}
cbind(contrast=names(wilcox.rolwin20.p),perc_sig=wilcox.rolwin20.contrast)

# Save percent of significant Wilcoxon rank sum p-values across iterations
write.table(cbind(contrast=names(wilcox.rolwin20.p),perc_sig=wilcox.rolwin20.contrast), file="rudflies_2023_redo.glm.all.wilcox.rolwin20.p.percsig.txt", quote = FALSE, sep = "\t", row.names = F)


