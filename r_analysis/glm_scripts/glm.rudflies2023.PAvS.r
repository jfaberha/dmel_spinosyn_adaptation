########################################################################
## This R script will run GLMs for all SNPs in Rudman 2023 sequences. ##
## Specifically, the model examines all pairwise contrasts between    ## 
## PA and S samples and extracts p-values with emmeans.               ##
########################################################################

#configure r environment
library(emmeans)
library(matrixStats)
setwd("/scratch/user/jfaberha/20251020_145408/admera/gp_analysis/rudflies_2023_redo/r")

#load required input files
haf.meta <- read.table("rudflies_2023_meta.tsv", header=TRUE)
haf.freq <- read.delim("rudflies_2023_hafpipe.csv", header=TRUE, sep = ",")
names(haf.freq) <- gsub("[.]af","",names(haf.freq))
names(haf.freq) <- gsub("X","",names(haf.freq))
names(haf.freq)[1]="CHROM"
names(haf.freq)[2]="POS"
names(haf.freq) <- gsub("[.]af","",names(haf.freq))
names(haf.freq) <- gsub("X","",names(haf.freq))

#filter for founder samples for PA and S
haf.meta.filt <- haf.meta[haf.meta$batch == "a" & (haf.meta$treat.fix == "PA" | haf.meta$treat.fix == "S"),]
haf.freq.filt <- haf.freq[, which((names(haf.freq) %in% haf.meta.filt$samp)==TRUE)]
haf.meta.filt <- haf.meta.filt[which((haf.meta.filt$samp %in% names(haf.freq.filt))==TRUE),]

#additional filter to remove extreme low variance loci
haf.freq.filt.loci <- na.omit(haf.freq.filt[rowVars(as.matrix(haf.freq.filt))>0.001,c(1,2)])
haf.freq.filt <- na.omit(haf.freq.filt[rowVars(as.matrix(haf.freq.filt))>0.001,])

#tables need a bit of reformatting for glm to run
haf.meta.filt$treat.fix <- as.factor(haf.meta.filt$treat.fix)
haf.meta.filt$tpt <- as.factor(haf.meta.filt$tpt)

glmres <- c() #initialize glm object
contrastout.treat <- c() #initialize p-val summary table
contrastout.tpt <- c() #initialize p-val summary table

#loop through all rows of allele frequency table
for(g in 1:nrow(haf.freq.filt)) { 
  tryCatch({
    glmdf <- cbind(haf.meta.filt,t(round(haf.freq.filt[g,]*100)))
    colnames(glmdf)[13] <- "count"
    #sapply(glmdf, class)
    glmres[[g]] <- glm(count~treat.fix + tpt, data = glmdf, family = quasipoisson, na.action = na.exclude)
    #summary(emtemp)$p.value
    contrastout.treat <- rbind(contrastout.treat,joint_tests(glmres[[g]], by = "tpt")$p.value)
    contrastout.tpt <- rbind(contrastout.tpt,as.data.frame(emmeans(glmres[[g]], pairwise ~ tpt)$contrasts)$p.value)
  }, error=function(e){})
}

#save results so we don't have to do this again
write.table(contrastout.treat, file="rudflies_2023_PAvS_treat.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)
write.table(contrastout.tpt, file="rudflies_2023_PAvS_tpt.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)

#append loci to table
contrast.treat.table <- cbind(haf.freq[rowVars(as.matrix(haf.freq.filt))>0.001,c(1:2)],af.mean=rowMeans(haf.freq.filt),contrastout.treat)
contrast.tpt.table <- cbind(haf.freq[rowVars(as.matrix(haf.freq.filt))>0.001,c(1:2)],af.mean=rowMeans(haf.freq.filt),contrastout.tpt)
write.table(contrastout.treat.table, file="rudflies_2023_PAvS_treat.table.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)
write.table(contrastout.tpt.table, file="rudflies_2023_PAvS_tpt.table.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)
