########################################################################
## This R script will run GLMs for all SNPs in Rudman 2023 sequences. ##
## Specifically, the model examines all pairwise contrasts between    ## 
## timepoints for PA samples and extracts p-values with emmeans.      ##
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
haf.meta.filt.PA <- haf.meta[haf.meta$batch == "a" & haf.meta$treat.fix == "PA",]
haf.freq.filt.PA <- haf.freq[, which((names(haf.freq) %in% haf.meta.filt.PA$samp)==TRUE)]
haf.meta.filt.PA <- haf.meta.filt.PA[which((haf.meta.filt.PA$samp %in% names(haf.freq.filt.PA))==TRUE),]

#additional filter to remove extreme low variance loci
haf.freq.filt.PA.loci <- na.omit(haf.freq[rowVars(as.matrix(haf.freq.filt.PA))>0.001,c(1:2)])
haf.freq.filt.PA <- na.omit(haf.freq.filt.PA[rowVars(as.matrix(haf.freq.filt.PA))>0.001,])

#tables need a bit of reformatting for glm to run
haf.meta.filt.PA$treat.fix <- as.factor(haf.meta.filt.PA$treat.fix)
haf.meta.filt.PA$tpt <- as.factor(haf.meta.filt.PA$tpt)

glmres.PA <- c() #initialize glm object
contrastout.tpt.PA <- c() #initialize p-val summary table

#loop through all rows of allele frequency table
for(g in 1:nrow(haf.freq.filt.PA)) { #nrow(haf.freq.filt.PA)
  tryCatch({
    glmdf <- cbind(haf.meta.filt.PA,t(round(haf.freq.filt.PA[g,]*100)))
    colnames(glmdf)[13] <- "count"
    #sapply(glmdf, class)
    glmres.PA[[g]] <- glm(count~ tpt, data = glmdf, family = quasipoisson, na.action = na.exclude)
    #summary(emtemp)$p.value
    contrastout.tpt.PA <- rbind(contrastout.tpt.PA,as.data.frame(emmeans(glmres.PA[[g]], pairwise ~ tpt)$contrasts)$p.value)
  }, error=function(e){})
}

#save results so we don't have to do this again
write.table(contrastout.tpt.PA, file="rudflies_2023_PA_tpt.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)

#append loci to table
contrastout.tpt.PA.table <- cbind(haf.freq.filt.PA.loci,af.mean=rowMeans(haf.freq.filt.PA),contrastout.tpt.PA)
write.table(contrastout.tpt.PA.table, file="rudflies_2023_PA_tpt.table.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)

