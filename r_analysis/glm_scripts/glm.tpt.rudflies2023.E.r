########################################################################
## This R script will run GLMs for all SNPs in Rudman 2023 sequences. ##
## Specifically, the model examines all pairwise contrasts between    ## 
## timepoints for E samples and extracts p-values with emmeans.       ##
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

#filter for founder samples for S and S
haf.meta.filt.E <- haf.meta[haf.meta$batch == "a" & haf.meta$treat.fix == "E",]
haf.meta.filt.E <- haf.meta.filt.E[which((haf.meta.filt.E$samp %in% c(1:9,30:171))==TRUE),]
haf.freq.filt.E <- haf.freq[, which((names(haf.freq) %in% haf.meta.filt.E$samp)==TRUE)]
haf.meta.filt.E <- haf.meta.filt.E[which((haf.meta.filt.E$samp %in% names(haf.freq.filt.E))==TRUE),]

#additional filter to remove extreme low variance loci
haf.freq.filt.E.loci <- na.omit(haf.freq[rowVars(as.matrix(haf.freq.filt.E))>0.001,c(1:2)])
haf.freq.filt.E <- na.omit(haf.freq.filt.E[rowVars(as.matrix(haf.freq.filt.E))>0.001,])

#tables need a bit of reformatting for glm to run
haf.meta.filt.E$treat.fix <- as.factor(haf.meta.filt.E$treat.fix)
haf.meta.filt.E$tpt <- as.factor(haf.meta.filt.E$tpt)

glmres.E <- c() #initialize glm object
contrastout.tpt.E <- c() #initialize p-val summary table

#loop through all rows of allele frequency table
for(g in 1:nrow(haf.freq.filt.E)) { #nrow(haf.freq.filt.E)
  tryCatch({
    glmdf <- cbind(haf.meta.filt.E,t(round(haf.freq.filt.E[g,]*100)))
    colnames(glmdf)[13] <- "count"
    #sapply(glmdf, class)
    glmres.E[[g]] <- glm(count~ tpt, data = glmdf, family = quasipoisson, na.action = na.exclude)
    #summary(emtemp)$p.value
    contrastout.tpt.E <- rbind(contrastout.tpt.E,as.data.frame(emmeans(glmres.E[[g]], pairwise ~ tpt)$contrasts)$p.value)
  }, error=function(e){})
}

#save results so we don't have to do this again
write.table(contrastout.tpt.E, file="rudflies_2023_E_tpt.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)

#append loci to table
contrastout.tpt.E.table <- cbind(haf.freq.filt.E.loci,af.mean=rowMeans(haf.freq.filt.E),contrastout.tpt.E)
write.table(contrastout.tpt.E.table, file="rudflies_2023_E_tpt.table.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)

