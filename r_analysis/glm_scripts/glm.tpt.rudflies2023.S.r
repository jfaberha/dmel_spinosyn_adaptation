########################################################################
## This R script will run GLMs for all SNPs in Rudman 2023 sequences. ##
## Specifically, the model examines all pairwise contrasts between    ## 
## timepoints for S samples and extracts p-values with emmeans.       ##
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
haf.meta.filt.S <- haf.meta[haf.meta$batch == "a" & haf.meta$treat.fix == "S",]
haf.freq.filt.S <- haf.freq[, which((names(haf.freq) %in% haf.meta.filt.S$samp)==TRUE)]
haf.meta.filt.S <- haf.meta.filt.S[which((haf.meta.filt.S$samp %in% names(haf.freq.filt.S))==TRUE),]

#additional filter to remove extreme low variance loci
haf.freq.filt.S.loci <- na.omit(haf.freq[rowVars(as.matrix(haf.freq.filt.S))>0.001,c(1:2)])
haf.freq.filt.S <- na.omit(haf.freq.filt.S[rowVars(as.matrix(haf.freq.filt.S))>0.001,])

#tables need a bit of reformatting for glm to run
haf.meta.filt.S$treat.fix <- as.factor(haf.meta.filt.S$treat.fix)
haf.meta.filt.S$tpt <- as.factor(haf.meta.filt.S$tpt)

glmres.S <- c() #initialize glm object
contrastout.tpt.S <- c() #initialize p-val summary table

#loop through all rows of allele frequency table
for(g in 1:nrow(haf.freq.filt.S)) { #nrow(haf.freq.filt.S)
  tryCatch({
    glmdf <- cbind(haf.meta.filt.S,t(round(haf.freq.filt.S[g,]*100)))
    colnames(glmdf)[13] <- "count"
    #sapply(glmdf, class)
    glmres.S[[g]] <- glm(count~ tpt, data = glmdf, family = quasipoisson, na.action = na.exclude)
    #summary(emtemp)$p.value
    contrastout.tpt.S <- rbind(contrastout.tpt.S,as.data.frame(emmeans(glmres.S[[g]], pairwise ~ tpt)$contrasts)$p.value)
  }, error=function(e){})
}

#save results so we don't have to do this again
write.table(contrastout.tpt.S, file="rudflies_2023_S_tpt.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)

#append loci to table
contrastout.tpt.S.table <- cbind(haf.freq.filt.S.loci,af.mean=rowMeans(haf.freq.filt.S),contrastout.tpt.S)
write.table(contrastout.tpt.S.table, file="rudflies_2023_S_tpt.table.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)

