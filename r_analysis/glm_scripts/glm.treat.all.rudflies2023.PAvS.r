########################################################################
## This R script will run GLMs for all SNPs in Rudman 2023 sequences. ##
## Specifically, the model examines all the contrast between PA and S ## 
## across all time-points and extracts p-values with emmeans.         ##
########################################################################

#configure r environment
library(emmeans)
library(matrixStats)
setwd("/scratch/user/jfaberha/20251208_111954/admera/gp_analysis/rudflies_2023_redo/r")

#load required input files
haf.meta <- read.table("rudflies_2023_meta.tsv", header=TRUE)
haf.freq <- read.delim("rudflies_2023_hafpipe.csv", header=TRUE, sep = ",")
names(haf.freq) <- gsub("[.]af","",names(haf.freq))
names(haf.freq) <- gsub("X","",names(haf.freq))
names(haf.freq)[1]="CHROM"
names(haf.freq)[2]="POS"
names(haf.freq) <- gsub("[.]af","",names(haf.freq))
names(haf.freq) <- gsub("X","",names(haf.freq))

#filter for founder samples for S and PA
haf.meta.filt.PA.S <- haf.meta[haf.meta$batch == "a" & (haf.meta$treat.fix == "S" | haf.meta$treat.fix == "PA"),]
haf.meta.filt.PA.S <- haf.meta.filt.PA.S[which((haf.meta.filt.PA.S$samp %in% c(1:9,30:171))==TRUE),]
haf.freq.filt.PA.S <- haf.freq[, which((names(haf.freq) %in% haf.meta.filt.PA.S$samp)==TRUE)]
haf.meta.filt.PA.S <- haf.meta.filt.PA.S[which((haf.meta.filt.PA.S$samp %in% names(haf.freq.filt.PA.S))==TRUE),]

#additional filter to remove extreme low variance loci
haf.freq.filt.PA.S.loci <- na.omit(haf.freq[rowVars(as.matrix(haf.freq.filt.PA.S))>0.001,c(1:2)])
haf.freq.filt.PA.S <- na.omit(haf.freq.filt.PA.S[rowVars(as.matrix(haf.freq.filt.PA.S))>0.001,])

#tables need a bit of reformatting for glm to run
haf.meta.filt.PA.S$treat.fix <- as.factor(haf.meta.filt.PA.S$treat.fix)
haf.meta.filt.PA.S$tpt <- as.factor(haf.meta.filt.PA.S$tpt)

glmres.PA.S <- c() #initialize glm object
contrastout.treat.PA.S <- c() #initialize p-val summary table

#loop through all rows of allele frequency table
for(g in 1:nrow(haf.freq.filt.PA.S)) { #nrow(haf.freq.filt.PA.S)
  tryCatch({
    glmdf <- cbind(haf.meta.filt.PA.S,t(round(haf.freq.filt.PA.S[g,]*100)))
    colnames(glmdf)[12] <- "count"
    #sapply(glmdf, class)
    glmres.PA.S[[g]] <- glm(count~ treat.fix, data = glmdf, family = quasipoisson, na.action = na.exclude)
    #summary(emtemp)$p.value
    contrastout.treat.PA.S <- rbind(contrastout.treat.PA.S,as.data.frame(emmeans(glmres.PA.S[[g]], pairwise ~ treat.fix)$contrasts)$p.value)
  }, error=function(e){})
}

#save results so we don't have to do this again
write.table(contrastout.treat.PA.S, file="rudflies_2023_PAvS_treat.all.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)

#append loci to table
contrastout.treat.PA.S.table <- cbind(haf.freq.filt.PA.S.loci,af.mean=rowMeans(haf.freq.filt.PA.S),contrastout.treat.PA.S)
write.table(contrastout.treat.PA.S.table, file="rudflies_2023_PAvS_treat.all.table.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)
