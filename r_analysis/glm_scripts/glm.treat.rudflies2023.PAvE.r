########################################################################
## This R script will run GLMs for all SNPs in Rudman 2023 sequences. ##
## Specifically, the model examines all the contrast between PA and E ## 
## for each time-point and extracts p-values with emmeans.            ##
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

#filter for founder samples for S and S
haf.meta.filt.PA.E <- haf.meta[haf.meta$batch == "a" & (haf.meta$treat.fix == "E" | haf.meta$treat.fix == "PA"),]
haf.meta.filt.PA.E <- haf.meta.filt.PA.E[which((haf.meta.filt.PA.E$samp %in% c(1:9,30:171))==TRUE),]
haf.freq.filt.PA.E <- haf.freq[, which((names(haf.freq) %in% haf.meta.filt.PA.E$samp)==TRUE)]
haf.meta.filt.PA.E <- haf.meta.filt.PA.E[which((haf.meta.filt.PA.E$samp %in% names(haf.freq.filt.PA.E))==TRUE),]

#additional filter to remove extreme low variance loci
haf.freq.filt.PA.E.loci <- na.omit(haf.freq[rowVars(as.matrix(haf.freq.filt.PA.E))>0.001,c(1:2)])
haf.freq.filt.PA.E <- na.omit(haf.freq.filt.PA.E[rowVars(as.matrix(haf.freq.filt.PA.E))>0.001,])

#tables need a bit of reformatting for glm to run
haf.meta.filt.PA.E$treat.fix <- as.factor(haf.meta.filt.PA.E$treat.fix)
haf.meta.filt.PA.E$tpt <- as.factor(haf.meta.filt.PA.E$tpt)

glmres.PA.E <- c() #initialize glm object
contrastout.treat.PA.E <- c() #initialize p-val summary table

#loop through all rows of allele frequency table
for(g in 1:nrow(haf.freq.filt.PA.E)) { #nrow(haf.freq.filt.PA.E)
  tryCatch({
    glmdf <- cbind(haf.meta.filt.PA.E,t(round(haf.freq.filt.PA.E[g,]*100)))
    colnames(glmdf)[13] <- "count"
    #sapply(glmdf, class)
    glmres.PA.E[[g]] <- glm(count~treat.fix * tpt, data = glmdf, family = quasipoisson, na.action = na.exclude)
    #summary(emtemp)$p.value
    contrastout.treat.PA.E <- rbind(contrastout.treat.PA.E,joint_tests(glmres.PA.E[[g]], by = "tpt")$p.value)
  }, error=function(e){})
}

#save results so we don't have to do this again
write.table(contrastout.treat.PA.E, file="rudflies_2023_PAvE_treat.GLMcontrast2.txt", sep="\t", quote = FALSE, row.names = F)

#append loci to table
contrastout.treat.PA.E.table <- cbind(haf.freq.filt.PA.E.loci,af.mean=rowMeans(haf.freq.filt.PA.E),contrastout.treat.PA.E)
write.table(contrastout.treat.PA.E.table, file="rudflies_2023_PAvE_treat.table.GLMcontrast2.txt", sep="\t", quote = FALSE, row.names = F)

