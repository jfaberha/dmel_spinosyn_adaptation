########################################################################
## This R script will run GLMs for all SNPs in Rudman 2023 sequences. ##
## Specifically, the model examines all the contrast between S and E ## 
## across all time-points and extracts p-values with emmeans.         ##
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
haf.meta.filt.S.E <- haf.meta[haf.meta$batch == "a" & (haf.meta$treat.fix == "E" | haf.meta$treat.fix == "S"),]
haf.meta.filt.S.E <- haf.meta.filt.S.E[which((haf.meta.filt.S.E$samp %in% c(1:9,30:171))==TRUE),]
haf.freq.filt.S.E <- haf.freq[, which((names(haf.freq) %in% haf.meta.filt.S.E$samp)==TRUE)]
haf.meta.filt.S.E <- haf.meta.filt.S.E[which((haf.meta.filt.S.E$samp %in% names(haf.freq.filt.S.E))==TRUE),]

#additional filter to remove extreme low variance loci
haf.freq.filt.S.E.loci <- na.omit(haf.freq[rowVars(as.matrix(haf.freq.filt.S.E))>0.001,c(1:2)])
haf.freq.filt.S.E <- na.omit(haf.freq.filt.S.E[rowVars(as.matrix(haf.freq.filt.S.E))>0.001,])

#tables need a bit of reformatting for glm to run
haf.meta.filt.S.E$treat.fix <- as.factor(haf.meta.filt.S.E$treat.fix)
haf.meta.filt.S.E$tpt <- as.factor(haf.meta.filt.S.E$tpt)

glmres.S.E <- c() #initialize glm object
#contrastout.treat.S.E <- c() #initialize p-val summary table
contrastout.int.S.E <- c() #initialize p-val summary table

#loop through all rows of allele frequency table
for(g in 1:nrow(haf.freq.filt.S.E)) { #nrow(haf.freq.filt.S.E)
  tryCatch({
    glmdf <- cbind(haf.meta.filt.S.E,t(round(haf.freq.filt.S.E[g,]*100)))
    colnames(glmdf)[13] <- "count"
    #sapply(glmdf, class)
    glmres.S.E[[g]] <- glm(count~treat.fix * tpt, data = glmdf, family = quasipoisson, na.action = na.exclude)
    #summary(emtemp)$p.value
    #contrastout.treat.S.E <- rbind(contrastout.treat.S.E,joint_tests(glmres.S.E[[g]], by = "tpt")$p.value)
	emm <- emmeans(glmres.S.E[[g]],~tpt*treat.fix) 
    contrastout.int.S.E <- rbind(contrastout.int.S.E,as.data.frame(contrast(emm, interaction = c("pairwise", "consec"), by = NULL))$p.value)
  }, error=function(e){})
}

#save results so we don't have to do this again
#write.table(contrastout.treat.S.E, file="rudflies_2023_SvE_treat.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)
write.table(contrastout.int.S.E, file="rudflies_2023_SvE_int.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)

#append loci to table
#contrastout.treat.S.E.table <- cbind(haf.freq.filt.S.E.loci,af.mean=rowMeans(haf.freq.filt.S.E),contrastout.treat.S.E)
#write.table(contrastout.treat.S.E.table, file="rudflies_2023_SvE_treat.table.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)
contrastout.int.S.E.table <- cbind(haf.freq.filt.S.E.loci,af.mean=rowMeans(haf.freq.filt.S.E),contrastout.int.S.E)
write.table(contrastout.int.S.E.table, file="rudflies_2023_SvE_int.table.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)
