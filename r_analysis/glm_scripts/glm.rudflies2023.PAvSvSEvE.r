########################################################################
## This R script will run GLMs for all SNPs in Rudman 2023 sequences. ##
## Specifically, the model examines all pairwise contrasts between    ## 
## PA, E, S, and S-extinct sample from T1 and extracts p-values.      ##
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

#filter for S T1 samples and append survival condition

filt.table <- cbind(cage=c("3","7","11","15","21","27","33","37","41","45","2","6","10","14","20","24","31","36","40","44","1","5","9","13","19","23","29","35","39","43","1","5","9","13","19","23","29","35","39","43"),tpt=c("1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1"),condition=c("SP","SP","SE","SE","SE","SE","SP","SP","SE","SE","PA","PA","PA","PA","PA","PA","PA","PA","PA","PA","E","E","E","E","E","E","E","E","E","E","E","E","E","E","E","E","E","E","E","E"))
haf.meta.filt <- merge(haf.meta, filt.table, by=c("cage","tpt"))

haf.meta.filt <- unique(haf.meta.filt[haf.meta.filt$batch == "a" & (haf.meta.filt$treat == "S" | haf.meta.filt$treat == "PA" | haf.meta.filt$treat == "E"),])
haf.meta.filt <- haf.meta.filt[which((haf.meta.filt$samp %in% c(1:9,30:171))==TRUE),]
haf.freq.filt <- haf.freq[, which((names(haf.freq) %in% haf.meta.filt$samp)==TRUE)]

#additional filter to remove extreme low variance loci
haf.freq.filt.loci <- na.omit(haf.freq[rowVars(as.matrix(haf.freq.filt))>0.001,c(1:2)])
haf.freq.filt <- na.omit(haf.freq.filt[rowVars(as.matrix(haf.freq.filt))>0.001,])

#tables need a bit of reformatting for glm to run
haf.meta.filt$treat <- as.factor(haf.meta.filt$treat)
haf.meta.filt$tpt <- as.factor(haf.meta.filt$tpt)
haf.meta.filt$condition <- as.factor(haf.meta.filt$condition)
haf.meta.filt <- haf.meta.filt[order(haf.meta.filt$samp),]

glmres <- c() #initialize glm object
contrastout <- c() #initialize p-val summary table

#loop through all rows of allele frequency table
for(g in 1:nrow(haf.freq.filt)) {
  tryCatch({
    glmdf <- cbind(haf.meta.filt,t(round(haf.freq.filt[g,]*100)))
    colnames(glmdf)[13] <- "count"
    #sapply(glmdf, class)
    glmres[[g]] <- glm(count~condition, data = glmdf, family = quasipoisson, na.action = na.exclude)
    emtemp <- emmeans(glmres[[g]], pairwise ~ condition)$contrasts
    #summary(emtemp)$p.value
    contrastout <- rbind(contrastout,summary(emtemp)$p.value)
  }, error=function(e){})
}

#save results so we don't have to do this again
write.table(contrastout, file="rudflies_2023_PAvSvSEvE.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)
write.table(cbind(haf.freq.filt.loci,contrastout), file="rudflies_2023_PAvSvSEvE.wLoci.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)
