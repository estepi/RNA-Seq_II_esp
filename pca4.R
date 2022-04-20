
library(GenomicFeatures)

TxDb <- makeTxDbFromGFF(
  file="~/Documents/BioApps/ssh/chr14.gtf",
  format="gtf")

library(ASpli)
features <- binGenome(TxDb) 

library(RNAseqData.HNRNPC.bam.chr14)
genome <- TxDb
targets <- data.frame(bam=RNAseqData.HNRNPC.bam.chr14_BAMFILES,
                      condition=c(rep("CT",4),rep("KD",4)))

gbcounts <- gbCounts(features=features, targets=targets,
                     minReadLength = 100, maxISize = 50000)
asd <- jCounts(counts=gbcounts, features=features, minReadLength=100)
gb <-  gbDUreport(gbcounts, contrast = c(-1,1))
jdur <- jDUreport(asd, contrast=c(-1,1), strongFilter = FALSE)
jdur

###########################
mBAMs <- data.frame( bam = sub("_[012]","",targets$bam[c(1,4)]),
                       condition = c("CT","KD"))
mBAMs
###########################
sr <- splicingReport(gb, jdur, counts=gbcounts)
is <- integrateSignals(sr,asd)
exportIntegratedSignals(is,sr=sr,
                          output.dir = "aspliExample",
                          counts=gbcounts,features=features,asd=asd,
                          mergedBams = mBAMs)