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

#######################################
library(ASpli)
library(GenomicFeatures)

# gtf preprocessing ----
gtfFileName <- aspliExampleGTF()
genomeTxDb  <- makeTxDbFromGFF( gtfFileName )

# feature extraction ----
features    <- binGenome( genomeTxDb )
#bams and target file ----
BAMFiles <- aspliExampleBamList()
targets  <- data.frame(row.names = paste0('Sample',c(1:6)),
                       bam = BAMFiles[1:6],
                       f1  = c( 'control','control','control','treatment','treatment','treatment'),
                       stringsAsFactors = FALSE)
mBAMs <- data.frame( bam = sub("_[012]","",targets$bam[c(1,4)]),
                     condition = c("control","treatment"))

gbcounts <- gbCounts(features=features, targets=targets,
                     minReadLength = 100, maxISize = 50000)
gbcounts
head(gbcounts@gene.counts)
head(gbcounts@exon.intron.counts)


asd <- jCounts(counts=gbcounts, features=features, minReadLength=100)
head(asd@esPSI)
head(asd@irPIR)
head(asd@altPSI)



gb  <- gbDUreport(gbcounts, contrast = c(-1,1))
gb@genes
gb@bins
jdur <- jDUreport(asd, contrast=c(-1,1))
jdur

sr <- splicingReport(gb, jdur, counts=gbcounts)
is <- integrateSignals(sr,asd)
exportIntegratedSignals(is,sr=sr,
                        output.dir = "aspliExample",
                        counts=gbcounts,features=features,asd=asd,
                        mergedBams = mBAMs)

getwd()

head(gbcounts@junction.counts)
