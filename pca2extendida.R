library(ggplot2)
library(reshape2) 
library(edgeR)
library(ggrepel)

# NSC-->ENB-->LNB
## ----targets, echo=TRUE, eval=TRUE----------------------------------
counts <-
    read.table(
        "~/Documents/BioApps/data/mouse/STAR/geneCountsMmu.tab",
        row.names = 1,
        header = T,
        stringsAsFactors = F
    )


# Ejemplo NSC-->ENB
class(counts)
dim(counts)
head(counts)
summary(counts)

dataset1<-counts[,c(7,8,9,1,2,3)]
head(dataset1)
condition<-factor(rep(c("NSC","ENB"), each=3), 
                  levels = c("NSC","ENB"))
levels(condition)

y <- DGEList(counts=dataset1, group=condition)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
dim(y)
condition
design <- model.matrix(~ condition)
design

y <- calcNormFactors(y)
y

plotMDS(y, labels=condition,
col=c("darkgreen","blue")[factor(condition)])

y <- estimateCommonDisp(y,verbose = T)
y <- estimateTagwiseDisp(y, verbose = T)

y <- estimateDisp(y, verbose=TRUE)
plotBCV(y)

de<-exactTest(y)
str(de)
head(de)

## ----toptagsDefault, echo=TRUE, eval=TRUE---------------------------
tt <- topTags(de)
tt

## ----toptags, echo=TRUE, eval=TRUE----------------------------------
tt <- topTags(de, n = nrow(de))
tt1000 <- topTags(de, n=1000, sort.by = "logFC")
head(tt1000)
####################
countsTC<-counts[, c(7,8,9,1,2,3,4,5,6)]
colnames(countsTC)
conditionTC<-factor(rep(c("NSC","ENB", "LNB"), each=3), 
                  levels = c("NSC","ENB", "LNB"))
levels(conditionTC)
conditionTC

yTC <- DGEList(counts=countsTC, group=conditionTC)
keepTC <- filterByExpr(yTC)
summary(keepTC)
yTC <- yTC[keepTC, , keep.lib.sizes=FALSE]
yTC<-calcNormFactors(yTC)
yTC
designTC <- model.matrix(~yTC$samples$group)
rownames(designTC) <- colnames(yTC)
designTC

yTC <- estimateDisp(yTC)
yTC$common.dispersion
sqrt(yTC$common.dispersion)
summary(yTC$trended.dispersion)
summary(yTC$tagwise.dispersion)
#coeficiente de variación biológica
plotBCV(yTC)

plotMDS(yTC,
        labels = conditionTC,
        col = c("darkgreen", "blue", "red")[factor(conditionTC)])

fitTC<-glmFit(yTC)
head(fitTC$coefficients)
class(rowMeans(yTC$counts))
fitTC<-glmFit(yTC,designTC)
lrtENB<-glmLRT(fitTC, coef=2)
?topTags
head(lrtENB)
ENBtt1000 <- topTags(lrtENB, n=1000, sort.by = "logFC")
head(ENBtt1000)
###############COMPARACION#################

tti <- match(rownames(tt1000), rownames(ENBtt1000))
ttiCleanTC <- tti[!is.na(tti)]
length(ttiCleanTC)

ttiCleanPW<-which(!is.na(tti))
ttiCleanPW

toPlot <- data.frame(pwNames = 
                     rownames(tt1000$table)[ttiCleanPW],
                     pwFC = 
                     tt1000$table$logFC[ttiCleanPW],
                     pwFDR=tt1000$table$FDR[ttiCleanPW],
                     tcNames = 
                     rownames(ENBtt1000$table)[ttiCleanTC],
                     tcFC = ENBtt1000$table$logFC[ttiCleanTC],
                     tcFDR = ENBtt1000$table$FDR[ttiCleanTC]
                     )

plot(toPlot$pwFC, toPlot$tcFC)
abline(h=c(-3,0,3), v=c(-3,0,3))
plot(toPlot$pwFDR, toPlot$tcFDR)
