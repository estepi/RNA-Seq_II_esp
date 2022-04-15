library(ggplot2)
library(reshape2) 
library(edgeR)
library(pheatmap)
# NSC-->ENB-->LNB
## ----targets, echo=TRUE, eval=TRUE----------------------------------
counts <-
    read.table(
        "~/Documents/BioApps/data/mouse/STAR/geneCountsMmu.tab",
        row.names = 1,
        header = T,
        stringsAsFactors = F
    )

class(counts)
dim(counts)
head(counts)
summary(counts)
colnames(counts)

counts<-counts[, c(7,8,9,1,2,3,4,5,6)]

df<- melt(counts)
df$condition<-df$variable
df$condition<-gsub("_1", "",df$condition)
df$condition<-gsub("_2", "",df$condition)
df$condition<-gsub("_3", "",df$condition)
df$condition <- factor(df$condition, 
                       levels = c("NSC","ENB", "LNB"))

ggplot(df, aes(x=condition, y=value, fill=condition)) + 
  geom_boxplot()+
  theme_classic()+
  labs(title="Read distribution")


condition<-factor(rep(c("NSC","ENB", "LNB"), each=3), 
                  levels = c("NSC","ENB", "LNB"))
levels(condition)
condition
y <- DGEList(counts=counts, group=condition)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y<-calcNormFactors(y)

design <- model.matrix(~condition)
y <- estimateDisp(y, design, robust = T, verbose=TRUE)
plotBCV(y)

y$common.dispersion
y$common.dispersion^1/2

sqrt(0.159)

fit<-glmFit(y,design)

design

lrtENB<-glmLRT(fit, coef=2)
head(topTags(lrtENB))
summary(decideTests(lrtENB))

lrtLNB<-glmLRT(fit, coef=3)
head(lrtLNB)
summary(decideTests(lrtLNB))

lrtCT<-glmLRT(fit, coef=1)
head(topTags(lrtCT))

lrt<-glmLRT(fit, coef=2:3); 
head(topTags(lrt))
summary(decideTests(lrt))

#########################split the pipeline#########################
tt<-topTags(lrt, n = nrow(lrt$fitted.values))
head(tt)
tt$table$FDR
topis<-which(tt$table$FDR < 0.00001); length(topis)##856
fc<-tt$table[topis,][,1:2]
#colors
paletteLength <- 50
myColor <- colorRampPalette(c("blue","white", "red"))(paletteLength)

myBreaks <- c(seq(min(fc), 0, 
              length.out=ceiling(paletteLength/2) + 1), 
              seq(max(fc)/paletteLength, 
              max(fc), 
              length.out=floor(paletteLength/2)))
head(fc)
pheatmap(fc, 
         fontsize_col= 1,
         fontsize_row=2, 
         myColor, 
         breaks=myBreaks, 
         cluster_rows = TRUE,
         cutree_rows=10        )

heat <-
  pheatmap(fc, 
           cluster_rows = TRUE,  
           myColor           )
hc <-heat$tree_row
heat$tree_row
lbl <- cutree(hc, 10) # split gene dendrogram in 5 groups
table(lbl)
which(lbl==2) # grab genes of first group
#################################################
#podriamos seleccionar ahora los genes por grupos para analizar cosas concretas cmo GE, motivos de union TF, etc, etc