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
head(rowMeans(counts))

counts<-counts[, c(7,8,9,1,2,3,4,5,6)]
colnames(counts)

df<- melt(counts)
df$condition<-df$variable
df$condition<-gsub("_1", "",df$condition)
df$condition<-gsub("_2", "",df$condition)
df$condition<-gsub("_3", "",df$condition)
df$condition <- factor(df$condition, 
                       levels = c("NSC","ENB", "LNB"))

ggplot(df, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot()+
  theme_classic()+
  labs(title="Read distribution")

ggplot(df, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot()+
  theme_classic()+
  labs(title="Read distribution")

plot(density(rowMeans(counts)))
hist(rowMeans(counts))

condition<-factor(rep(c("NSC","ENB", "LNB"), each=3), 
                  levels = c("NSC","ENB", "LNB"))
levels(condition)
condition

y <- DGEList(counts=counts, group=condition)
keep <- filterByExpr(y)
summary(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
y<-calcNormFactors(y)

design <- model.matrix(~y$samples$group)
rownames(design) <- colnames(y)
design

y <- estimateDisp(y)
y$common.dispersion
sqrt(y$common.dispersion)
summary(y$trended.dispersion)
summary(y$tagwise.dispersion)
#coeficiente de variación biológica
plotBCV(y)

plotMDS(y, labels=condition,
        col=c("darkgreen","blue","red")[factor(condition)])

fit<-glmFit(y)
head(fit$coefficients)
class(rowMeans(y$counts))
plotMDS(y, labels=condition,
        col=c("darkgreen","blue","red")[factor(condition)])
fit<-glmFit(y,design)
lrtENB<-glmLRT(fit, coef=2)
head(topTags(lrtENB))
summary(decideTests(lrtENB))

lrtLNB<-glmLRT(fit, coef=3)
head(topTags(lrtLNB))
summary(decideTests(lrtLNB))

lrtCT<-glmLRT(fit, coef=1)
head(topTags(lrtCT))
summary(decideTests(lrtCT))

lrt<-glmLRT(fit, coef=2:3); 
head(topTags(lrt))
summary(decideTests(lrt))

tt<-topTags(lrt, n = nrow(lrt$fitted.values))
head(tt)
topis<-which(tt$table$FDR < 0.00001); length(topis)
fc<-tt$table[topis,][,1:2]


paletteLength <- 50
myColor <- colorRampPalette(c("blue","white", "red"))(paletteLength)
myBreaks <- c(seq(min(fc), 0, 
                  length.out=ceiling(paletteLength/2) + 1), 
              seq(max(fc)/paletteLength, 
                  max(fc), 
                  length.out=floor(paletteLength/2)))

heat<-pheatmap(fc, 
               fontsize_col= 10,
               fontsize_row=1, 
               myColor, 
               breaks=myBreaks, 
               cluster_rows = TRUE,
               cutree_rows=10)

hc <-heat$tree_row
lbl <- cutree(hc, 10) 
table(lbl)
