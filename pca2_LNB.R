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



# Ejemplo NSC-->LNB
class(counts)
dim(counts)
head(counts)
summary(counts)

dataset2<-counts[,c(7,8,9,4,5,6)]
head(dataset2)

df<- melt(dataset2)

ggplot(df, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot()+
  theme_classic()+
  labs(title="Read distribution")


highlyExpress<-rowMeans(dataset2)>50000
dataset2[which(highlyExpress),]


# quienes son los outliers ?ggplot(df, aes(x=variable, y=value, fill=variable)) + 
cpms<-cpm(dataset2)
dim(dataset2)

keep<-rowSums(cpms>1)>3
dim(dataset2)
countsf<-dataset2[keep,]
dim(countsf)# 35%
summary(countsf)
#####################
#primero hay que armar el objeto:
# Podemos replotear la dist:
dff<- melt(countsf)

ggplot(dff, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot()+
  theme_classic()+
  labs(title="Read distribution")

highlyExpressF<-rowMeans(countsf)>50000
countsf[which(highlyExpressF),]


condition<-factor(rep(c("NSC","LNB"), each=3))
levels(condition)

# The levels of a factor are re-ordered so that 
# the level specified by ref is first and the others are moved down. 
# This is useful for contr.treatment contrasts 
# which take the first level as the reference.

condition<-factor(rep(c("NSC","LNB"), each=3), 
                  levels = c("NSC","LNB"))
levels(condition)

head(dataset2)

y <- DGEList(counts=dataset2, group=condition)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
dim(y)
condition
# ojo porque se factoriza por orden alfabetico !
# 
design <- model.matrix(~ condition)
design

y <- calcNormFactors(y)
y
factor(condition)
plotMDS(y, labels=condition,
col=c("darkgreen","blue")[factor(condition)])


y <- estimateDisp(y, verbose=TRUE)
y <- estimateCommonDisp(y,verbose = T)
y <- estimateTagwiseDisp(y, verbose = T)

plotBCV(y)

de<-exactTest(y, pair = c("NSC","LNB"))
str(de)
head(de)

## ----toptagsDefault, echo=TRUE, eval=TRUE---------------------------
tt <- topTags(de)
tt

## ----toptags, echo=TRUE, eval=TRUE----------------------------------
tt <- topTags(de, n = nrow(de))

## ----toptags2, echo=TRUE, eval=TRUE---------------------------------
table(tt$table$FDR <0.05)

df <- data.frame(
  exp="LNB",
  fdr=tt$table$FDR)
# Change colors
ggplot(df, aes(x=fdr)) + 
geom_histogram(color="black", fill="white", binwidth=0.01)+
geom_vline(xintercept=0.05, linetype="dashed")+
theme_classic()+
labs(title="FDR distribution")

# ----deg, echo=TRUE, eval=TRUE--------------------------------------
deg<-rownames(tt)[tt$table$FDR <.05 &   
                  abs(tt$table$logFC )>1 ]
plotSmear(y, de.tags=deg)
abline(h=c(-1,0,1))

## ----degGGPLOT, echo=TRUE, eval=TRUE, warning=FALSE-----------------
forPlot <-tt$table
tt10 <- topTags(de, n=10)

forPlot$gene_color <- rep("grey", nrow(forPlot))
forPlot$gene_color[forPlot$logFC>1] <-"red"   
forPlot$gene_color[forPlot$logFC< (-1)]<-"green"
forPlot$imp_genes<-NA

ii <- match(rownames(tt10), rownames(forPlot))
forPlot$imp_genes[ii]<-rownames(forPlot)[ii]


ggplot(forPlot, aes(x=logFC, y=-log10(FDR))) +
  geom_point(aes(col=gene_color), cex= 1.2) +
  scale_color_manual(values=c("dark green","dark grey", "dark red")) +
  labs(title="DEG LNB", x="log2(FC)", y="-log10(FDR)") +
  geom_vline(xintercept= c(-1, 1), colour= 'black', linetype= 'dashed') +
  geom_hline(yintercept= 1.30103, colour= 'black', linetype= 'dashed') +
  theme_minimal()+
  theme(legend.position = "none",
        plot.title = element_text(size = 12, face="italic", hjust=0.4),
        axis.title.x = element_text(color = "black", size=12, hjust = 0.4),   
        axis.title.y = element_text(size =12, hjust = 0.5)) +
  geom_text_repel(data=forPlot,
                  aes(x=logFC, y=-log10(FDR)), 
                  label =forPlot$imp_genes,
                  box.padding = unit(0.25, "lines"),
                  hjust =1,
                  max.overlaps = 50)

## ----tt2, echo=TRUE, eval=FALSE-------------------------------------
write.csv(tt$table, file="LNB_edgeR.csv")
tt500 <- topTags(de, n =500)
write.csv(tt500$table, file="top500_LNB_edgeR.csv")

# symbol retrieving

# gene ontology
