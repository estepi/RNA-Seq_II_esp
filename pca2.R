library(ggplot2)
library(reshape2) 
library(edgeR)

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

df<- melt(dataset1)

ggplot(df, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot()+
  theme_classic()+
  labs(title="Read distribution")


highlyExpress<-rowMeans(dataset1)>50000
dataset1[which(highlyExpress),]


# quienes son los outliers ?ggplot(df, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot()+
  theme_classic()+
  labs(title="Read distribution")+ 
  ylim(0, 100000)


cpms<-cpm(dataset1)

dim(dataset1)

keep<-rowSums(cpms>1)>3

dim(dataset1)
countsf<-dataset1[keep,]
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


condition<-factor(rep(c("NSC","ENB"), each=3))
levels(condition)

# The levels of a factor are re-ordered so that 
# the level specified by ref is first and the others are moved down. 
# This is useful for contr.treatment contrasts 
# which take the first level as the reference.

condition<-factor(rep(c("NSC","ENB"), each=3), 
                  levels = c("NSC","ENB"))
levels(condition)

y <- DGEList(counts=dataset1, group=condition)
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

plotMDS(y, labels=condition,
col=c("darkgreen","blue")[factor(condition)])


y <- estimateDisp(y, verbose=TRUE)
y <- estimateCommonDisp(y,verbose = T)
y <- estimateTagwiseDisp(y, verbose = T)

plotBCV(y)

de<-exactTest(y, pair = c("NSC","ENB"))
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
  exp="ENB",
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
y <-tt$table
tt10 <- topTags(de, n=20)

y$gene_color <- rep("grey", nrow(y))
y$gene_color[y$logFC>1] <-"red"   
y$gene_color[y$logFC< (-1)]<-"green"
y$imp_genes<-NA

ii <- match(rownames(tt10), rownames(y))
y$imp_genes[ii]<-rownames(y)[ii]

library(ggrepel)
ggplot(y, aes(x=logFC, y=-log10(FDR))) +
  geom_point(aes(col=gene_color), cex= 1.2) +
  scale_color_manual(values=c("dark green","dark grey", "dark red")) +
  labs(title="DEG ENB", x="log2(FC)", y="-log10(FDR)") +
  geom_vline(xintercept= c(-1, 1), colour= 'black', linetype= 'dashed') +
  geom_hline(yintercept= 1.30103, colour= 'black', linetype= 'dashed') +
  theme_minimal()+
  theme(legend.position = "none",
        plot.title = element_text(size = 12, face="italic", hjust=0.4),
        axis.title.x = element_text(color = "black", size=12, hjust = 0.4),   
        axis.title.y = element_text(size =12, hjust = 0.5)) +
  geom_text_repel(data=y,
                  aes(x=logFC, y=-log10(FDR)), 
                  label =y$imp_genes,
                  box.padding = unit(0.25, "lines"),
                  hjust =1,
                  max.overlaps = 50)

## ----tt2, echo=TRUE, eval=FALSE-------------------------------------
 write.csv(tt$table, file="ENB_edgeR.csv")

# symbol retrieving

# gene ontology
