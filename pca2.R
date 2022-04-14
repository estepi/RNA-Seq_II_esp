library(ggplot2)
library(reshape2) 
library(edgeR)

# NSC-->ENB-->LNB
## ----targets, echo=TRUE, eval=TRUE----------------------------------
targets <- read.table("~/Documents/BioApps/data/mouse/targets.txt", header = T, row.names = 1)

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


ggplot(df, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot()+
  theme_classic()+
  labs(title="Read distribution")+ 
  ylim(0, 150000)

# quienes son los outliers ?
# podriamos pedir reads > 50?
# o replotear con un eje Y diferente


## ----filtros, echo=TRUE, eval=TRUE----------------------------------
cpms<-cpm(counts[,1:ncol(counts)])


## ----keep, echo=TRUE, eval=TRUE-------------------------------------
keep<-rowSums(cpms>1)>4


## ----keep2, echo=TRUE, eval=TRUE------------------------------------
dim(counts)
countsf<-counts[keep,]
dim(countsf)


## ----keep3, echo=TRUE, eval=TRUE------------------------------------
summary(countsf)


## ----DGEList1, echo=TRUE, eval=TRUE---------------------------------
d<-DGEList(counts=countsf, group=targets$condition)
str(d)


## ----calcNormF, echo=TRUE, eval=TRUE--------------------------------
d<-calcNormFactors(d)
d


## ----names, echo=TRUE, eval=TRUE------------------------------------
shortNames<-paste(targets$condition, rep(1:4, 2), sep=".")
targets<-cbind(targets,shortNames)

plotMDS(d, labels=targets$shortNames,
col=c("darkgreen","blue")[factor(targets$condition)])



## ----names2, echo=TRUE, eval=TRUE-----------------------------------
d<-estimateCommonDisp(d, verbose=TRUE)


## ----disp, echo=TRUE, eval=TRUE-------------------------------------
d<-estimateTagwiseDisp(d)
d


## ----plotBCV, echo=TRUE, eval=TRUE----------------------------------
plotBCV(d)


## ----exactTest, echo=TRUE, eval=TRUE--------------------------------
de<-exactTest(d, pair=c("CT","KD"))
str(de)



## ----toptagsDefault, echo=TRUE, eval=TRUE---------------------------
tt <- topTags(de)
tt


## ----toptags, echo=TRUE, eval=TRUE----------------------------------
tt <- topTags(de, n = nrow(de))


## ----toptags2, echo=TRUE, eval=TRUE---------------------------------
table(tt$table$FDR <0.05)

df <- data.frame(
  exp="KD",
  fdr=tt$table$FDR)
# Change colors
ggplot(df, aes(x=fdr)) + 
geom_histogram(color="black", fill="white", binwidth=0.01)+
geom_vline(xintercept=0.05, linetype="dashed")+
theme_classic()+
labs(title="FDR distribution")



## ----deg, echo=TRUE, eval=TRUE--------------------------------------
deg<-rownames(tt)[tt$table$FDR <.05 &   
                  abs(tt$table$logFC )>1 ]
plotSmear(d, de.tags=deg)
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
  labs(title="DEG", x="log2(FC)", y="-log10(FDR)") +
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
## write.table() o write.csv():
## write.csv(tt$table, file="red_edgeR.csv")

