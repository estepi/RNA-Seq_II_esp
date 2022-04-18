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
?plotSmear
plotSmear(y, de.tags=deg)
abline(h=c(-1,0,1))

## ----degGGPLOT, echo=TRUE, eval=TRUE, warning=FALSE-----------------
ttable <-tt$table
dim(ttable)
tt10 <- topTags(de, n=20)

ttable$gene_color <- rep("grey", nrow(ttable))
ttable$gene_color[ttable$logFC>1] <-"red"   
ttable$gene_color[ttable$logFC< (-1)]<-"green"
ttable$imp_genes<-NA

ii <- match(rownames(tt10), rownames(ttable))
ttable$imp_genes[ii]<-rownames(ttable)[ii]

ggplot(ttable, aes(x=logFC, y=-log10(FDR))) +
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
  geom_text_repel(data=ttable,
                  aes(x=logFC, y=-log10(FDR)), 
                  label =ttable$imp_genes,
                  box.padding = unit(0.25, "lines"),
                  hjust =1,
                  max.overlaps = 50)

## ----tt2, echo=TRUE, eval=FALSE-------------------------------------
library(org.Mm.eg.db)
# primero tenemos que tener los ENTREX IDs
# nosotros tenemos los ENSEMBL Ids.
# usamos el paquete org.Mm.eg.db para hacer el mapeo
# org.Mm.egENSEMBL 
# Map Ensembl gene accession numbers with Entrez Gene identifiers
# Convert to a list
MmuENS <- toTable(org.Mm.egENSEMBL2EG)
head(MmuENS)
gi <- match(rownames(ttable), MmuENS$ensembl_id)
length(which(is.na(gi)))#983
ttable$ensemble_id <- rownames(ttable)
ttable$EntrezId <- MmuENS$gene_id[gi]
head(ttable)

MmuSYMBOL <- toTable(org.Mm.egSYMBOL)
head(MmuSYMBOL)
si <- match(ttable$EntrezId, MmuSYMBOL$gene_id)
length(which(is.na(si)))#983
ttable$Symbol <- MmuSYMBOL$symbol[si]
head(ttable)

# Podemos repetir ahora el plot un poco más entendible:
ii <- match(rownames(tt10), rownames(ttable))
ttable$imp_genes[ii]<-ttable$Symbol[ii]

ggplot(ttable, aes(x=logFC, y=-log10(FDR))) +
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
  geom_text_repel(data=ttable,
                  aes(x=logFC, y=-log10(FDR)), 
                  label =ttable$imp_genes,
                  box.padding = unit(0.25, "lines"),
                  hjust =1,
                  max.overlaps = 50)

###########################################
rownames(de)
entrezi <- match(rownames(de), MmuENS$ensembl_id)
length(which(is.na(entrezi)))#983
de$ensemble_id <- rownames(de)
de$EntrezId <- MmuENS$gene_id[entrezi]

go <- goana(de, species = "Mm", geneid=de$EntrezId)
topBPUp<-topGO(go, ont="BP", n=30, truncate=30, sort = "Up")
head(topBPUp)
topBPDown<-topGO(go, ont="BP", n=30, truncate=30, sort = "Down")
head(topBPDown)

topMF<-topGO(go, ont="MF", n=30, truncate=30)
head(topMF)
topCC<-topGO(go, ont="CC", n=30, truncate=30)
head(topCC)

keg <- kegga(de, species="Mm",  geneid=de$EntrezId)
topKEGG(keg)

