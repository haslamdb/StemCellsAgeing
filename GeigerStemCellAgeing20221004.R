

setwd("C:/Users/dbhas/OneDrive/Documents/Code/Metagenomics/GeigerData/StemCellsAgeing")
setwd("C:/Users/HASI9S/OneDrive/Documents/Code/Metagenomics/GeigerData/StemCellsAgeing")
setwd("~/Documents/Code/Metagenomics/GeigerData/StemCellsAgeing")
# setwd("C:/Users/HASI9S/Documents/Code/Metagenomics/GeigerData/StemCellsAgeing")



list.of.packages<-c("reshape2","stringdist","stringr","plyr","vegan","labdsv","pvclust","gplots","RColorBrewer",
                    "Heatplus","plyr","reshape2","fossil","ade4","scales",
                    "ggbiplot","MASS","ggthemes","tidyverse","pheatmap","magrittr","readxl","robustbase","cowplot", "sda","locfdr",
                    "FactoMineR","factoextra","dunn.test","FSA","NBZIMM","DataCombine","gtools","heatmap3","glmmTMB","DescTools","CompQuadForm",
                    "dirmult","ecodist","GUniFrac","lme4","MASS","Matrix","permute", "GMPR", "bbmle", "VennDiagram", "vegan"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# devtools::install_version("ggplot2", "3.3.5")

lapply(list.of.packages, require, character.only = TRUE)
library(ggplot2)



# Set some RAG1-Old_HSCs colors
metaphlan.colors <- colorRampPalette(c("#000033",  "#007FFF", "cyan", "red",   "yellow"))
scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(50)
col=c("gray32", "royalblue4","firebrick")
diverging.colors<-colorRampPalette(c("#000033", "#4575b4","#74add1","#abd9e9","#e0f3f8","#fee090","#fdae61", "#f46d43","#d73027", "red"))
tableau.colors<-c("#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#59a14f", "#edc948", "#b07aa1",
                  "#ff9da7", "#9c755f", "#bab0ac")

## Set up some parameters and functions here

# Set RAG1-Old_HSCs parameters for ggplot graphs:
params<- function(x){theme(axis.text.x= element_text(size= 14, color="black")) +
    theme(axis.text.y = element_text(size= 14, color="black")) +
    theme(plot.title = element_text(size= 22, color="black")) +
    theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=18)) +
    theme(legend.title = element_text(size=18)) +
    theme(legend.text = element_text(size = 14))
}

paramsAngled <- function(x){theme(axis.text.x= element_text(size= 12, , angle=45, vjust = 1, hjust= 1, color="black")) +
    theme(axis.text.y = element_text(size= 14, color="black")) +
    theme(plot.title = element_text(size= 22, color="black")) +
    theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=18)) +
    theme(legend.title = element_text(size=18)) +
    theme(legend.position="none") +
    theme(legend.text = element_text(size = 14))
}


paramsBox<- function(x){theme(axis.text.x= element_text(size= 14, color="black")) +
    theme(axis.text.y = element_text(size= 18, color="black")) +
    theme(plot.title = element_text(size= 22, color="black")) +
    theme(axis.title.x = element_text(size=20), axis.title.y = element_text(size=22)) +
    theme(legend.title = element_text(size=18)) +
    theme(legend.text = element_text(size = 14)) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(size= 22, color="black")) +
    theme(axis.text.y = element_text(size= 16, angle=0, color="black"))
}

## Regression plot

ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

# Set some RAG1-Old_HSCs colors
metaphlan.colors <- colorRampPalette(c("#000033",  "#007FFF", "cyan", "red",   "yellow"))
scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(50)
col=c("gray32", "royalblue4","firebrick")
diverging.colors<-colorRampPalette(c("#000033", "#4575b4","#74add1","#abd9e9","#e0f3f8","#fee090","#fdae61", "#f46d43","#d73027", "red"))

# PAM clustering function

pam.clustering<- function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}

# JSD distance calculator - from  http://enterotype.embl.de/enterotypes.html
dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) {
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix)
}


# Noise removal function - note this can be increased to focus on the organisms
# most relevant to distance and PCoA calculation (e.g. set to 1% to get a few species)
noise.removal <- function(dataframe, percent=0.001, top=NULL){
  dataframe->Matrix
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent
  Matrix_1 <- Matrix[bigones,]
  print(percent)
  return(Matrix_1)
}



compute_ef_Background <- function(d, g, min_shrink = 0.3){
  
  raw_scores <-  sda.ranking(d,  g, diagonal = FALSE, verbose = TRUE,
                             fdr = FALSE)
  
  # undo scaling by class frequencies in order to get effect sizes
  g_summaries <- sda:::pvt.groups(g)
  freqs <- freqs.shrink(g_summaries$samples)
  m <- sqrt((1 - freqs)/freqs/length(g))
  
  # compute effect sizes
  # notice that the "cat.Yes" and "cat.No" have to be changed based on the factor
  ef <- raw_scores[, "cat.C57BL"]*m["C57BL"] - raw_scores[, "cat.RAG1"]*m["RAG1"]
  
  # compute tests statistics, first, we need the scaling factor
  # for the cat score between the two classes
  n0 <- sum(g_summaries$idx[, "C57BL"])
  n1 <- sum(g_summaries$idx[, "RAG1"])
  
  m_stats <- 1/sqrt(1/n0 + 1/n1)
  stats <- m_stats * ef
  
  # now get the local FDR of the stats
  # Here, we use locfdr as it handles the two tails separately
  # and seems to give better results than fdrtool.
  
  lfdr <-  locfdr(stats, nulltype = 1)$fdr
  
  # impose minimum shrinkage
  ef_shrunk <- ef * pmax(1 - min_shrink, 1 - lfdr)
  
  res <- tibble(Species = names(ef),
                ef = ef,
                ef_shrunk = ef_shrunk,
                stat = stats,
                lfdr = lfdr)
  
  return(res)
}


compute_ef_YoungOldHSCs <- function(d, g, min_shrink = 0.3){
  
  raw_scores <-  sda.ranking(d,  g, diagonal = FALSE, verbose = TRUE,
                             fdr = FALSE)
  
  # undo scaling by class frequencies in order to get effect sizes
  g_summaries <- sda:::pvt.groups(g)
  freqs <- freqs.shrink(g_summaries$samples)
  m <- sqrt((1 - freqs)/freqs/length(g))
  
  # compute effect sizes
  # notice that the "cat.Yes" and "cat.No" have to be changed based on the factor
  ef <- raw_scores[, "cat.Young_HSCs"]*m["Young_HSCs"] - raw_scores[, "cat.Old_HSCs"]*m["Old_HSCs"]
  
  # compute tests statistics, first, we need the scaling factor
  # for the cat score between the two classes
  n0 <- sum(g_summaries$idx[, "Young_HSCs"])
  n1 <- sum(g_summaries$idx[, "Old_HSCs"])
  
  m_stats <- 1/sqrt(1/n0 + 1/n1)
  stats <- m_stats * ef
  
  # now get the local FDR of the stats
  # Here, we use locfdr as it handles the two tails separately
  # and seems to give better results than fdrtool.
  
  lfdr <-  locfdr(stats, nulltype = 1)$fdr
  
  # impose minimum shrinkage
  ef_shrunk <- ef * pmax(1 - min_shrink, 1 - lfdr)
  
  res <- tibble(Species = names(ef),
                ef = ef,
                ef_shrunk = ef_shrunk,
                stat = stats,
                lfdr = lfdr)
  
  return(res)
}


compute_ef_Recipient <- function(d, g, min_shrink = 0.3){
  
  raw_scores <-  sda.ranking(d,  g, diagonal = FALSE, verbose = TRUE,
                             fdr = FALSE)
  
  # undo scaling by class frequencies in order to get effect sizes
  g_summaries <- sda:::pvt.groups(g)
  freqs <- freqs.shrink(g_summaries$samples)
  m <- sqrt((1 - freqs)/freqs/length(g))
  
  # compute effect sizes
  # notice that the "cat.Yes" and "cat.No" have to be changed based on the factor
  ef <- raw_scores[, "cat.Young"]*m["Young"] - raw_scores[, "cat.Old"]*m["Old"]
  
  # compute tests statistics, first, we need the scaling factor
  # for the cat score between the two classes
  n0 <- sum(g_summaries$idx[, "Young"])
  n1 <- sum(g_summaries$idx[, "Old"])
  
  m_stats <- 1/sqrt(1/n0 + 1/n1)
  stats <- m_stats * ef
  
  # now get the local FDR of the stats
  # Here, we use locfdr as it handles the two tails separately
  # and seems to give better results than fdrtool.
  
  lfdr <-  locfdr(stats, nulltype = 1)$fdr
  
  # impose minimum shrinkage
  ef_shrunk <- ef * pmax(1 - min_shrink, 1 - lfdr)
  
  res <- tibble(Species = names(ef),
                ef = ef,
                ef_shrunk = ef_shrunk,
                stat = stats,
                lfdr = lfdr)
  
  return(res)
}

glog2 <- function(x) ((asinh(x) - log(2))/log(2))

# function to compute a robust mean
robustMean <- function(x){
  huberM(x)$mu
}


# Create SampleListFormatted to use for renaming samples.
Metadata<-read.csv("GeigerSampleKeyRevised20220910_Corrected.csv", header = TRUE, stringsAsFactors = FALSE, sep = ",")
names(Metadata)[1]<-"SampleID"
Metadata<-Metadata[,1:7]

Metadata$Background<-factor(Metadata$Background, levels =  c("C57BL", "RAG1"))
Metadata$Transplant<-factor(Metadata$Transplant, levels =  c("None", "Young_HSCs", "Old_HSCs", "Casin_treated_Old_HSCs", "RAG1"))
Metadata$Recipient<-factor(Metadata$Recipient, levels = c("Young", "Old"))
Metadata$ImmuneSystem<-factor(Metadata$ImmuneSystem, levels = c("Young", "Rejuvenated", "Old"))
Metadata$Cage<-factor(Metadata$Cage)
Metadata$Experiment<-factor(Metadata$Experiment)

Metadata$Groups<-paste(Metadata$Recipient, Metadata$Background, Metadata$Transplant, Metadata$ImmuneSystem, sep = "-")
Metadata$CageGroups<-paste(Metadata$Groups, Metadata$Cage, sep = "-")

Metadata$Groups<-gsub("Young-RAG1-Casin_treated_Old_HSCs-Rejuvenated", "RAG1-Rejuvenated_HSCs", Metadata$Groups)
Metadata$Groups<-gsub("Young-RAG1-Young_HSCs-Young", "RAG1-Young_HSCs", Metadata$Groups)
Metadata$Groups<-gsub("Young-RAG1-Old_HSCs-Old", "RAG1-Old_HSCs", Metadata$Groups)
Metadata$Groups<-gsub("Young-RAG1-None-Young", "Young_RAG1-Untransplanted", Metadata$Groups)
Metadata$Groups<-gsub("Old-C57BL-None-Old", "Old_C57BL-Untransplanted", Metadata$Groups)
Metadata$Groups<-gsub("Young-C57BL-None-Young", "Young_C57BL-Untransplanted", Metadata$Groups)
Metadata$Groups<-factor(Metadata$Groups, levels = c("Young_C57BL-Untransplanted", "Old_C57BL-Untransplanted", 
            "Young_RAG1-Untransplanted", "RAG1-Young_HSCs", "RAG1-Old_HSCs", "RAG1-Rejuvenated_HSCs", "Young-RAG1-RAG1-Young"))
Metadata$GroupsExp<-paste(Metadata$Groups, Metadata$Experiment, sep = "-")


GeigerSamples<-Metadata$SampleID

# Get the files
setwd("C:/Users/dbhas/OneDrive/Documents/Alignments/KrakenAlignments/Kraken2")
setwd("C:/Users/HASI9S/OneDrive/Documents/Alignments/KrakenAlignments/Kraken2")
setwd("~/Documents/Alignments/KrakenAlignments/Kraken2")

AllKrakenFiles<-list.files()
SpeciesFileList<-grep("_species_abundance.txt", AllKrakenFiles)
SpeciesFiles<-AllKrakenFiles[SpeciesFileList]
FileList<-gsub("_species_abundance.txt", "", SpeciesFiles)

NewSpeciesFileList<-subset(FileList, FileList %in% GeigerSamples)

#setwd("~/Documents/Code/Metagenomics/HD")

# get the files
for(f in 1:length(NewSpeciesFileList)){
  fnr = NewSpeciesFileList[f]
  x =NewSpeciesFileList[f]
  
  # assign(fnr, read.csv(paste("/home/david/Databases/FecalABDR/NewKrakenAlignments/",x, "_Species_abundance.txt",sep=""),
  #                      sep="\t", header = TRUE, stringsAsFactors = FALSE))
  
  assign(fnr, read.csv(paste(x, "_species_abundance.txt",sep=""),
                       sep="\t", header = TRUE, stringsAsFactors = FALSE))
}


# NewSpeciesNR<-read.csv(paste("/home/david/Databases/FecalABDR/NewKrakenAlignments/",NewSpeciesFileList[1], "_Species_abundance.txt", sep=""), sep = "\t",header = TRUE)
NewSpeciesNR<-read.csv(paste(NewSpeciesFileList[1], "_species_abundance.txt", sep=""), sep = "\t",header = TRUE)

names(NewSpeciesNR)[1]<-"Species"
NewSpeciesNR$Species<-gsub(" ", ".", NewSpeciesNR$Species, fixed = TRUE)
NewSpeciesNR$Species<-gsub("..", ".", NewSpeciesNR$Species, fixed = TRUE)
NewSpeciesNR$Species<-gsub("_", ".", NewSpeciesNR$Species, fixed = TRUE)
NewSpeciesNR$Species<-gsub("-", ".", NewSpeciesNR$Species, fixed = TRUE)
NewSpeciesNR$Species<-gsub("/", ".", NewSpeciesNR$Species, fixed = TRUE)
NewSpeciesNR$Species<-gsub("X,", "", NewSpeciesNR$Species, fixed = TRUE)
NewSpeciesNR$Species<-gsub("[", "", NewSpeciesNR$Species, fixed = TRUE)
NewSpeciesNR$Species<-gsub("]", "", NewSpeciesNR$Species, fixed = TRUE)
NewSpeciesNR<-as.data.frame(NewSpeciesNR[,c(1,4)])
names(NewSpeciesNR)<-c("Species",paste(paste(NewSpeciesFileList[1])))
NewSpeciesNR<-subset(NewSpeciesNR, ! duplicated(NewSpeciesNR$Species))

for ( i in 2:length(NewSpeciesFileList)){
  x <- get( NewSpeciesFileList[i] )
  x<-as.data.frame(x[,c(1,4)])
  names( x ) <- c("Species", paste(NewSpeciesFileList[i]))
  x$Species<-gsub(" ", ".", x$Species, fixed = TRUE)
  x$Species<-gsub("..", ".", x$Species, fixed = TRUE)
  x$Species<-gsub("_", ".", x$Species, fixed = TRUE)
  x$Species<-gsub("-", ".", x$Species, fixed = TRUE)
  x$Species<-gsub("/", ".", x$Species, fixed = TRUE)
  x$Species<-gsub("X,", "", x$Species, fixed = TRUE)
  x$Species<-gsub("[", "", x$Species, fixed = TRUE)
  x$Species<-gsub("]", "", x$Species, fixed = TRUE)
  x<-subset(x, ! duplicated(x$Species))
  NewSpeciesNR<-merge(x, NewSpeciesNR, by = "Species", all = TRUE)
}

row.names(NewSpeciesNR)<-NewSpeciesNR$Species
NewSpeciesNR$Species<-NULL
NewSpeciesNR[is.na(NewSpeciesNR)]<-0

rm(list = FileList)

setwd("C:/Users/HASI9S/OneDrive/Documents/Code/Metagenomics/GeigerData/StemCellsAgeing")
setwd("C:/Users/dbhas/OneDrive/Documents/Code/Metagenomics/GeigerData/StemCellsAgeing")
setwd("~/Documents/Code/Metagenomics/GeigerData/StemCellsAgeing")

# Remove samples in less than 10% of samples and less than 0.005% of microbiome

NewSpeciesNR<-as.data.frame(t(NewSpeciesNR))
SaliniCol<-grep("Salinibacter", colnames(NewSpeciesNR))
# NewSpeciesNR<-NewSpeciesNR[, -SaliniCol]
# HumanRow<-grep("Homo.sapiens", colnames(NewSpeciesNR))
# NewSpeciesNR<-NewSpeciesNR[, -HumanRow]
# MouseRow<-grep("Mus.musculus", colnames(NewSpeciesNR))
# NewSpeciesNR<-NewSpeciesNR[, -MouseRow]


GeigerSpecies<-NewSpeciesNR

# HumanRow<-grep("Homo.sapiens", colnames(NewSpeciesNR))
# GeigerSpeciesNoHuman<-NewSpeciesNR[, -HumanRow]


#### Ok this is without human DNA
TenPercentCutoff<-floor(nrow(GeigerSpecies)/20)

NonZeroCounts<-list()
for (i in 1:ncol(GeigerSpecies)){
  NonZeroCounts[i]<-length(which(GeigerSpecies[,i] > 0))
  
}

TenPercentNotZero<-which(NonZeroCounts >= TenPercentCutoff)
GeigerSpeciesNR<-GeigerSpecies[,TenPercentNotZero]

# Ok let's not remove species that are missing in 10% of samples since the sample number are low
GeigerSpeciesNR<-as.data.frame(t(noise.removal(t(GeigerSpeciesNR), 0.001)))
#GeigerSpeciesNR<-as.data.frame(t(noise.removal(t(GeigerSpecies), 0.01)))
LowSamples<-which(rowSums(GeigerSpeciesNR) <= 750000)
# GeigerSpeciesNR<-GeigerSpeciesNR[-LowSamples,] 
minCount<-min(rowSums(GeigerSpeciesNR))
GeigerSpeciesNR<-data.frame(rrarefy(GeigerSpeciesNR, 750000))
GeigerSpeciesNR$SampleID<-row.names(GeigerSpeciesNR)

BackupGeigerSpecies<-GeigerSpeciesNR
#GeigerSpeciesNR<-BackupGeigerSpecies


GeigerSpeciesNR<-merge(Metadata, GeigerSpeciesNR, by="SampleID", all.y = TRUE)
row.names(GeigerSpeciesNR)<-GeigerSpeciesNR$SampleID


Species<-as.data.frame(t(GeigerSpeciesNR[,c(11:ncol(GeigerSpeciesNR))]))

# Calculate sample diversity using 'diversity' function in Vegan
H <- data.frame(diversity(t(Species)))
simpson <- data.frame(diversity(t(Species), "simpson"))
shannon<-data.frame(diversity(t(Species), "shannon"))
invsimp <- data.frame(diversity(t(Species), "inv"))
alpha <- data.frame(fisher.alpha(t(Species)))
## Species richness (S) and Pielou's evenness (J):
S <- data.frame(specnumber(t(Species)))
J <- data.frame(H/log(S))
Diversity<-cbind(simpson, shannon, invsimp,  S, J)
Diversity$SampleID<-row.names(GeigerSpeciesNR)
#Diversity<-Diversity[,c(1,3,5,7,8,9,10)]
names(Diversity)<-c("Simpson", "Shannon",  "InvSimpson",  "SpeciesNo", "Evenness", "SampleID")

pairs(cbind(simpson, shannon,  J, S), pch="+", col="blue")


Diversity<-merge(Metadata, Diversity, by = "SampleID", all.y = TRUE)

col=tableau.colors


Diversity_Groups<-ggplot(Diversity,
                         aes(x=Groups, y=Shannon, fill=Groups)) + geom_boxplot(lwd=1,aes(color=factor(Groups),fill = NA), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Groups))) + xlab(NULL) +  
  ylab("Shannon Diversity Index")  +  paramsAngled() + scale_color_tableau() 
ggsave(filename = "Diversity_NewGroups.pdf", plot = Diversity_Groups, width = 10,  height = 10, limitsize = FALSE)


Diversity_NewGroups_Experiment<-ggplot(Diversity,
                         aes(x=Groups, y=Shannon, fill=Groups)) + geom_boxplot(lwd=1,aes(color=factor(Groups),fill = NA), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Groups))) + xlab(NULL) +  facet_grid(facets = . ~ Experiment) +
  ylab("Shannon Diversity Index")  +  paramsAngled() + scale_color_tableau() 
ggsave(filename = "Diversity_NewGroups_Experiment.pdf", plot = Diversity_NewGroups_Experiment, width = 10,  height = 10, limitsize = FALSE)

pairwise.wilcox.test(Diversity$Shannon, Diversity$Groups) # Null-DSS versus WT-RAG1-Old_HSCs is significant 


# Multiresponse permutation procedures (MRPP) 
Groups.mrpp<-mrpp(GeigerSpeciesNR[,11:length(colnames(GeigerSpeciesNR))], GeigerSpeciesNR[,8], distance = "bray")
Groups.mrpp # p 0.001


GeigerSpeciesNR2<-as.data.frame(t(noise.removal(t(GeigerSpeciesNR[,11:ncol(GeigerSpeciesNR)]), 0.05)))

SpeciesLog<-log(GeigerSpeciesNR2)
SpeciesCutLog<-log2(GeigerSpeciesNR2)
SpeciesCutLog[SpeciesCutLog==-Inf]<-0

HellingerSpecies<- tran(SpeciesCutLog, method="hellinger")
Cluster <- vegdist(SpeciesCutLog, method = "euclidean", diag = FALSE, upper = FALSE)

row.clus <- hclust(Cluster, "ward.D2")
cluster2<-vegdist(t(SpeciesCutLog), method="bray")
col.culst<-hclust(cluster2, "ward.D2")

heatmap.2(as.matrix(SpeciesCutLog), Rowv = as.dendrogram(row.clus), Colv=as.dendrogram(col.culst),col=metaphlan.colors(100),  margins=c(15,10),trace = "none",
          xlab = NULL, main = NULL, lhei= c(5,20))



metadata <- GeigerSpeciesNR[, 1:11]
cts <- as.matrix(GeigerSpeciesNR[, -(1:11)])
rownames(cts) <- metadata$SampleID
grps<-metadata$GroupsExp

cts_l2 <- glog2(cts)

TransplantedInd<-GeigerSpeciesNR$SampleID[GeigerSpeciesNR$Transplant %in% c("Young_HSCs", "Old_HSCs")]

col<-rep(tableau.colors, 2)
col=c("deepskyblue1", "deepskyblue2", "deepskyblue3", "deeppink1", "deeppink2", "deeppink3")
newPCA<-PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point",
                      pointsize = 2,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Species Abundance Colored by Mouse Group",
                      col.ind = grps,
                      addEllipses = TRUE,
                      ellipse.alpha = 0.01,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Mouse Group",
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank",
                      select.ind = list(name = TransplantedInd)
)
TransplantedGroupPCA<- pcaPlot 

pdf("TransplantedGroupPCA.pdf")
print(TransplantedGroupPCA)
dev.off()

UntransplantedInd<-GeigerSpeciesNR$SampleID[GeigerSpeciesNR$Transplant %in% c("None")& !GeigerSpeciesNR$Background == "RAG1"]

metadata <- GeigerSpeciesNR[, 1:11]
cts <- as.matrix(GeigerSpeciesNR[, -(1:11)])
rownames(cts) <- metadata$SampleID
grps<-metadata$GroupsExp

cts_l2 <- glog2(cts)

col<-rep(tableau.colors, 2)
newPCA<-PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point",
                      pointsize = 2,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Species Abundance Colored by Mouse Group",
                      col.ind = grps,
                      addEllipses = TRUE,
                      ellipse.alpha = 0.01,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Mouse Group",
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank",
                      select.ind = list(name = UntransplantedInd)
)
UntransplantedGroupPCA<- pcaPlot 

pdf("UntransplantedGroupPCA.pdf")
print(UntransplantedGroupPCA)
dev.off()<- pcaPlot 

save.image(file = "GeigerData20221006")

pcaPlot<-fviz_pca_ind(newPCA, axes = c(2,3),
                      geom.ind = "point",
                      pointsize = 2,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Species Abundance Colored by Mouse Group",
                      col.ind = grps,
                      addEllipses = TRUE,
                      ellipse.alpha = 0.01,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Mouse Group",
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank"
)
RotatedGroupPCA<- pcaPlot 

pdf("RotatedGroupPCA.pdf")
print(RotatedGroupPCA)
dev.off()


pcaPlot<-fviz_pca_ind(newPCA, axes = c(3,1),
                      geom.ind = "point",
                      pointsize = 2,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Species Abundance Colored by Mouse Group",
                      col.ind = grps,
                      addEllipses = TRUE,
                      ellipse.alpha = 0.01,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Mouse Group",
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank"
)
RotatedGroupPCA2<- pcaPlot 

pdf("RotatedGroupPCA2.pdf")
print(RotatedGroupPCA2)
dev.off()


metadata <- GeigerSpeciesNR[, 1:10]
cts <- as.matrix(GeigerSpeciesNR[, -(1:10)])
rownames(cts) <- metadata$SampleID
grps<-metadata$Cage

CageGroups<-metadata$SampleID[! metadata$Cage == "SRS0001-00440"] 

col=rep(tableau.colors,2)
newPCA<-PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point",
                      pointsize = 2,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Species Abundance Colored by Mouse Group",
                      col.ind = grps,
                      addEllipses = TRUE,
                      ellipse.alpha = 0.01,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Mouse Group",
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank",
                      select.ind = list(name = CageGroups)
)
CagePCA<- pcaPlot 

pdf("CagePCA.pdf")
print(CagePCA)
dev.off()


metadata <- GeigerSpeciesNR[, 1:10]
cts <- as.matrix(GeigerSpeciesNR[, -(1:10)])
rownames(cts) <- metadata$SampleID
grps<-metadata$CageGroup

CageGroups2<-metadata$SampleID[! metadata$Cage %in% c("SRS0001-00440", "SRS0001-00427",  "SHG0155-00026", "SHG0155-00023")] 

col=rep(tableau.colors,2)
newPCA<-PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point",
                      pointsize = 2,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Species Abundance Colored by Mouse Group",
                      col.ind = grps,
                      addEllipses = TRUE,
                      ellipse.alpha = 0.01,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Mouse and Cage Group",
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank",
                      select.ind = list(name = CageGroups2)
)
CageGroupPCA<- pcaPlot 

NotRag<-metadata$SampleID[! metadata$Background == "RAG1"]

pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "text",
                      pointsize = 2,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Species Abundance Colored by Mouse Group",
                      col.ind = grps,
                      addEllipses = FALSE,
                      ellipse.alpha = 0.01,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Mouse and Cage Group",
                      legend.size = 11,
                      mean.point = FALSE,
                      palette = col,
                      axes.linetype = "blank",
                      select.ind = list(name = NotRag),
)
CageGroupPCA<- pcaPlot 

pdf("CageGroupLabeled3PCA.pdf")
print(CagePCA)
dev.off()



metadata <- GeigerSpeciesNR[, 1:10]
cts <- as.matrix(GeigerSpeciesNR[, -(1:10)])
rownames(cts) <- metadata$SampleID
grps<-metadata$CageGroups


NotTransplanted<-metadata$SampleID[metadata$Transplant == "None"]
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point",
                      pointsize = 2,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Species Abundance Colored by Mouse Group",
                      col.ind = grps,
                      addEllipses = TRUE,
                      ellipse.alpha = 0.01,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Non-transplanted Group",
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank",
                      select.ind = list(name = NotTransplanted),
)
NoTransplantCagesPCA<- pcaPlot 

pdf("NoTransplantCagesPCA.pdf")
print(NoTransplantCagesPCA)
dev.off()


metadata <- GeigerSpeciesNR[, 1:10]
cts <- as.matrix(GeigerSpeciesNR[, -(1:10)])
rownames(cts) <- metadata$SampleID
grps<-metadata$Transplant


Transplanted<-metadata$SampleID[! metadata$Transplant == "None"]
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "text",
                      pointsize = 2,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Species Abundance Colored by Mouse Group",
                      col.ind = grps,
                      addEllipses = FALSE,
                      ellipse.alpha = 0.01,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Non-transplanted Group",
                      legend.size = 11,
                      mean.point = FALSE,
                      palette = col,
                      axes.linetype = "blank",
                      select.ind = list(name = Transplanted),
)
TransplantCagesPCA<- pcaPlot 

pdf("TransplantPCA.pdf")
print(TransplantCagesPCA)
dev.off()

save.image(file="GeigerData20220927")




geom<-"point"
fviz_pca_biplot(newPCA, axes = c(1, 2), geom = c("point"),
                geom.ind = geom, geom.var = c("text"), col.ind = "black",
                fill.ind = "white", col.var = "steelblue", fill.var = "white",
                gradient.cols = NULL, label = "all", invisible = "none", select.var =list(cos2 = 50),
                repel = FALSE, habillage = "none", palette = NULL, 
                addEllipses = FALSE, title = "PCA - Biplot")



#### alright now look at difference in species between young and old in non-transplanted mice  

YoungOldSamples<-subset(GeigerSpeciesNR, GeigerSpeciesNR$Transplant == "None" & GeigerSpeciesNR$Background == "C57BL")
YoungOldSamples<-YoungOldSamples[,c(1,2,11:ncol(YoungOldSamples))]

GeigerSpeciesNR3<-as.data.frame(t(noise.removal(t(YoungOldSamples[,11:ncol(YoungOldSamples)]), 0.035)))

SpeciesLog<-log(GeigerSpeciesNR3)
SpeciesCutLog<-log2(GeigerSpeciesNR3)
SpeciesCutLog[SpeciesCutLog==-Inf]<-0

HellingerSpecies<- tran(SpeciesCutLog, method="hellinger")
Cluster <- vegdist(SpeciesCutLog, method = "euclidean", diag = FALSE, upper = FALSE)

row.clus <- hclust(Cluster, "ward.D2")
cluster2<-vegdist(t(SpeciesCutLog), method="bray")
col.culst<-hclust(cluster2, "ward.D2")

heatmap.2(as.matrix(SpeciesCutLog), Rowv = as.dendrogram(row.clus), Colv=as.dendrogram(col.culst),col=metaphlan.colors(100),  margins=c(15,10),trace = "none",
          xlab = NULL, main = NULL, lhei= c(5,20))




Wilcox<-pairwise.wilcox.test(YoungOldSamples[,3], YoungOldSamples[,2])
WilcoxTable<-as.data.frame(Wilcox$p.value)
names(WilcoxTable)<-names(YoungOldSamples)[3]
row.names(WilcoxTable)<-"Wilcox"


for (i in 4:length(colnames(YoungOldSamples))){
  x<-pairwise.wilcox.test(YoungOldSamples[,i], YoungOldSamples[,2])
  y<-as.data.frame(x$p.value)
  names(y)<-names(YoungOldSamples)[i]
  WilcoxTable<-cbind(WilcoxTable, y)
}

WilcoxTable<-as.data.frame(t(WilcoxTable))
names(WilcoxTable)<-c("Unadjusted_p")
WilcoxTable$FDR<-p.adjust(WilcoxTable$Unadjusted_p, method = "fdr")
WilcoxTable$Species<-row.names(WilcoxTable)
WilcoxTable<-WilcoxTable[order(WilcoxTable$Unadjusted_p),]


SpeciesCounts<-as.data.frame(colSums(YoungOldSamples[3:length(colnames(YoungOldSamples))]))
names(SpeciesCounts)<-"SpeciesTotal"
SpeciesCountsBigTable<-as.data.frame(colSums(YoungOldSamples[3:length(colnames(YoungOldSamples))]))
names(SpeciesCountsBigTable)<-"SpeciesTotal"
TotalCountsAll<-sum(SpeciesCountsBigTable$SpeciesTotal)
SpeciesCounts$Species<-row.names(SpeciesCounts)

for (i in 1:length(rownames(SpeciesCounts))){
  SpeciesCounts$Fraction[i]<-SpeciesCounts[i,1] / TotalCountsAll
}

SpeciesCountMeans<-aggregate(.~YoungOldSamples$Recipient, mean, data=YoungOldSamples[3:ncol(YoungOldSamples)])
row.names(SpeciesCountMeans)<-SpeciesCountMeans$`YoungOldSamples$Recipient`
SpeciesCountMeans<-as.data.frame(t(SpeciesCountMeans[,2:ncol(SpeciesCountMeans)]))
SpeciesCountMeans$Species<-row.names(SpeciesCountMeans)
OverallMeans<-as.data.frame(colMeans(YoungOldSamples[3:ncol(YoungOldSamples)]))
OverallMeans$Species<-row.names(OverallMeans)
names(OverallMeans)<-c("OverallMean", "Species")
SpeciesCountMeans<-merge(SpeciesCountMeans, OverallMeans, by = "Species", all.x = TRUE)
SpeciesCountMeans$Young<-as.numeric(SpeciesCountMeans$Young)
SpeciesCountMeans$Old<-as.numeric(SpeciesCountMeans$Old)
SpeciesCountMeans$FoldChange<-(SpeciesCountMeans$Young-SpeciesCountMeans$Old) / SpeciesCountMeans$OverallMean


SpeciesCountMeans$Log2<-log2(SpeciesCountMeans$Young/SpeciesCountMeans$OverallMean)
SpeciesCountTable<-merge(SpeciesCountMeans, SpeciesCounts, by = "Species", all.x = TRUE)
SpeciesCountTable<-merge(SpeciesCountTable, WilcoxTable, by = "Species", all.x = TRUE)
SpeciesCountTable<-SpeciesCountTable[order(SpeciesCountTable$Unadjusted_p),]
SpeciesCountTable<-SpeciesCountTable[,c(1,8,4, 2,3,5,9,10)]
names(SpeciesCountTable)<-c("Species", "Abundance", "OverallMean", "Young.Mean", "Old.Mean", "Fold.Difference", "p_Unadjusted", "FDR")
SpeciesCountTable$Abundance<-SpeciesCountTable$Abundance * 100
SpeciesCountTable$Fold.Difference[is.na(SpeciesCountTable$Fold.Difference)]<-0
SpeciesCountTable<-SpeciesCountTable[order(SpeciesCountTable$Fold.Difference, decreasing = TRUE),]

SigSpeciesCountTable<-subset(SpeciesCountTable, SpeciesCountTable$p_Unadjusted < 0.05)
SigSpeciesCountTable$HigherIn<-ifelse(SigSpeciesCountTable$Fold.Difference > 0, "Young", "Old")
SigSpeciesCountTable$HigherIn<-factor(SigSpeciesCountTable$HigherIn, levels = c("Young", "Old"))
write.csv(SigSpeciesCountTable, file = "SignificantSpeciesUntransplantedYoungOld20221004.csv")
SigYoungOldNonTransplantedSpecies<-SigSpeciesCountTable$Species


metadata <- YoungOldSamples[, 1:2]
cts <- as.matrix(YoungOldSamples[, -(1:2)])
rownames(cts) <- metadata$SampleID
grps<-metadata$Recipient

cts_l2 <- glog2(cts)

res <- compute_ef_Recipient(d = cts_l2, g = YoungOldSamples$Recipient, min_shrink = 0.3)
names(res)[1]<-"Species"

res <- arrange(left_join(res, SigSpeciesCountTable), desc(abs(ef_shrunk))) %>%
  mutate(Species = as_factor(Species))

res<-subset(res, ! is.na(res$FDR))

res$HigherIn<-ifelse(res$ef_shrunk > 0, "Young", "Old")
res$HigherIn<-factor(res$HigherIn, levels = c("Young", "Old"))

resRecipient<-res

write.table(SigSpeciesCountTable, file = "SignificantRecipientSpecies.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)
# write.table(SpeciesCountTable, file = "GroupGeigerSpeciesNRK2.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)

col=tableau.colors[1:2]

Recipient_Effect<- ggplot(filter(resRecipient, resRecipient$ef_shrunk > 1.25 | resRecipient$ef_shrunk < -1.25),
                      aes(x = reorder(Species, ef_shrunk), y = ef_shrunk, fill = factor(HigherIn))) +
  geom_bar(colour="black", stat="identity") +
  #guides(fill=FALSE) + scale_fill_hue(l=40) +
  guides(fill=guide_legend(override.aes=list(colour=col))) +
  coord_flip() +
  labs( x = NULL) +
  ylab(expression(atop("Effect Size"))) +
  guides(fill=guide_legend(title="More abundant in:")) +
  #geom_hline(yintercept = c(-0.2, 0.2), color = "firebrick", lty = 2) +
  theme(axis.text.y = element_text(size= 12, color="black")) +
  theme(legend.title = element_text(size=12)) +
  theme(legend.text = element_text(size = 11)) 

col<-tableau.colors
Recipient_Effect<- Recipient_Effect + scale_fill_manual(values = col)

ggsave(filename = "RecipientYoungvsOld_Effect20221004.pdf", plot = Recipient_Effect, width = 12,
       height = 14, limitsize = FALSE)


DummySpecies<-GeigerSpeciesNR

col<-tableau.colors
Akkermansia.muciniphila<-ggplot(DummySpecies,
                              aes(x=Groups, y=Akkermansia.muciniphila, fill=NULL)) + geom_boxplot(lwd=1,aes(color=factor(Groups),fill = NULL), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
geom_point(size=4,aes(color = factor(Groups))) + xlab(NULL) + theme_bw() + facet_grid(facets = . ~ Experiment) +
  ylab("Akkermansia.muciniphila\n")  +  paramsAngled() + scale_color_tableau() + scale_y_log10()
ggsave(filename = "Akkermansia.muciniphila.pdf", plot = Akkermansia.muciniphila, width = 10,  height = 8, limitsize = FALSE)




### Ok now look at difference between young C57BL and RAG

GenotypeSamples<-subset(GeigerSpeciesNR, GeigerSpeciesNR$Background %in% c("C57BL", "RAG1") & GeigerSpeciesNR$Transplant == "None")
GenotypeSamples<-GenotypeSamples[,c(1,3,11:ncol(GenotypeSamples))]


Wilcox<-pairwise.wilcox.test(GenotypeSamples[,3], GenotypeSamples[,2])
WilcoxTable<-as.data.frame(Wilcox$p.value)
names(WilcoxTable)<-names(GenotypeSamples)[3]
row.names(WilcoxTable)<-"Wilcox"


for (i in 4:length(colnames(GenotypeSamples))){
  x<-pairwise.wilcox.test(GenotypeSamples[,i], GenotypeSamples[,2])
  y<-as.data.frame(x$p.value)
  names(y)<-names(GenotypeSamples)[i]
  WilcoxTable<-cbind(WilcoxTable, y)
}

WilcoxTable<-as.data.frame(t(WilcoxTable))
names(WilcoxTable)<-c("Unadjusted_p")
WilcoxTable$FDR<-p.adjust(WilcoxTable$Unadjusted_p, method = "fdr")
WilcoxTable$Species<-row.names(WilcoxTable)
WilcoxTable<-WilcoxTable[order(WilcoxTable$Unadjusted_p),]


SpeciesCounts<-as.data.frame(colSums(GenotypeSamples[3:length(colnames(GenotypeSamples))]))
names(SpeciesCounts)<-"SpeciesTotal"
SpeciesCountsBigTable<-as.data.frame(colSums(GenotypeSamples[3:length(colnames(GenotypeSamples))]))
names(SpeciesCountsBigTable)<-"SpeciesTotal"
TotalCountsAll<-sum(SpeciesCountsBigTable$SpeciesTotal)
SpeciesCounts$Species<-row.names(SpeciesCounts)

for (i in 1:length(rownames(SpeciesCounts))){
  SpeciesCounts$Fraction[i]<-SpeciesCounts[i,1] / TotalCountsAll
}

SpeciesCountMeans<-aggregate(.~GenotypeSamples$Background, mean, data=GenotypeSamples[3:ncol(GenotypeSamples)])
row.names(SpeciesCountMeans)<-SpeciesCountMeans$`GenotypeSamples$Background`
SpeciesCountMeans<-as.data.frame(t(SpeciesCountMeans[,2:ncol(SpeciesCountMeans)]))
SpeciesCountMeans$Species<-row.names(SpeciesCountMeans)
OverallMeans<-as.data.frame(colMeans(GenotypeSamples[3:ncol(GenotypeSamples)]))
OverallMeans$Species<-row.names(OverallMeans)
names(OverallMeans)<-c("OverallMean", "Species")
SpeciesCountMeans<-merge(SpeciesCountMeans, OverallMeans, by = "Species", all.x = TRUE)
SpeciesCountMeans$C57BL<-as.numeric(SpeciesCountMeans$C57BL)
SpeciesCountMeans$RAG1<-as.numeric(SpeciesCountMeans$RAG1)
SpeciesCountMeans$FoldChange<-(SpeciesCountMeans$C57BL-SpeciesCountMeans$RAG1) / SpeciesCountMeans$OverallMean
SpeciesCountMeans$Log2<-log2(SpeciesCountMeans$C57BL/SpeciesCountMeans$OverallMean)
SpeciesCountTable<-merge(SpeciesCountMeans, SpeciesCounts, by = "Species", all.x = TRUE)
SpeciesCountTable<-merge(SpeciesCountTable, WilcoxTable, by = "Species", all.x = TRUE)
SpeciesCountTable<-SpeciesCountTable[order(SpeciesCountTable$Unadjusted_p),]
SpeciesCountTable<-SpeciesCountTable[,c(1,8,4, 2,3,5,9,10)]
names(SpeciesCountTable)<-c("Species", "Abundance", "OverallMean", "C57BL.Mean", "RAG1.Mean", "Fold.Difference", "p_Unadjusted", "FDR")
SpeciesCountTable$Abundance<-SpeciesCountTable$Abundance * 100
SpeciesCountTable$Fold.Difference[is.na(SpeciesCountTable$Fold.Difference)]<-0
SpeciesCountTable<-SpeciesCountTable[order(SpeciesCountTable$Fold.Difference, decreasing = TRUE),]

SigSpeciesCountTable<-subset(SpeciesCountTable, SpeciesCountTable$p_Unadjusted < 0.05)
SigSpeciesCountTable$HigherIn<-ifelse(SigSpeciesCountTable$Fold.Difference > 0, "C57BL", "RAG1")
SigSpeciesCountTable$HigherIn<-factor(SigSpeciesCountTable$HigherIn, levels = c("C57BL", "RAG1"))
write.csv(SigSpeciesCountTable, file = "SignificantSpeciesC57BLvsRAG20221004.csv")
SigRAGC57BLSpecies<-SigSpeciesCountTable$Species


metadata <- GenotypeSamples[, 1:2]
cts <- as.matrix(GenotypeSamples[, -(1:2)])
rownames(cts) <- metadata$SampleID
grps<-metadata$Background

cts_l2 <- glog2(cts)

res <- compute_ef_Background(d = cts_l2, g = GenotypeSamples$Background, min_shrink = 0.3)
names(res)[1]<-"Species"

res <- arrange(left_join(res, SigSpeciesCountTable), desc(abs(ef_shrunk))) %>%
  mutate(Species = as_factor(Species))

res<-subset(res, ! is.na(res$FDR))

res$HigherIn<-ifelse(res$ef_shrunk > 0, "c57BL", "RAG1")
res$HigherIn<-factor(res$HigherIn, levels = c("c57BL", "RAG1"))

resBackground<-res

# write.table(SigSpeciesCountTable, file = "SignificantBackgroundSpecies.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)
# # write.table(SpeciesCountTable, file = "GroupGeigerSpeciesNRK2.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)

col=tableau.colors[1:2]

Background_Effect<- ggplot(filter(resBackground, resBackground$ef_shrunk > 10 | resBackground$ef_shrunk < -5),
                         aes(x = reorder(Species, ef_shrunk), y = ef_shrunk, fill = factor(HigherIn))) +
  geom_bar(colour="black", stat="identity") +
  #guides(fill=FALSE) + scale_fill_hue(l=40) +
  guides(fill=guide_legend(override.aes=list(colour=col))) +
  coord_flip() +
  labs( x = NULL) +
  ylab(expression(atop("Effect Size"))) +
  guides(fill=guide_legend(title="More abundant in:")) +
  #geom_hline(yintercept = c(-0.2, 0.2), color = "firebrick", lty = 2) +
  theme(axis.text.y = element_text(size= 12, color="black")) +
  theme(legend.title = element_text(size=12)) +
  theme(legend.text = element_text(size = 11)) 


Background_Effect<- Background_Effect + scale_fill_manual(values = col)

ggsave(filename = "RAGvsC57BL_Species_Effect20221004.pdf", plot = Background_Effect, width = 12,
       height = 16, limitsize = FALSE)


col<-tableau.colors
Alistipes.shahii<-ggplot(DummySpecies,
                                aes(x=Groups, y=Alistipes.shahii, fill=NULL)) + geom_boxplot(lwd=1,aes(color=factor(Groups),fill = NULL), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Groups))) + xlab(NULL) + theme_bw() + facet_grid(facets = . ~ Experiment) +
  ylab("Alistipes.shahii\n")  +  paramsAngled() + scale_color_tableau() 
ggsave(filename = "Alistipes.shahii.pdf", plot = Alistipes.shahii, width = 10,  height = 8, limitsize = FALSE)

Alistipes.obesi<-ggplot(DummySpecies,
                         aes(x=Groups, y=Alistipes.obesi, fill=NULL)) + geom_boxplot(lwd=1,aes(color=factor(Groups),fill = NULL), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Groups))) + xlab(NULL) + theme_bw() + facet_grid(facets = . ~ Experiment) +
  ylab("Alistipes.obesi\n")  +  paramsAngled() + scale_color_tableau() 
ggsave(filename = "Alistipes.obesi.pdf", plot = Alistipes.obesi, width = 10,  height = 8, limitsize = FALSE)


col<-tableau.colors
Lactobacillus.murinus<-ggplot(DummySpecies,
                         aes(x=Groups, y=Lactobacillus.murinus, fill=NULL)) + geom_boxplot(lwd=1,aes(color=factor(Groups),fill = NULL), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Groups))) + xlab(NULL) + theme_bw() + facet_grid(facets = . ~ Experiment) +
  ylab("Lactobacillus.murinus\n")  +  paramsAngled() + scale_color_tableau() 
ggsave(filename = "Lactobacillus.murinus.pdf", plot = Lactobacillus.murinus, width = 10,  height = 8, limitsize = FALSE)

save.image(file = "GeigerData20220927")


### alright look at Young_HSCs versus Old_HSCs ####

YoungOldHSC<-subset(GeigerSpeciesNR, GeigerSpeciesNR$Transplant %in% c("Young_HSCs", "Old_HSCs"))
YoungOldHSC<-YoungOldHSC[, c(1,4,11:ncol(YoungOldHSC))]

Wilcox<-pairwise.wilcox.test(YoungOldHSC[,3], YoungOldHSC[,2])
WilcoxTable<-as.data.frame(Wilcox$p.value)
names(WilcoxTable)<-names(YoungOldHSC)[3]
row.names(WilcoxTable)<-"Wilcox"


for (i in 4:length(colnames(YoungOldHSC))){
  x<-pairwise.wilcox.test(YoungOldHSC[,i], YoungOldHSC[,2])
  y<-as.data.frame(x$p.value)
  names(y)<-names(YoungOldHSC)[i]
  WilcoxTable<-cbind(WilcoxTable, y)
}

WilcoxTable<-as.data.frame(t(WilcoxTable))
names(WilcoxTable)<-c("Unadjusted_p")
WilcoxTable$FDR<-p.adjust(WilcoxTable$Unadjusted_p, method = "fdr")
WilcoxTable$Species<-row.names(WilcoxTable)
WilcoxTable<-WilcoxTable[order(WilcoxTable$Unadjusted_p),]


SpeciesCounts<-as.data.frame(colSums(YoungOldHSC[3:length(colnames(YoungOldHSC))]))
names(SpeciesCounts)<-"SpeciesTotal"
SpeciesCountsBigTable<-as.data.frame(colSums(YoungOldHSC[3:length(colnames(YoungOldHSC))]))
names(SpeciesCountsBigTable)<-"SpeciesTotal"
TotalCountsAll<-sum(SpeciesCountsBigTable$SpeciesTotal)
SpeciesCounts$Species<-row.names(SpeciesCounts)

for (i in 1:length(rownames(SpeciesCounts))){
  SpeciesCounts$Fraction[i]<-SpeciesCounts[i,1] / TotalCountsAll
}

SpeciesCountMeans<-aggregate(.~YoungOldHSC$Transplant, mean, data=YoungOldHSC[3:ncol(YoungOldHSC)])
row.names(SpeciesCountMeans)<-SpeciesCountMeans$`YoungOldHSC$Transplant`
SpeciesCountMeans<-as.data.frame(t(SpeciesCountMeans[,2:ncol(SpeciesCountMeans)]))
SpeciesCountMeans$Species<-row.names(SpeciesCountMeans)
OverallMeans<-as.data.frame(colMeans(YoungOldHSC[3:ncol(YoungOldHSC)]))
OverallMeans$Species<-row.names(OverallMeans)
names(OverallMeans)<-c("OverallMean", "Species")
SpeciesCountMeans<-merge(SpeciesCountMeans, OverallMeans, by = "Species", all.x = TRUE)
SpeciesCountMeans$Young_HSCs<-as.numeric(SpeciesCountMeans$Young_HSCs)
SpeciesCountMeans$Old_HSCs<-as.numeric(SpeciesCountMeans$Old_HSCs)
SpeciesCountMeans$FoldChange<-(SpeciesCountMeans$Young_HSCs-SpeciesCountMeans$Old_HSCs) / SpeciesCountMeans$OverallMean
SpeciesCountMeans$Log2<-log2(SpeciesCountMeans$Young_HSCs/SpeciesCountMeans$OverallMean)
SpeciesCountTable<-merge(SpeciesCountMeans, SpeciesCounts, by = "Species", all.x = TRUE)
SpeciesCountTable<-merge(SpeciesCountTable, WilcoxTable, by = "Species", all.x = TRUE)
SpeciesCountTable<-SpeciesCountTable[order(SpeciesCountTable$Unadjusted_p),]
SpeciesCountTable<-SpeciesCountTable[,c(1,8,4, 2,3,5,9,10)]
names(SpeciesCountTable)<-c("Species", "Abundance", "OverallMean", "Young_HSCs.Mean", "Old_HSCs.Mean", "Fold.Difference", "p_Unadjusted", "FDR")
SpeciesCountTable$Abundance<-SpeciesCountTable$Abundance * 100
SpeciesCountTable$Fold.Difference[is.na(SpeciesCountTable$Fold.Difference)]<-0
SpeciesCountTable<-SpeciesCountTable[order(SpeciesCountTable$Fold.Difference, decreasing = TRUE),]

SigSpeciesCountTable<-subset(SpeciesCountTable, SpeciesCountTable$p_Unadjusted < 0.05)
SigSpeciesCountTable$HigherIn<-ifelse(SigSpeciesCountTable$Fold.Difference > 0, "Young_HSCs", "Old_HSCs")
SigSpeciesCountTable$HigherIn<-factor(SigSpeciesCountTable$HigherIn, levels = c("Young_HSCs", "Old_HSCs"))
write.csv(SigSpeciesCountTable, file = "SignificantSpeciesYoungvsOldHSC20221004.csv")
SigYoungOldHSCSpecies<-SigSpeciesCountTable$Species

intersect(SigYoungOldNonTransplantedSpecies, SigYoungOldHSCSpecies)


metadata <- YoungOldHSC[, 1:2]
cts <- as.matrix(YoungOldHSC[, -(1:2)])
rownames(cts) <- metadata$SampleID
grps<-metadata$Transplant

cts_l2 <- glog2(cts)

res <- compute_ef_YoungOldHSCs(d = cts_l2, g = YoungOldHSC$Transplant, min_shrink = 0.3)
names(res)[1]<-"Species"

res <- arrange(left_join(res, SigSpeciesCountTable), desc(abs(ef_shrunk))) %>%
  mutate(Species = as_factor(Species))

res<-subset(res, ! is.na(res$FDR))

res$HigherIn<-ifelse(res$ef_shrunk > 0, "Young_HSCs", "Old_HSCs")
res$HigherIn<-factor(res$HigherIn, levels = c("Young_HSCs", "Old_HSCs"))

resYoungOldHSCs<-res

# write.table(SigSpeciesCountTable, file = "SignificantBackgroundSpecies.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)
# write.table(SpeciesCountTable, file = "GroupGeigerSpeciesNRK2.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)

col=tableau.colors[1:2]

YoungOldHSCEffect<- ggplot(filter(resYoungOldHSCs, resYoungOldHSCs$ef_shrunk > 0.3 | resYoungOldHSCs$ef_shrunk < -0.3),
                           aes(x = reorder(Species, ef_shrunk), y = ef_shrunk, fill = factor(HigherIn))) +
  geom_bar(colour="black", stat="identity") +
  #guides(fill=FALSE) + scale_fill_hue(l=40) +
  guides(fill=guide_legend(override.aes=list(colour=col))) +
  coord_flip() +
  labs( x = NULL) +
  ylab(expression(atop("Effect Size"))) +
  guides(fill=guide_legend(title="More abundant in:")) +
  #geom_hline(yintercept = c(-0.2, 0.2), color = "firebrick", lty = 2) +
  theme(axis.text.y = element_text(size= 12, color="black")) +
  theme(legend.title = element_text(size=12)) +
  theme(legend.text = element_text(size = 11)) 


YoungOldHSCEffect<- YoungOldHSCEffect + scale_fill_manual(values = col)

ggsave(filename = "YoungvsOldHSC_Species_Effect20221004.pdf", plot = YoungOldHSCEffect, width = 12,
       height = 6, limitsize = FALSE)

OverlappingSpecies<-intersect(SigYoungOldNonTransplantedSpecies, SigYoungOldHSCSpecies)


## Ok going to look at Young vs Old HSC in the two experiements and see which species are moving in same direction


### alright look at Young_HSCs versus Old_HSCs ####

YoungOldHSCExp1<-subset(GeigerSpeciesNR, GeigerSpeciesNR$Transplant %in% c("Young_HSCs", "Old_HSCs") & GeigerSpeciesNR$Experiment == 1 )
YoungOldHSCExp1<-YoungOldHSCExp1[, c(1,4,11:ncol(YoungOldHSCExp1))]

Wilcox<-pairwise.wilcox.test(YoungOldHSCExp1[,3], YoungOldHSCExp1[,2])
WilcoxTable<-as.data.frame(Wilcox$p.value)
names(WilcoxTable)<-names(YoungOldHSCExp1)[3]
row.names(WilcoxTable)<-"Wilcox"


for (i in 4:length(colnames(YoungOldHSCExp1))){
  x<-pairwise.wilcox.test(YoungOldHSCExp1[,i], YoungOldHSCExp1[,2])
  y<-as.data.frame(x$p.value)
  names(y)<-names(YoungOldHSCExp1)[i]
  WilcoxTable<-cbind(WilcoxTable, y)
}

WilcoxTable<-as.data.frame(t(WilcoxTable))
names(WilcoxTable)<-c("Unadjusted_p")
WilcoxTable$FDR<-p.adjust(WilcoxTable$Unadjusted_p, method = "fdr")
WilcoxTable$Species<-row.names(WilcoxTable)
WilcoxTable<-WilcoxTable[order(WilcoxTable$Unadjusted_p),]


SpeciesCounts<-as.data.frame(colSums(YoungOldHSCExp1[3:length(colnames(YoungOldHSCExp1))]))
names(SpeciesCounts)<-"SpeciesTotal"
SpeciesCountsBigTable<-as.data.frame(colSums(YoungOldHSCExp1[3:length(colnames(YoungOldHSCExp1))]))
names(SpeciesCountsBigTable)<-"SpeciesTotal"
TotalCountsAll<-sum(SpeciesCountsBigTable$SpeciesTotal)
SpeciesCounts$Species<-row.names(SpeciesCounts)

for (i in 1:length(rownames(SpeciesCounts))){
  SpeciesCounts$Fraction[i]<-SpeciesCounts[i,1] / TotalCountsAll
}

SpeciesCountMeans<-aggregate(.~YoungOldHSCExp1$Transplant, mean, data=YoungOldHSCExp1[3:ncol(YoungOldHSCExp1)])
row.names(SpeciesCountMeans)<-SpeciesCountMeans$`YoungOldHSCExp1$Transplant`
SpeciesCountMeans<-as.data.frame(t(SpeciesCountMeans[,2:ncol(SpeciesCountMeans)]))
SpeciesCountMeans$Species<-row.names(SpeciesCountMeans)
OverallMeans<-as.data.frame(colMeans(YoungOldHSCExp1[3:ncol(YoungOldHSCExp1)]))
OverallMeans$Species<-row.names(OverallMeans)
names(OverallMeans)<-c("OverallMean", "Species")
SpeciesCountMeans<-merge(SpeciesCountMeans, OverallMeans, by = "Species", all.x = TRUE)
SpeciesCountMeans$Young_HSCs<-as.numeric(SpeciesCountMeans$Young_HSCs)
SpeciesCountMeans$Old_HSCs<-as.numeric(SpeciesCountMeans$Old_HSCs)
SpeciesCountMeans$FoldChange<-(SpeciesCountMeans$Young_HSCs-SpeciesCountMeans$Old_HSCs) / SpeciesCountMeans$OverallMean
SpeciesCountMeans$Log2<-log2(SpeciesCountMeans$Young_HSCs/SpeciesCountMeans$OverallMean)
SpeciesCountTable<-merge(SpeciesCountMeans, SpeciesCounts, by = "Species", all.x = TRUE)
SpeciesCountTable<-merge(SpeciesCountTable, WilcoxTable, by = "Species", all.x = TRUE)
SpeciesCountTable<-SpeciesCountTable[order(SpeciesCountTable$Unadjusted_p),]
SpeciesCountTable<-SpeciesCountTable[,c(1,8,4, 2,3,5,9,10)]
names(SpeciesCountTable)<-c("Species", "Abundance", "OverallMean", "Young_HSCs.Mean", "Old_HSCs.Mean", "Fold.Difference", "p_Unadjusted", "FDR")
SpeciesCountTable$Abundance<-SpeciesCountTable$Abundance * 100
SpeciesCountTable$Fold.Difference[is.na(SpeciesCountTable$Fold.Difference)]<-0
SpeciesCountTable<-SpeciesCountTable[order(SpeciesCountTable$Fold.Difference, decreasing = TRUE),]

SigSpeciesCountTable<-subset(SpeciesCountTable, SpeciesCountTable$p_Unadjusted < 0.05)
SigSpeciesCountTable$HigherIn<-ifelse(SigSpeciesCountTable$Fold.Difference > 0, "Young_HSCs", "Old_HSCs")
SigSpeciesCountTable$HigherIn<-factor(SigSpeciesCountTable$HigherIn, levels = c("Young_HSCs", "Old_HSCs"))
write.csv(SigSpeciesCountTable, file = "SignificantSpeciesYoungvsOldHSC.csv")
SigYoungOldHSCExp1Species<-SigSpeciesCountTable$Species


metadata <- YoungOldHSCExp1[, 1:2]
cts <- as.matrix(YoungOldHSCExp1[, -(1:2)])
rownames(cts) <- metadata$SampleID
grps<-metadata$Transplant

cts_l2 <- glog2(cts)

res <- compute_ef_YoungOldHSCs(d = cts_l2, g = YoungOldHSCExp1$Transplant, min_shrink = 0.3)
names(res)[1]<-"Species"

res <- arrange(left_join(res, SigSpeciesCountTable), desc(abs(ef_shrunk))) %>%
  mutate(Species = as_factor(Species))

res<-subset(res, ! is.na(res$FDR))

res$HigherIn<-ifelse(res$ef_shrunk > 0, "Young_HSCs", "Old_HSCs")
res$HigherIn<-factor(res$HigherIn, levels = c("Young_HSCs", "Old_HSCs"))

resYoungOldHSCExp1s<-res

# write.table(SigSpeciesCountTable, file = "SignificantBackgroundSpecies.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)
# write.table(SpeciesCountTable, file = "GroupGeigerSpeciesNRK2.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)

col=tableau.colors[1:2]

YoungOldHSCExp1Effect<- ggplot(filter(resYoungOldHSCExp1s, resYoungOldHSCExp1s$ef_shrunk > 1 | resYoungOldHSCExp1s$ef_shrunk < -1),
                           aes(x = reorder(Species, ef_shrunk), y = ef_shrunk, fill = factor(HigherIn))) +
  geom_bar(colour="black", stat="identity") +
  #guides(fill=FALSE) + scale_fill_hue(l=40) +
  guides(fill=guide_legend(override.aes=list(colour=col))) +
  coord_flip() +
  labs( x = NULL) +
  ylab(expression(atop("Effect Size"))) +
  guides(fill=guide_legend(title="More abundant in:")) +
  #geom_hline(yintercept = c(-0.2, 0.2), color = "firebrick", lty = 2) +
  theme(axis.text.y = element_text(size= 12, color="black")) +
  theme(legend.title = element_text(size=12)) +
  theme(legend.text = element_text(size = 11)) 


YoungOldHSCExp1Effect<- YoungOldHSCExp1Effect + scale_fill_manual(values = col)

ggsave(filename = "YoungvsOldHSCLExpt1_Species_Effect.pdf", plot = YoungOldHSCExp1Effect, width = 12,
       height = 16, limitsize = FALSE)

## Now expt 2


YoungOldHSCExp2<-subset(GeigerSpeciesNR, GeigerSpeciesNR$Transplant %in% c("Young_HSCs", "Old_HSCs") & GeigerSpeciesNR$Experiment == 2 )
YoungOldHSCExp2<-YoungOldHSCExp2[, c(1,4,11:ncol(YoungOldHSCExp2))]

Wilcox<-pairwise.wilcox.test(YoungOldHSCExp2[,3], YoungOldHSCExp2[,2])
WilcoxTable<-as.data.frame(Wilcox$p.value)
names(WilcoxTable)<-names(YoungOldHSCExp2)[3]
row.names(WilcoxTable)<-"Wilcox"


for (i in 4:length(colnames(YoungOldHSCExp2))){
  x<-pairwise.wilcox.test(YoungOldHSCExp2[,i], YoungOldHSCExp2[,2])
  y<-as.data.frame(x$p.value)
  names(y)<-names(YoungOldHSCExp2)[i]
  WilcoxTable<-cbind(WilcoxTable, y)
}

WilcoxTable<-as.data.frame(t(WilcoxTable))
names(WilcoxTable)<-c("Unadjusted_p")
WilcoxTable$FDR<-p.adjust(WilcoxTable$Unadjusted_p, method = "fdr")
WilcoxTable$Species<-row.names(WilcoxTable)
WilcoxTable<-WilcoxTable[order(WilcoxTable$Unadjusted_p),]


SpeciesCounts<-as.data.frame(colSums(YoungOldHSCExp2[3:length(colnames(YoungOldHSCExp2))]))
names(SpeciesCounts)<-"SpeciesTotal"
SpeciesCountsBigTable<-as.data.frame(colSums(YoungOldHSCExp2[3:length(colnames(YoungOldHSCExp2))]))
names(SpeciesCountsBigTable)<-"SpeciesTotal"
TotalCountsAll<-sum(SpeciesCountsBigTable$SpeciesTotal)
SpeciesCounts$Species<-row.names(SpeciesCounts)

for (i in 1:length(rownames(SpeciesCounts))){
  SpeciesCounts$Fraction[i]<-SpeciesCounts[i,1] / TotalCountsAll
}

SpeciesCountMeans<-aggregate(.~YoungOldHSCExp2$Transplant, mean, data=YoungOldHSCExp2[3:ncol(YoungOldHSCExp2)])
row.names(SpeciesCountMeans)<-SpeciesCountMeans$`YoungOldHSCExp2$Transplant`
SpeciesCountMeans<-as.data.frame(t(SpeciesCountMeans[,2:ncol(SpeciesCountMeans)]))
SpeciesCountMeans$Species<-row.names(SpeciesCountMeans)
OverallMeans<-as.data.frame(colMeans(YoungOldHSCExp2[3:ncol(YoungOldHSCExp2)]))
OverallMeans$Species<-row.names(OverallMeans)
names(OverallMeans)<-c("OverallMean", "Species")
SpeciesCountMeans<-merge(SpeciesCountMeans, OverallMeans, by = "Species", all.x = TRUE)
SpeciesCountMeans$Young_HSCs<-as.numeric(SpeciesCountMeans$Young_HSCs)
SpeciesCountMeans$Old_HSCs<-as.numeric(SpeciesCountMeans$Old_HSCs)
SpeciesCountMeans$FoldChange<-(SpeciesCountMeans$Young_HSCs-SpeciesCountMeans$Old_HSCs) / SpeciesCountMeans$OverallMean
SpeciesCountMeans$Log2<-log2(SpeciesCountMeans$Young_HSCs/SpeciesCountMeans$OverallMean)
SpeciesCountTable<-merge(SpeciesCountMeans, SpeciesCounts, by = "Species", all.x = TRUE)
SpeciesCountTable<-merge(SpeciesCountTable, WilcoxTable, by = "Species", all.x = TRUE)
SpeciesCountTable<-SpeciesCountTable[order(SpeciesCountTable$Unadjusted_p),]
SpeciesCountTable<-SpeciesCountTable[,c(1,8,4, 2,3,5,9,10)]
names(SpeciesCountTable)<-c("Species", "Abundance", "OverallMean", "Young_HSCs.Mean", "Old_HSCs.Mean", "Fold.Difference", "p_Unadjusted", "FDR")
SpeciesCountTable$Abundance<-SpeciesCountTable$Abundance * 100
SpeciesCountTable$Fold.Difference[is.na(SpeciesCountTable$Fold.Difference)]<-0
SpeciesCountTable<-SpeciesCountTable[order(SpeciesCountTable$Fold.Difference, decreasing = TRUE),]

SigSpeciesCountTable<-subset(SpeciesCountTable, SpeciesCountTable$p_Unadjusted < 0.05)
SigSpeciesCountTable$HigherIn<-ifelse(SigSpeciesCountTable$Fold.Difference > 0, "Young_HSCs", "Old_HSCs")
SigSpeciesCountTable$HigherIn<-factor(SigSpeciesCountTable$HigherIn, levels = c("Young_HSCs", "Old_HSCs"))
write.csv(SigSpeciesCountTable, file = "SignificantSpeciesYoungvsOldHSC.csv")
SigYoungOldHSCExp2Species<-SigSpeciesCountTable$Species


metadata <- YoungOldHSCExp2[, 1:2]
cts <- as.matrix(YoungOldHSCExp2[, -(1:2)])
rownames(cts) <- metadata$SampleID
grps<-metadata$Transplant

cts_l2 <- glog2(cts)

res <- compute_ef_YoungOldHSCs(d = cts_l2, g = YoungOldHSCExp2$Transplant, min_shrink = 0.3)
names(res)[1]<-"Species"

res <- arrange(left_join(res, SigSpeciesCountTable), desc(abs(ef_shrunk))) %>%
  mutate(Species = as_factor(Species))

res<-subset(res, ! is.na(res$FDR))

res$HigherIn<-ifelse(res$ef_shrunk > 0, "Young_HSCs", "Old_HSCs")
res$HigherIn<-factor(res$HigherIn, levels = c("Young_HSCs", "Old_HSCs"))

resYoungOldHSCExp2s<-res

# write.table(SigSpeciesCountTable, file = "SignificantBackgroundSpecies.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)
# write.table(SpeciesCountTable, file = "GroupGeigerSpeciesNRK2.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)

col=tableau.colors[1:2]

YoungOldHSCExp2Effect<- ggplot(filter(resYoungOldHSCExp2s, resYoungOldHSCExp2s$ef_shrunk > 1 | resYoungOldHSCExp2s$ef_shrunk < -1),
                               aes(x = reorder(Species, ef_shrunk), y = ef_shrunk, fill = factor(HigherIn))) +
  geom_bar(colour="black", stat="identity") +
  #guides(fill=FALSE) + scale_fill_hue(l=40) +
  guides(fill=guide_legend(override.aes=list(colour=col))) +
  coord_flip() +
  labs( x = NULL) +
  ylab(expression(atop("Effect Size"))) +
  guides(fill=guide_legend(title="More abundant in:")) +
  #geom_hline(yintercept = c(-0.2, 0.2), color = "firebrick", lty = 2) +
  theme(axis.text.y = element_text(size= 12, color="black")) +
  theme(legend.title = element_text(size=12)) +
  theme(legend.text = element_text(size = 11)) 


YoungOldHSCExp2Effect<- YoungOldHSCExp2Effect + scale_fill_manual(values = col)

ggsave(filename = "YoungvsOldHSCLExpt2_Species_Effect.pdf", plot = YoungOldHSCExp2Effect, width = 12,
       height = 16, limitsize = FALSE)




OverlappingYOHSCSpecies<-intersect(resYoungOldHSCExp1s$Species, resYoungOldHSCExp2s$Species)




# can make a Venn diagram
SpeciesVennDiagram<-draw.pairwise.venn(area1 = length(SigYoungOldNonTransplantedSpecies), area2 = length(SigYoungOldHSCSpecies),
                                       cross.area = length(OverlappingSpecies), category = c("Young.Old.Mice",  "Young.Old.HSC"), lty = rep("blank", 2),
                                       fill = c("light blue", "pink"),  cat.pos =c(0, 0),
                                       fontfamily = "Helvetica", fontface = "plain", cat.default.pos = "outer", cat.fontfamily = "Helvetica")

ggsave(filename = "SpeciesVennDiagram.pdf", plot = SpeciesVennDiagram, width = 6,
       height = 6, limitsize = FALSE)

## How did I loop through graphs again...? ####
# this works!#

col=tableau.colors

for(i in 1:length(OverlappingSpecies)){
  
  Yaxis<-grep(paste(OverlappingSpecies[i]), names(DummySpecies))

Fig<-ggplot(DummySpecies,aes(x=Groups, y=DummySpecies[,Yaxis], fill=NULL)) + geom_boxplot(lwd=1,aes(color=factor(Groups),fill = NULL), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Groups))) + xlab(NULL) + theme_bw() + facet_grid(facets = . ~ Experiment) +
  ylab(paste(OverlappingSpecies[i]))  +  paramsAngled() + scale_color_tableau() + scale_y_log10()
ggsave(filename = paste(OverlappingSpecies[i], ".pdf", sep=""), plot = Fig, width = 10,  height = 8, limitsize = FALSE)

}

save.image(file = "GeigerData20220927")


## Now expt 3


YoungOldHSCExp3<-subset(GeigerSpeciesNR, GeigerSpeciesNR$Transplant %in% c("Young_HSCs", "Old_HSCs") & GeigerSpeciesNR$Experiment == 3 )
YoungOldHSCExp3<-YoungOldHSCExp3[, c(1,4,11:ncol(YoungOldHSCExp3))]

Wilcox<-pairwise.wilcox.test(YoungOldHSCExp3[,3], YoungOldHSCExp3[,2])
WilcoxTable<-as.data.frame(Wilcox$p.value)
names(WilcoxTable)<-names(YoungOldHSCExp3)[3]
row.names(WilcoxTable)<-"Wilcox"


for (i in 4:length(colnames(YoungOldHSCExp3))){
  x<-pairwise.wilcox.test(YoungOldHSCExp3[,i], YoungOldHSCExp3[,2])
  y<-as.data.frame(x$p.value)
  names(y)<-names(YoungOldHSCExp3)[i]
  WilcoxTable<-cbind(WilcoxTable, y)
}

WilcoxTable<-as.data.frame(t(WilcoxTable))
names(WilcoxTable)<-c("Unadjusted_p")
WilcoxTable$FDR<-p.adjust(WilcoxTable$Unadjusted_p, method = "fdr")
WilcoxTable$Species<-row.names(WilcoxTable)
WilcoxTable<-WilcoxTable[order(WilcoxTable$Unadjusted_p),]


SpeciesCounts<-as.data.frame(colSums(YoungOldHSCExp3[3:length(colnames(YoungOldHSCExp3))]))
names(SpeciesCounts)<-"SpeciesTotal"
SpeciesCountsBigTable<-as.data.frame(colSums(YoungOldHSCExp3[3:length(colnames(YoungOldHSCExp3))]))
names(SpeciesCountsBigTable)<-"SpeciesTotal"
TotalCountsAll<-sum(SpeciesCountsBigTable$SpeciesTotal)
SpeciesCounts$Species<-row.names(SpeciesCounts)

for (i in 1:length(rownames(SpeciesCounts))){
  SpeciesCounts$Fraction[i]<-SpeciesCounts[i,1] / TotalCountsAll
}

SpeciesCountMeans<-aggregate(.~YoungOldHSCExp3$Transplant, mean, data=YoungOldHSCExp3[3:ncol(YoungOldHSCExp3)])
row.names(SpeciesCountMeans)<-SpeciesCountMeans$`YoungOldHSCExp3$Transplant`
SpeciesCountMeans<-as.data.frame(t(SpeciesCountMeans[,2:ncol(SpeciesCountMeans)]))
SpeciesCountMeans$Species<-row.names(SpeciesCountMeans)
OverallMeans<-as.data.frame(colMeans(YoungOldHSCExp3[3:ncol(YoungOldHSCExp3)]))
OverallMeans$Species<-row.names(OverallMeans)
names(OverallMeans)<-c("OverallMean", "Species")
SpeciesCountMeans<-merge(SpeciesCountMeans, OverallMeans, by = "Species", all.x = TRUE)
SpeciesCountMeans$Young_HSCs<-as.numeric(SpeciesCountMeans$Young_HSCs)
SpeciesCountMeans$Old_HSCs<-as.numeric(SpeciesCountMeans$Old_HSCs)
SpeciesCountMeans$FoldChange<-(SpeciesCountMeans$Young_HSCs-SpeciesCountMeans$Old_HSCs) / SpeciesCountMeans$OverallMean
SpeciesCountMeans$Log2<-log2(SpeciesCountMeans$Young_HSCs/SpeciesCountMeans$OverallMean)
SpeciesCountTable<-merge(SpeciesCountMeans, SpeciesCounts, by = "Species", all.x = TRUE)
SpeciesCountTable<-merge(SpeciesCountTable, WilcoxTable, by = "Species", all.x = TRUE)
SpeciesCountTable<-SpeciesCountTable[order(SpeciesCountTable$Unadjusted_p),]
SpeciesCountTable<-SpeciesCountTable[,c(1,8,4, 2,3,5,9,10)]
names(SpeciesCountTable)<-c("Species", "Abundance", "OverallMean", "Young_HSCs.Mean", "Old_HSCs.Mean", "Fold.Difference", "p_Unadjusted", "FDR")
SpeciesCountTable$Abundance<-SpeciesCountTable$Abundance * 100
SpeciesCountTable$Fold.Difference[is.na(SpeciesCountTable$Fold.Difference)]<-0
SpeciesCountTable<-SpeciesCountTable[order(SpeciesCountTable$Fold.Difference, decreasing = TRUE),]

SigSpeciesCountTable<-subset(SpeciesCountTable, SpeciesCountTable$p_Unadjusted < 0.05)
SigSpeciesCountTable$HigherIn<-ifelse(SigSpeciesCountTable$Fold.Difference > 0, "Young_HSCs", "Old_HSCs")
SigSpeciesCountTable$HigherIn<-factor(SigSpeciesCountTable$HigherIn, levels = c("Young_HSCs", "Old_HSCs"))
write.csv(SigSpeciesCountTable, file = "SignificantSpeciesYoungvsOldHSC.csv")
SigYoungOldHSCExp2Species<-SigSpeciesCountTable$Species


metadata <- YoungOldHSCExp3[, 1:2]
cts <- as.matrix(YoungOldHSCExp3[, -(1:2)])
rownames(cts) <- metadata$SampleID
grps<-metadata$Transplant

cts_l2 <- glog2(cts)

res <- compute_ef_YoungOldHSCs(d = cts_l2, g = YoungOldHSCExp3$Transplant, min_shrink = 0.3)
names(res)[1]<-"Species"

res <- arrange(left_join(res, SigSpeciesCountTable), desc(abs(ef_shrunk))) %>%
  mutate(Species = as_factor(Species))

res<-subset(res, ! is.na(res$FDR))

res$HigherIn<-ifelse(res$ef_shrunk > 0, "Young_HSCs", "Old_HSCs")
res$HigherIn<-factor(res$HigherIn, levels = c("Young_HSCs", "Old_HSCs"))

resYoungOldHSCExp3s<-res

# write.table(SigSpeciesCountTable, file = "SignificantBackgroundSpecies.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)
# write.table(SpeciesCountTable, file = "GroupGeigerSpeciesNRK2.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)

col=tableau.colors[1:2]

YoungOldHSCExp3Effect<- ggplot(filter(resYoungOldHSCExp3s, resYoungOldHSCExp3s$ef_shrunk > 1 | resYoungOldHSCExp3s$ef_shrunk < -1),
                               aes(x = reorder(Species, ef_shrunk), y = ef_shrunk, fill = factor(HigherIn))) +
  geom_bar(colour="black", stat="identity") +
  #guides(fill=FALSE) + scale_fill_hue(l=40) +
  guides(fill=guide_legend(override.aes=list(colour=col))) +
  coord_flip() +
  labs( x = NULL) +
  ylab(expression(atop("Effect Size"))) +
  guides(fill=guide_legend(title="More abundant in:")) +
  #geom_hline(yintercept = c(-0.2, 0.2), color = "firebrick", lty = 2) +
  theme(axis.text.y = element_text(size= 12, color="black")) +
  theme(legend.title = element_text(size=12)) +
  theme(legend.text = element_text(size = 11)) 


YoungOldHSCExp3Effect<- YoungOldHSCExp3Effect + scale_fill_manual(values = col)

ggsave(filename = "YoungvsOldHSCLExpt3_Species_Effect.pdf", plot = YoungOldHSCExp3Effect, width = 12,
       height = 16, limitsize = FALSE)



#########################################################


col<-tableau.colors
Prevotella.sp.CAG.1124<-ggplot(DummySpecies,
                              aes(x=Groups, y=Prevotella.sp.CAG.1124, fill=NULL)) + geom_boxplot(lwd=1,aes(color=factor(Groups),fill = NULL), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Groups))) + xlab(NULL) + theme_bw() +facet_grid(facets = . ~ Experiment) +
  ylab("Prevotella.sp.CAG.1124\n")  +  paramsAngled() + scale_color_tableau() 
ggsave(filename = "Prevotella.sp.CAG.1124_not_good_result.pdf", plot = Prevotella.sp.CAG.1124, width = 10,  height = 8, limitsize = FALSE)

Prevotella.sp.CAG.487<-ggplot(DummySpecies,
                               aes(x=Groups, y=Prevotella.sp.CAG.487, fill=NULL)) + geom_boxplot(lwd=1,aes(color=factor(Groups),fill = NULL), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Groups))) + xlab(NULL) + theme_bw() +facet_grid(facets = . ~ Experiment) +
  ylab("Prevotella.sp.CAG.487\n")  +  paramsAngled() + scale_color_tableau() #+ scale_y_log10()
ggsave(filename = "Prevotella.sp.CAG.487_not_good_result.pdf", plot = Prevotella.sp.CAG.487, width = 10,  height = 8, limitsize = FALSE)



Odoribacter.sp.AF15.53<-ggplot(DummySpecies,
                               aes(x=Groups, y=Odoribacter.sp.AF15.53, fill=NULL)) + geom_boxplot(lwd=1,aes(color=factor(Groups),fill = NULL), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Groups))) + xlab(NULL) + theme_bw() +facet_grid(facets = . ~ Experiment) +
  ylab("Odoribacter.sp.AF15.53\n")  +  paramsAngled() + scale_color_tableau() #+ #scale_y_log10()
ggsave(filename = "Odoribacter.sp.AF15.53_pretty_good_result.pdf", plot = Odoribacter.sp.AF15.53, width = 10,  height = 8, limitsize = FALSE)

Alistipes.sp.Marseille.P2431<-ggplot(DummySpecies,
                               aes(x=Groups, y=Alistipes.sp.Marseille.P2431, fill=NULL)) + geom_boxplot(lwd=1,aes(color=factor(Groups),fill = NULL), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Groups))) + xlab(NULL) + theme_bw() +facet_grid(facets = . ~ Experiment) +
  ylab("Alistipes.sp.Marseille.P2431\n")  +  paramsAngled() + scale_color_tableau() + scale_y_log10()
ggsave(filename = "Alistipes.sp.Marseille.P2431.pdf_ok_result.pdf", plot = Alistipes.sp.Marseille.P2431, width = 10,  height = 8, limitsize = FALSE)

Butyricimonas.sp.Marseille.P4593<-ggplot(DummySpecies,
                                     aes(x=Groups, y=Butyricimonas.sp.Marseille.P4593, fill=NULL)) + geom_boxplot(lwd=1,aes(color=factor(Groups),fill = NULL), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Groups))) + xlab(NULL) + theme_bw() +facet_grid(facets = . ~ Experiment) +
  ylab("Butyricimonas.sp.Marseille.P4593\n")  +  paramsAngled() + scale_color_tableau() 
ggsave(filename = "Butyricimonas.sp.Marseille.P4593_not_good_result.pdf", plot = Butyricimonas.sp.Marseille.P4593, width = 10,  height = 8, limitsize = FALSE)

Bacteroides.coprophilus<-ggplot(DummySpecies,
                                         aes(x=Groups, y=Bacteroides.coprophilus, fill=NULL)) + geom_boxplot(lwd=1,aes(color=factor(Groups),fill = NULL), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Groups))) + xlab(NULL) + theme_bw() +facet_grid(facets = . ~ Experiment) +
  ylab("Bacteroides.coprophilus\n")  +  paramsAngled() + scale_color_tableau() #+ scale_y_log10()
ggsave(filename = "Bacteroides.coprophilus.pdf", plot = Bacteroides.coprophilus, width = 10,  height = 8, limitsize = FALSE)

Bacteroides fragilis

Bacteroides.fragilis<-ggplot(DummySpecies,
                                aes(x=Groups, y=Bacteroides.fragilis, fill=NULL)) + geom_boxplot(lwd=1,aes(color=factor(Groups),fill = NULL), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Groups))) + xlab(NULL) + theme_bw() +facet_grid(facets = . ~ Experiment) +
  ylab("Bacteroides.fragilis\n")  +  paramsAngled() + scale_color_tableau() #+ scale_y_log10()
ggsave(filename = "Bacteroides.fragilis.pdf", plot = Bacteroides.fragilis, width = 10,  height = 8, limitsize = FALSE)

Prevotella.copri<-ggplot(DummySpecies,
                             aes(x=Groups, y=Prevotella.copri, fill=NULL)) + geom_boxplot(lwd=1,aes(color=factor(Groups),fill = NULL), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Groups))) + xlab(NULL) + theme_bw() +facet_grid(facets = . ~ Experiment) +
  ylab("Prevotella.copri\n")  +  paramsAngled() + scale_color_tableau() #+ scale_y_log10()
ggsave(filename = "Prevotella.copri.pdf", plot = Prevotella.copri, width = 10,  height = 8, limitsize = FALSE)

Prevotella.copri<-ggplot(DummySpecies,
                         aes(x=Groups, y=Prevotella.copri, fill=NULL)) + geom_boxplot(lwd=1,aes(color=factor(Groups),fill = NULL), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Groups))) + xlab(NULL) + theme_bw() +facet_grid(facets = . ~ Experiment) +
  ylab("Prevotella.copri\n")  +  paramsAngled() + scale_color_tableau() #+ scale_y_log10()
ggsave(filename = "Prevotella.copri.pdf", plot = Prevotella.copri, width = 10,  height = 8, limitsize = FALSE)

Bifidobacterium.longum<-ggplot(DummySpecies,
                         aes(x=Groups, y=Bifidobacterium.longum, fill=NULL)) + geom_boxplot(lwd=1,aes(color=factor(Groups),fill = NULL), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Groups))) + xlab(NULL) + theme_bw() +facet_grid(facets = . ~ Experiment) +
  ylab("Bifidobacterium.longum\n")  +  paramsAngled() + scale_color_tableau() #+ scale_y_log10()
ggsave(filename = "Bifidobacterium.longum.pdf", plot = Bifidobacterium.longum, width = 10,  height = 8, limitsize = FALSE)

Bifidobacterium.longum<-ggplot(DummySpecies,
                               aes(x=Groups, y=Bifidobacterium.longum, fill=NULL)) + geom_boxplot(lwd=1,aes(color=factor(Groups),fill = NULL), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Groups))) + xlab(NULL) + theme_bw() +facet_grid(facets = . ~ Experiment) +
  ylab("Bifidobacterium.longum\n")  +  paramsAngled() + scale_color_tableau() #+ scale_y_log10()
ggsave(filename = "Bifidobacterium.longum.pdf", plot = Bifidobacterium.longum, width = 10,  height = 8, limitsize = FALSE)


Collinsella.aerofaciens<-ggplot(DummySpecies,
                               aes(x=Groups, y=Collinsella.aerofaciens, fill=NULL)) + geom_boxplot(lwd=1,aes(color=factor(Groups),fill = NULL), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Groups))) + xlab(NULL) + theme_bw() +facet_grid(facets = . ~ Experiment) +
  ylab("Collinsella.aerofaciens\n")  +  paramsAngled() + scale_color_tableau() #+ scale_y_log10()
ggsave(filename = "Collinsella.aerofaciens.pdf", plot = Collinsella.aerofaciens, width = 10,  height = 8, limitsize = FALSE)


Prevotella.copri<-ggplot(DummySpecies,
                                 aes(x=Groups, y=Prevotella.copri, fill=NULL)) + geom_boxplot(lwd=1,aes(color=factor(Groups),fill = NULL), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Groups))) + xlab(NULL) + theme_bw()  +facet_grid(facets = . ~ Experiment) +
  ylab("Prevotella.copri\n")  +  paramsAngled() + scale_color_tableau() #+ scale_y_log10()
ggsave(filename = "Prevotella.copri.pdf", plot = Prevotella.copri, width = 10,  height = 8, limitsize = FALSE)


Prevotella.copri.allExpt<-ggplot(DummySpecies,
                         aes(x=Groups, y=Prevotella.copri, fill=NULL)) + geom_boxplot(lwd=1,aes(color=factor(Groups),fill = NULL), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Groups))) + xlab(NULL) + theme_bw()  +
  ylab("Prevotella.copri\n")  +  paramsAngled() + scale_color_tableau() #+ scale_y_log10()
ggsave(filename = "Prevotella.copri.allExpt.pdf", plot = Prevotella.copri.allExpt, width = 5,  height = 8, limitsize = FALSE)


Helicobacter.pylori<-ggplot(DummySpecies,
                                aes(x=Groups, y=Helicobacter.pylori, fill=NULL)) + geom_boxplot(lwd=1,aes(color=factor(Groups),fill = NULL), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Groups))) + xlab(NULL) + theme_bw() +facet_grid(facets = . ~ Experiment) +
  ylab("Helicobacter.pylori\n")  +  paramsAngled() + scale_color_tableau() #+ scale_y_log10()
ggsave(filename = "Helicobacter.pylori.pdf", plot = Helicobacter.pylori, width = 10,  height = 8, limitsize = FALSE)


Collinsella.aerofaciens.allExpt<-ggplot(DummySpecies,
                                aes(x=Groups, y=Collinsella.aerofaciens, fill=NULL)) + geom_boxplot(lwd=1,aes(color=factor(Groups),fill = NULL), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Groups))) + xlab(NULL) + theme_bw()  +
  ylab("Collinsella.aerofaciens\n")  +  paramsAngled() + scale_color_tableau() #+ scale_y_log10()
ggsave(filename = "Collinsella.aerofaciens.allExp.pdf", plot = Collinsella.aerofaciens.allExpt, width = 5,  height = 8, limitsize = FALSE)



save.image(file = "GeigerData20220926")

### Alrighty, we're going to do difference between young and old using all groups. First untransplanted C57BL


UntransplantedSpecies<-subset(GeigerSpeciesNR, GeigerSpeciesNR$Groups %in% c("Young_C57BL-Untransplanted", "Old_C57BL-Untransplanted"))


UntransplantedSpeciesTable<-as.data.frame(t(aggregate( . ~ UntransplantedSpecies$Groups, data =UntransplantedSpecies[,11:ncol(UntransplantedSpecies)], FUN = mean)))
colnames(UntransplantedSpeciesTable)<-UntransplantedSpeciesTable[1,]
UntransplantedSpeciesTable<-UntransplantedSpeciesTable[-1,]
UntransplantedSpeciesTable<-as.data.frame(data.matrix(UntransplantedSpeciesTable))
UntransplantedSpeciesTable$Species<-row.names(UntransplantedSpeciesTable)

UntransplantedSpeciesTable$Difference<-(UntransplantedSpeciesTable$`Young_C57BL-Untransplanted` - UntransplantedSpeciesTable$`Old_C57BL-Untransplanted`)
UntransplantedSpeciesTable$FC<-foldchange(UntransplantedSpeciesTable$`Young_C57BL-Untransplanted`, UntransplantedSpeciesTable$`Old_C57BL-Untransplanted`)
UntransplantedSpeciesTable$logratio<-foldchange2logratio(UntransplantedSpeciesTable$FC)
# UntransplantedSpeciesTable$log2Delta<-ifelse(sign(UntransplantedSpeciesTable$Difference) == -1, -UntransplantedSpeciesTable$log2Delta,
#                                       ifelse(UntransplantedSpeciesTable$Difference == 0, 0, UntransplantedSpeciesTable$log2Delta))
UntransplantedSpeciesTable<-UntransplantedSpeciesTable[order(UntransplantedSpeciesTable$logratio, decreasing = FALSE),]
UntransplantedSpeciesTable$Species<-row.names(UntransplantedSpeciesTable)

Kruskal<-kruskal.test(UntransplantedSpecies[,11], factor(UntransplantedSpecies[,8]))
KruskalTable<-as.data.frame(Kruskal$p.value)
names(KruskalTable)<-names(UntransplantedSpecies)[11]
row.names(KruskalTable)<-"Kruskal"

for (i in 12:length(colnames(UntransplantedSpecies))){
  x<-kruskal.test(UntransplantedSpecies[,i], factor(UntransplantedSpecies[,8]))
  y<-as.data.frame(x$p.value)
  names(y)<-names(UntransplantedSpecies)[i]
  KruskalTable<-cbind(KruskalTable, y)
}

KruskalTable<-as.data.frame(t(KruskalTable))
names(KruskalTable)<-c("Unadjusted_p")
KruskalTable$FDR<-p.adjust(KruskalTable$Unadjusted_p, method = "fdr")
KruskalTable$Species<-row.names(KruskalTable)
KruskalTable<-KruskalTable[order(KruskalTable$Unadjusted_p),]

UntransplantedSpeciesTable<-merge(UntransplantedSpeciesTable, KruskalTable, by = "Species")


UntransplantedSpeciesTable2<-subset(UntransplantedSpeciesTable, UntransplantedSpeciesTable$Unadjusted_p < 0.05)
UntransplantedSpeciesTable2<-UntransplantedSpeciesTable2[order(UntransplantedSpeciesTable2$logratio, decreasing = FALSE),]
UntransplantedSpeciesTable3<-subset(UntransplantedSpeciesTable, UntransplantedSpeciesTable$FDR < 0.15)
UntransplantedPathwayTable<-UntransplantedSpeciesTable[order(UntransplantedSpeciesTable$logratio, decreasing = FALSE),]


SigYoungOldNonTransplantedSpeciesKruskal<-UntransplantedSpeciesTable2$Species
intersect(SigYoungOldNonTransplantedSpecies, SigYoungOldNonTransplantedSpeciesKruskal)
union(SigYoungOldNonTransplantedSpecies, SigYoungOldNonTransplantedSpeciesKruskal)
setdiff(SigYoungOldNonTransplantedSpecies, SigYoungOldNonTransplantedSpeciesKruskal)
setdiff(SigYoungOldNonTransplantedSpeciesKruskal, SigYoungOldNonTransplantedSpecies)



setwd("C:/Users/dbhas/OneDrive/Documents/Code/Metagenomics/GeigerData/StemCellsAgeing")
setwd("/home/david/Documents/Code/Metagenomics/GeigerData/StemCellsAgeing")
setwd("C:/Users/HASI9S/Documents/Code/Metagenomics/GeigerData/StemCellsAgeing/")

write.table(UntransplantedSpeciesTable2, file = "SignificantSpeciesYoungVsOldUntransplantedKruskall20220910.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)

#### Now transplanted young vs old HSC

TransplantedSpecies<-subset(GeigerSpeciesNR, GeigerSpeciesNR$Groups %in% c("RAG1-Young_HSCs", "RAG1-Old_HSCs"))


TransplantedSpeciesTable<-as.data.frame(t(aggregate( . ~ TransplantedSpecies$Groups, data =TransplantedSpecies[,11:ncol(TransplantedSpecies)], FUN = mean)))
colnames(TransplantedSpeciesTable)<-TransplantedSpeciesTable[1,]
TransplantedSpeciesTable<-TransplantedSpeciesTable[-1,]
TransplantedSpeciesTable<-as.data.frame(data.matrix(TransplantedSpeciesTable))
TransplantedSpeciesTable$Species<-row.names(TransplantedSpeciesTable)

TransplantedSpeciesTable$Difference<-(TransplantedSpeciesTable$`RAG1-Young_HSCs` - TransplantedSpeciesTable$`RAG1-Old_HSCs`)
TransplantedSpeciesTable$FC<-foldchange(TransplantedSpeciesTable$`RAG1-Young_HSCs`, TransplantedSpeciesTable$`RAG1-Old_HSCs`)
TransplantedSpeciesTable$logratio<-foldchange2logratio(TransplantedSpeciesTable$FC)
# TransplantedSpeciesTable$log2Delta<-ifelse(sign(TransplantedSpeciesTable$Difference) == -1, -TransplantedSpeciesTable$log2Delta,
#                                       ifelse(TransplantedSpeciesTable$Difference == 0, 0, TransplantedSpeciesTable$log2Delta))
TransplantedSpeciesTable<-TransplantedSpeciesTable[order(TransplantedSpeciesTable$logratio, decreasing = FALSE),]
TransplantedSpeciesTable$Species<-row.names(TransplantedSpeciesTable)

Kruskal<-kruskal.test(TransplantedSpecies[,11], factor(TransplantedSpecies[,8]))
KruskalTable<-as.data.frame(Kruskal$p.value)
names(KruskalTable)<-names(TransplantedSpecies)[11]
row.names(KruskalTable)<-"Kruskal"

for (i in 12:length(colnames(TransplantedSpecies))){
  x<-kruskal.test(TransplantedSpecies[,i], factor(TransplantedSpecies[,8]))
  y<-as.data.frame(x$p.value)
  names(y)<-names(TransplantedSpecies)[i]
  KruskalTable<-cbind(KruskalTable, y)
}

KruskalTable<-as.data.frame(t(KruskalTable))
names(KruskalTable)<-c("Unadjusted_p")
KruskalTable$FDR<-p.adjust(KruskalTable$Unadjusted_p, method = "fdr")
KruskalTable$Species<-row.names(KruskalTable)
KruskalTable<-KruskalTable[order(KruskalTable$Unadjusted_p),]

TransplantedSpeciesTable<-merge(TransplantedSpeciesTable, KruskalTable, by = "Species")


TransplantedSpeciesTable2<-subset(TransplantedSpeciesTable, TransplantedSpeciesTable$Unadjusted_p < 0.05)
TransplantedSpeciesTable2<-TransplantedSpeciesTable2[order(TransplantedSpeciesTable2$logratio, decreasing = FALSE),]
TransplantedSpeciesTable3<-subset(TransplantedSpeciesTable, TransplantedSpeciesTable$FDR < 0.15)
TransplantedPathwayTable<-TransplantedSpeciesTable[order(TransplantedSpeciesTable$logratio, decreasing = FALSE),]



SigYoungOldHSCSpeciesKruskal<-TransplantedSpeciesTable2$Species
intersect(SigYoungOldHSCSpeciesKruskal, SigYoungOldNonTransplantedSpeciesKruskal)



setwd("C:/Users/dbhas/OneDrive/Documents/Code/Metagenomics/GeigerData/StemCellsAgeing")
setwd("/home/david/Documents/Code/Metagenomics/GeigerData/StemCellsAgeing")
setwd("C:/Users/HASI9S/Documents/Code/Metagenomics/GeigerData/StemCellsAgeing/")

write.table(TransplantedSpeciesTable2, file = "SignificantSpeciesYoungVsOldTransplantedKruskal20220910.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)




IntersectingSpecies<-intersect(TransplantedSpeciesTable2$Species, UntransplantedSpeciesTable2$Species)
write.csv(IntersectingSpecies, file = "IntersectingSpeciesYvsOHSC_YvsOTransplanted.csv")


## ok make some graphs

# DummySpecies<-GeigerSpeciesNR
# 
# col<-tableau.colors
# Paenibacillus_Groups<-ggplot(DummySpecies,
#                               aes(x=Groups, y=Paenibacillus.terrigena, fill=Groups)) + geom_boxplot(lwd=1,aes(color=factor(Groups),fill = NA), outlier.size = 3)  +
#   stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
#   geom_point(size=4,aes(color = factor(Groups))) + xlab(NULL) +
#   ylab("Paenibacillus.terrigena\n")  +  paramsAngled() + scale_color_tableau() 
# ggsave(filename = "Paenibacillus_Groups.pdf", plot = Paenibacillus_Groups, width = 5,  height = 8, limitsize = FALSE)
# 
# 
# Lactobacillus.paragasseri_Groups<-ggplot(DummySpecies,
#                              aes(x=Groups, y=Lactobacillus.paragasseri, fill=Groups)) + geom_boxplot(lwd=1,aes(color=factor(Groups),fill = NA), outlier.size = 3)  +
#   stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
#   geom_point(size=4,aes(color = factor(Groups))) + xlab(NULL) +
#   ylab("Lactobacillus.paragasseri\n")  +  paramsAngled() + scale_color_tableau() 
# ggsave(filename = "Lactobacillus.paragasseri_Groups.pdf", plot = Lactobacillus.paragasseri_Groups, width = 5,  height = 8, limitsize = FALSE)
# 
# Lactobacillus.hamsteri_Groups<-ggplot(DummySpecies,
#                                          aes(x=Groups, y=Lactobacillus.hamsteri, fill=Groups)) + geom_boxplot(lwd=1,aes(color=factor(Groups),fill = NA), outlier.size = 3)  +
#   stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
#   geom_point(size=4,aes(color = factor(Groups))) + xlab(NULL) +
#   ylab("Lactobacillus.hamsteri\n")  +  paramsAngled() + scale_color_tableau() 
# ggsave(filename = "Lactobacillus.hamsteri_Groups.pdf", plot = Lactobacillus.hamsteri_Groups, width = 5,  height = 8, limitsize = FALSE)
# 
# Rufibacter.immobilis_Groups<-ggplot(DummySpecies,
#                                       aes(x=Groups, y=Rufibacter.immobilis, fill=Groups)) + geom_boxplot(lwd=1,aes(color=factor(Groups),fill = NA), outlier.size = 3)  +
#   stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
#   geom_point(size=4,aes(color = factor(Groups))) + xlab(NULL) +
#   ylab("Rufibacter.immobilis\n")  +  paramsAngled() + scale_color_tableau() 
# ggsave(filename = "Rufibacter.immobilis_Groups.pdf", plot = Rufibacter.immobilis_Groups, width = 5,  height = 8, limitsize = FALSE)
# 
# 
# Prevotella.bivia_Groups<-ggplot(DummySpecies,
#                                     aes(x=Groups, y=Prevotella.bivia, fill=Groups)) + geom_boxplot(lwd=1,aes(color=factor(Groups),fill = NA), outlier.size = 3)  +
#   stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
#   geom_point(size=4,aes(color = factor(Groups))) + xlab(NULL) +
#   ylab("Prevotella.bivia\n")  +  paramsAngled() + scale_color_tableau() 
# ggsave(filename = "Prevotella.bivia_Groups.pdf", plot = Prevotella.bivia_Groups, width = 5,  height = 8, limitsize = FALSE)
# 

save.image(file= "GeigerData20201230")



## Ok here are the steps for getting metabolic associations:

# see the script: Humann3ScriptsGeiger.sh


# setwd("P:/Documents/Alignments/Humann2Alignments/PathwayAbundance/GeigerSamplesNew")
# setwd("~/Documents/Alignments/Humann2Alignments/PathwayAbundance/GeigerSamplesNew")
NormalizedPathway<-read.csv("GeigerFiles_pathabundance20220910-cpm_unstratified.tsv",sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# Using stratified now

# NormalizedPathway<-read.csv("GeigerAssociated_pathabundance_stratified.tsv",sep = "\t", header = TRUE, stringsAsFactors = FALSE)

colnames(NormalizedPathway)[1]<-"Pathway"
row.names(NormalizedPathway)<-NormalizedPathway$Pathway
Group<-as.character(NormalizedPathway[1,-1])
NormalizedPathway<-NormalizedPathway[-1,]
NormalizedPathway$Pathway<-NULL

PathwayTable<-data.matrix(NormalizedPathway)
PathwayTable<-as.data.frame(t(PathwayTable))
PathwayTable$SampleID<-row.names(PathwayTable)
PathwayTable$SampleID<-gsub(".paired_Abundance.CPM", "", PathwayTable$SampleID)
PathwayTable<-merge(Metadata, PathwayTable, by = "SampleID", all.x = TRUE)
PathwayTable<-na.omit(PathwayTable)
# SampleCol<-(ncol(PathwayTable)-1)
# GeigerCol<-ncol(PathwayTable)
# GeigerDataCol<-1:(ncol(PathwayTable)-2)
# PathwayTable<-PathwayTable[, c(SampleCol, GeigerCol,GeigerDataCol)]
head(names(PathwayTable),12)
# Get rid of unintegrated and unmapped:
PathwayTable$UNINTEGRATED<-NULL
PathwayTable$UNMAPPED<-NULL


# ok now have to normalize: how about doing normalized log versus rarefaction

GeigerPathwayNR<-PathwayTable[,11:ncol(PathwayTable)]
GeigerPathwayNR<-GeigerPathwayNR*10
GeigerPathwayNR<-floor(GeigerPathwayNR)
minrow<-min(rowSums(GeigerPathwayNR))
GeigerPathwayNR<-as.data.frame(rrarefy(GeigerPathwayNR, minrow))
GeigerPathwayNR<-cbind(PathwayTable[,1:10], GeigerPathwayNR)

GeigerPathway.mrpp<-mrpp(GeigerPathwayNR[11:ncol(GeigerPathwayNR)], GeigerPathwayNR[,8])
# Ok p = 0.185

# let's calculate fold change and also significance of difference between Young_HSV vs Old_HSC
TransplantPathways<-subset(GeigerPathwayNR, GeigerPathwayNR$Groups %in% c("RAG1-Young_HSCs", "RAG1-Old_HSCs"))

GeigerPathwayTable<-as.data.frame(t(aggregate( . ~ TransplantPathways$Groups, data =TransplantPathways[,11:ncol(TransplantPathways)], FUN = mean)))
colnames(GeigerPathwayTable)<-GeigerPathwayTable[1,]
GeigerPathwayTable<-GeigerPathwayTable[-1,]
GeigerPathwayTable<-as.data.frame(data.matrix(GeigerPathwayTable))
GeigerPathwayTable$Pathway<-row.names(GeigerPathwayTable)

GeigerPathwayTable$Difference<-(GeigerPathwayTable$`RAG1-Young_HSCs` -GeigerPathwayTable$`RAG1-Old_HSCs`)
GeigerPathwayTable$FC<-foldchange(GeigerPathwayTable$`RAG1-Young_HSCs`, GeigerPathwayTable$`RAG1-Old_HSCs`)
GeigerPathwayTable$logratio<-foldchange2logratio(GeigerPathwayTable$FC)
# GeigerPathwayTable$log2Delta<-ifelse(sign(GeigerPathwayTable$Difference) == -1, -GeigerPathwayTable$log2Delta,
#                                       ifelse(GeigerPathwayTable$Difference == 0, 0, GeigerPathwayTable$log2Delta))
GeigerPathwayTable<-GeigerPathwayTable[order(GeigerPathwayTable$logratio, decreasing = FALSE),]
GeigerPathwayTable$Pathway<-row.names(GeigerPathwayTable)

# significant pathways, RAG1-Old_HSC vs RAG1-Young_HSC
TransplantPathways<-subset(GeigerPathwayNR, GeigerPathwayNR$Groups %in% c("RAG1-Young_HSCs", "RAG1-Old_HSCs"))

Kruskal<-kruskal.test(TransplantPathways[,11], factor(TransplantPathways[,8]))
KruskalTable<-as.data.frame(Kruskal$p.value)
names(KruskalTable)<-names(TransplantPathways)[11]
row.names(KruskalTable)<-"Kruskal"

for (i in 12:length(colnames(TransplantPathways))){
  x<-kruskal.test(TransplantPathways[,i], factor(TransplantPathways[,8]))
  y<-as.data.frame(x$p.value)
  names(y)<-names(TransplantPathways)[i]
  KruskalTable<-cbind(KruskalTable, y)
}

KruskalTable<-as.data.frame(t(KruskalTable))
names(KruskalTable)<-c("Unadjusted_p")
KruskalTable$FDR<-p.adjust(KruskalTable$Unadjusted_p, method = "fdr")
KruskalTable$Pathway<-row.names(KruskalTable)
KruskalTable<-KruskalTable[order(KruskalTable$Unadjusted_p),]

GeigerPathwayTable<-merge(GeigerPathwayTable, KruskalTable, by = "Pathway", all.x = TRUE)

# remove weird pathways
Engineered<-grep("engineered", GeigerPathwayTable$Pathway)
Fungi<-grep("fungi", GeigerPathwayTable$Pathway)
Animal<-grep("animal", GeigerPathwayTable$Pathway)
Yeast<-grep("yeast", GeigerPathwayTable$Pathway)
Weird<-c(Engineered, Fungi, Animal, Yeast)

GeigerPathwayTable<-GeigerPathwayTable[-Weird,]


GeigerPathwayTable2<-subset(GeigerPathwayTable, GeigerPathwayTable$Unadjusted_p < 0.05)
GeigerPathwayTable2<-GeigerPathwayTable2[order(GeigerPathwayTable2$logratio, decreasing = FALSE),]
GeigerPathwayTable3<-subset(GeigerPathwayTable, GeigerPathwayTable$FDR < 0.15)
GeigerPathwayTable<-GeigerPathwayTable[order(GeigerPathwayTable$logratio, decreasing = FALSE),]

GeigerPathwayTable$Pathway<-gsub("\t", " ", GeigerPathwayTable$Pathway)
GeigerPathwayTable$Pathway<-gsub("   ", " ", GeigerPathwayTable$Pathway)
GeigerPathwayTable$Pathway<-gsub("  ", " ", GeigerPathwayTable$Pathway)
GeigerPathwayTable$Pathway<-gsub(",", "", GeigerPathwayTable$Pathway)

setwd("C:/Users/dbhas/OneDrive/Documents/Code/Metagenomics/GeigerData/StemCellsAgeing")
setwd("/home/david/Documents/Code/Metagenomics/GeigerData/StemCellsAgeing")
setwd("C:/Users/HASI9S/Documents/Code/Metagenomics/GeigerData/StemCellsAgeing/")

write.table(GeigerPathwayTable2, file = "SignificantPathwaysYoungVsOldHSC20220927.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)



# make a heatmap of logratio:

matrixtable<-GeigerPathwayTable2[,6:7]
row.names(matrixtable)<-GeigerPathwayTable2$Pathway
matrixtable<-as.matrix(matrixtable)
heatmap.2(matrixtable, Rowv = NULL, Colv = FALSE, col = diverging.colors(100)[50:100], trace = "none", 
          tracecol= "black",hline = NULL, vline = NULL, margins=c(2,25),
          xlab = NULL, main = NULL)


matrixtable<-GeigerPathwayTable2[,6:7]
row.names(matrixtable)<-GeigerPathwayTable2$Pathway
matrixtable<-as.matrix(matrixtable)
heatmap.2(matrixtable, Rowv = NULL, Colv = FALSE, col = diverging.colors(100), trace = "none", 
          tracecol= "black",hline = NULL, vline = NULL, margins=c(2,25),
          xlab = NULL, main = NULL)


# make a heatmap of logratio:
# 
# matrixtable<-TransplantedPathwayTable2[,6:7]
# row.names(matrixtable)<-TransplantedPathwayTable2$Pathway
# matrixtable<-as.matrix(matrixtable)
# heatmap.2(matrixtable, Rowv = NULL, Colv = FALSE, col = diverging.colors(100), trace = "none", 
#           tracecol= "black",hline = NULL, vline = NULL, margins=c(2,25),
#           xlab = NULL, main = NULL)
# ##################






#### Now Young versus Old Untransplanted ####### 

UntransplantedPathways<-subset(GeigerPathwayNR, GeigerPathwayNR$Groups %in% c("Old_C57BL-Untransplanted", "Young_C57BL-Untransplanted"))

UntransplantedPathwayTable<-as.data.frame(t(aggregate( . ~ UntransplantedPathways$Groups, data =UntransplantedPathways[,11:ncol(UntransplantedPathways)], FUN = mean)))
colnames(UntransplantedPathwayTable)<-UntransplantedPathwayTable[1,]
UntransplantedPathwayTable<-UntransplantedPathwayTable[-1,]
UntransplantedPathwayTable<-as.data.frame(data.matrix(UntransplantedPathwayTable))
UntransplantedPathwayTable$Pathway<-row.names(UntransplantedPathwayTable)

UntransplantedPathwayTable$Difference<-(UntransplantedPathwayTable$`Young_C57BL-Untransplanted` -UntransplantedPathwayTable$`Old_C57BL-Untransplanted`)
UntransplantedPathwayTable$FC<-foldchange(UntransplantedPathwayTable$`Young_C57BL-Untransplanted`, UntransplantedPathwayTable$`Old_C57BL-Untransplanted`)
UntransplantedPathwayTable$logratio<-foldchange2logratio(UntransplantedPathwayTable$FC)
# UntransplantedPathwayTable$log2Delta<-ifelse(sign(UntransplantedPathwayTable$Difference) == -1, -UntransplantedPathwayTable$log2Delta,
#                                       ifelse(UntransplantedPathwayTable$Difference == 0, 0, UntransplantedPathwayTable$log2Delta))
UntransplantedPathwayTable<-UntransplantedPathwayTable[order(UntransplantedPathwayTable$logratio, decreasing = FALSE),]
UntransplantedPathwayTable$Pathway<-row.names(UntransplantedPathwayTable)

# significant pathways, RAG1-Old_HSC vs RAG1-Young_HSC
UntransplantedPathways<-subset(GeigerPathwayNR, GeigerPathwayNR$Groups %in% c("Young_C57BL-Untransplanted", "Old_C57BL-Untransplanted"))

Kruskal<-kruskal.test(UntransplantedPathways[,11], factor(UntransplantedPathways[,8]))
KruskalTable<-as.data.frame(Kruskal$p.value)
names(KruskalTable)<-names(UntransplantedPathways)[11]
row.names(KruskalTable)<-"Kruskal"

for (i in 12:length(colnames(UntransplantedPathways))){
  x<-kruskal.test(UntransplantedPathways[,i], factor(UntransplantedPathways[,8]))
  y<-as.data.frame(x$p.value)
  names(y)<-names(UntransplantedPathways)[i]
  KruskalTable<-cbind(KruskalTable, y)
}

KruskalTable<-as.data.frame(t(KruskalTable))
names(KruskalTable)<-c("Unadjusted_p")
KruskalTable$FDR<-p.adjust(KruskalTable$Unadjusted_p, method = "fdr")
KruskalTable$Pathway<-row.names(KruskalTable)
KruskalTable<-KruskalTable[order(KruskalTable$Unadjusted_p),]

UntransplantedPathwayTable<-merge(UntransplantedPathwayTable, KruskalTable, by = "Pathway", all.x = TRUE)

# remove weird pathways
Engineered<-grep("engineered", UntransplantedPathwayTable$Pathway)
Fungi<-grep("fungi", UntransplantedPathwayTable$Pathway)
Animal<-grep("animal", UntransplantedPathwayTable$Pathway)
Yeast<-grep("yeast", UntransplantedPathwayTable$Pathway)
Weird<-c(Engineered, Fungi, Animal, Yeast)

UntransplantedPathwayTable<-UntransplantedPathwayTable[-Weird,]


UntransplantedPathwayTable2<-subset(UntransplantedPathwayTable, UntransplantedPathwayTable$Unadjusted_p < 0.05)
UntransplantedPathwayTable2<-UntransplantedPathwayTable2[order(UntransplantedPathwayTable2$logratio, decreasing = FALSE),]
UntransplantedPathwayTable3<-subset(UntransplantedPathwayTable, UntransplantedPathwayTable$FDR < 0.15)
UntransplantedPathwayTable<-UntransplantedPathwayTable[order(UntransplantedPathwayTable$logratio, decreasing = FALSE),]

UntransplantedPathwayTable$Pathway<-gsub("\t", " ", UntransplantedPathwayTable$Pathway)
UntransplantedPathwayTable$Pathway<-gsub("   ", " ", UntransplantedPathwayTable$Pathway)
UntransplantedPathwayTable$Pathway<-gsub("  ", " ", UntransplantedPathwayTable$Pathway)
UntransplantedPathwayTable$Pathway<-gsub(",", "", UntransplantedPathwayTable$Pathway)

setwd("C:/Users/dbhas/OneDrive/Documents/Code/Metagenomics/GeigerData/StemCellsAgeing")
setwd("/home/david/Documents/Code/Metagenomics/GeigerData/StemCellsAgeing")
setwd("C:/Users/HASI9S/Documents/Code/Metagenomics/GeigerData/StemCellsAgeing/")

write.table(UntransplantedPathwayTable2, file = "SignificantPathwaysYoungVsOldUntransplanted20220927.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)

IntersectingPathways<-intersect(UntransplantedPathwayTable2$Pathway, GeigerPathwayTable2$Pathway)
write.csv(IntersectingPathways, file = "IntersectingYvsOHSC_YvsOUntransplanted.csv")

# make a heatmap of logratio:

matrixtable<-UntransplantedPathwayTable2[,6:7]
row.names(matrixtable)<-UntransplantedPathwayTable2$Pathway
matrixtable<-as.matrix(matrixtable)
heatmap.2(matrixtable, Rowv = NULL, Colv = FALSE, col = diverging.colors(100), trace = "none", 
          tracecol= "black",hline = NULL, vline = NULL, margins=c(2,25),
          xlab = NULL, main = NULL)
# let's calculate fold change and also significance of difference:

# UntransplantedPathwayTable<-as.data.frame(t(aggregate( . ~ Group, data =GeigerPathwayNR[,3:ncol(GeigerPathwayNR)], FUN = mean)))
# UntransplantedPathwayTable<-UntransplantedPathwayTable[-1,]
# UntransplantedPathwayTable<-as.data.frame(data.matrix(UntransplantedPathwayTable))
# colnames(UntransplantedPathwayTable)<-c("Null", "WT")
# 
# UntransplantedPathwayTable$Difference<-(UntransplantedPathwayTable$RAG1 -UntransplantedPathwayTable$C57BL) 
# UntransplantedPathwayTable$FC<-foldchange(UntransplantedPathwayTable$RAG1, UntransplantedPathwayTable$C57BL)
# UntransplantedPathwayTable$logratio<-foldchange2logratio(UntransplantedPathwayTable$FC)
# # UntransplantedPathwayTable$log2Delta<-ifelse(sign(UntransplantedPathwayTable$Difference) == -1, -UntransplantedPathwayTable$log2Delta, 
# #                                       ifelse(UntransplantedPathwayTable$Difference == 0, 0, UntransplantedPathwayTable$log2Delta))
# UntransplantedPathwayTable<-UntransplantedPathwayTable[order(UntransplantedPathwayTable$logratio, decreasing = FALSE),]
# UntransplantedPathwayTable$Pathway<-row.names(UntransplantedPathwayTable)
# 
# # Ok that's cool. Now figure out which are statistically signficant from the GeigerPathwayNR table
# Kruskal<-kruskal.test(GeigerPathwayNR[,3], factor(GeigerPathwayNR[,2]))
# KruskalTable<-as.data.frame(Kruskal$p.value)
# names(KruskalTable)<-names(GeigerPathwayNR)[3]
# row.names(KruskalTable)<-"Kruskal"
# 
# for (i in 4:length(colnames(GeigerPathwayNR))){
#   x<-kruskal.test(GeigerPathwayNR[,i], factor(GeigerPathwayNR[,2]))
#   y<-as.data.frame(x$p.value)
#   names(y)<-names(GeigerPathwayNR)[i]
#   KruskalTable<-cbind(KruskalTable, y)
# }
# 
# KruskalTable<-as.data.frame(t(KruskalTable))
# names(KruskalTable)<-c("Unadjusted_p")
# KruskalTable$FDR<-p.adjust(KruskalTable$Unadjusted_p, method = "fdr")
# KruskalTable$Pathway<-row.names(KruskalTable)
# KruskalTable<-KruskalTable[order(KruskalTable$Unadjusted_p),]
# 
# UntransplantedPathwayTable<-merge(UntransplantedPathwayTable, KruskalTable, by = "Pathway", all.x = TRUE)
# 
# UntransplantedPathwayTable2<-subset(UntransplantedPathwayTable, UntransplantedPathwayTable$Unadjusted_p < 0.05)
# UntransplantedPathwayTable2<-UntransplantedPathwayTable2[order(UntransplantedPathwayTable2$logratio, decreasing = FALSE),]
# UntransplantedPathwayTable3<-subset(UntransplantedPathwayTable, UntransplantedPathwayTable$FDR < 0.2)
# UntransplantedPathwayTable<-UntransplantedPathwayTable[order(UntransplantedPathwayTable$logratio, decreasing = FALSE),]
# 
# UntransplantedPathwayTable$Pathway<-gsub("\t", " ", UntransplantedPathwayTable$Pathway)
# UntransplantedPathwayTable$Pathway<-gsub("   ", " ", UntransplantedPathwayTable$Pathway)
# UntransplantedPathwayTable$Pathway<-gsub("  ", " ", UntransplantedPathwayTable$Pathway)
# UntransplantedPathwayTable$Pathway<-gsub(",", "", UntransplantedPathwayTable$Pathway)
# 
# setwd("C:/Users/dbhas/OneDrive/Documents/Code/Metagenomics/Geiger")
# setwd("~/Documents/Code/Metagenomics/Geiger")
# setwd("C:/Users/HASI9S/Documents/Code/Metagenomics/Geiger/")
# 
# write.table(UntransplantedPathwayTable, file = "UnstratifiedUntransplantedPathwayTable20190822.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
# 
# 
# 
# # make a heatmap of logratio:
# 
# matrixtable<-UntransplantedPathwayTable2[,6:7]
# row.names(matrixtable)<-UntransplantedPathwayTable2$Pathway
# matrixtable<-as.data.frame(matrixtable)
# matrixtable<-subset(matrixtable, matrixtable$logratio < -1 | matrixtable$logratio > 1)
# matrixtable<-as.matrix(matrixtable)
# heatmap.2(matrixtable, Rowv = NULL, Colv = FALSE, col = diverging.colors(100), trace = "none", 
#           tracecol= "black",hline = NULL, vline = NULL, margins=c(2,25),
#           xlab = NULL, main = NULL)
# 
# 
# heatmap.2(t(matrixtable), Rowv = NULL, Colv = FALSE, col = diverging.colors(100), trace = "none", 
#           tracecol= "black",hline = NULL, vline = NULL, margins=c(25,2),
#           xlab = NULL, main = NULL)


################################################################################################


# Genus level:

setwd("C:/Users/dbhas/OneDrive/Documents/Alignments/KrakenAlignments/Kraken2")
setwd("~/Documents/Alignments/KrakenAlignments/Kraken2")

AllKrakenFiles<-list.files()
GenusFileList<-grep("_genus_abundance.txt", AllKrakenFiles)
GenusFiles<-AllKrakenFiles[GenusFileList]
FileList<-gsub("_genus_abundance.txt", "", GenusFiles)

NewGenusFileList<-subset(FileList, FileList %in% GeigerSamples)

#setwd("~/Documents/Code/Metagenomics/HD")

# get the files
for(f in 1:length(NewGenusFileList)){
  fnr = NewGenusFileList[f]
  x =NewGenusFileList[f]
  
  # assign(fnr, read.csv(paste("/home/david/Databases/FecalABDR/NewKrakenAlignments/",x, "_Genus_abundance.txt",sep=""),
  #                      sep="\t", header = TRUE, stringsAsFactors = FALSE))
  
  assign(fnr, read.csv(paste(x, "_genus_abundance.txt",sep=""),
                       sep="\t", header = TRUE, stringsAsFactors = FALSE))
}


# NewGenusNR<-read.csv(paste("/home/david/Databases/FecalABDR/NewKrakenAlignments/",NewGenusFileList[1], "_Genus_abundance.txt", sep=""), sep = "\t",header = TRUE)
NewGenusNR<-read.csv(paste(NewGenusFileList[1], "_genus_abundance.txt", sep=""), sep = "\t",header = TRUE)

names(NewGenusNR)[1]<-"Genus"
NewGenusNR$Genus<-gsub(" ", ".", NewGenusNR$Genus, fixed = TRUE)
NewGenusNR$Genus<-gsub("..", ".", NewGenusNR$Genus, fixed = TRUE)
NewGenusNR$Genus<-gsub("_", ".", NewGenusNR$Genus, fixed = TRUE)
NewGenusNR$Genus<-gsub("-", ".", NewGenusNR$Genus, fixed = TRUE)
NewGenusNR$Genus<-gsub("/", ".", NewGenusNR$Genus, fixed = TRUE)
NewGenusNR$Genus<-gsub("X,", "", NewGenusNR$Genus, fixed = TRUE)
NewGenusNR$Genus<-gsub("[", "", NewGenusNR$Genus, fixed = TRUE)
NewGenusNR$Genus<-gsub("]", "", NewGenusNR$Genus, fixed = TRUE)
NewGenusNR<-as.data.frame(NewGenusNR[,c(1,4)])
names(NewGenusNR)<-c("Genus",paste(paste(NewGenusFileList[1])))
NewGenusNR<-subset(NewGenusNR, ! duplicated(NewGenusNR$Genus))

for ( i in 2:length(NewGenusFileList)){
  x <- get( NewGenusFileList[i] )
  x<-as.data.frame(x[,c(1,4)])
  names( x ) <- c("Genus", paste(NewGenusFileList[i]))
  x$Genus<-gsub(" ", ".", x$Genus, fixed = TRUE)
  x$Genus<-gsub("..", ".", x$Genus, fixed = TRUE)
  x$Genus<-gsub("_", ".", x$Genus, fixed = TRUE)
  x$Genus<-gsub("-", ".", x$Genus, fixed = TRUE)
  x$Genus<-gsub("/", ".", x$Genus, fixed = TRUE)
  x$Genus<-gsub("X,", "", x$Genus, fixed = TRUE)
  x$Genus<-gsub("[", "", x$Genus, fixed = TRUE)
  x$Genus<-gsub("]", "", x$Genus, fixed = TRUE)
  x<-subset(x, ! duplicated(x$Genus))
  NewGenusNR<-merge(x, NewGenusNR, by = "Genus", all = TRUE)
}

row.names(NewGenusNR)<-NewGenusNR$Genus
NewGenusNR$Genus<-NULL
NewGenusNR[is.na(NewGenusNR)]<-0

rm(list = FileList)


setwd("C:/Users/dbhas/OneDrive/Documents/Code/Metagenomics/GeigerData/StemCellsAgeing")
setwd("C:/Users/HASI9S/OneDrive/Documents/Code/Metagenomics/GeigerData/StemCellsAgeing")
setwd("~/Documents/Code/Metagenomics/GeigerData/StemCellsAgeing")

# Remove samples in less than 10% of samples and less than 0.005% of microbiome

NewGenusNR<-as.data.frame(t(NewGenusNR))
SaliniCol<-grep("Salinibacter", colnames(NewGenusNR))
NewGenusNR<-NewGenusNR[, -SaliniCol]
HumanRow<-grep("Homo", colnames(NewGenusNR))
NewGenusNR<-NewGenusNR[, -HumanRow]
MouseRow<-grep("Mus", colnames(NewGenusNR))
# NewGenusNR<-NewGenusNR[, -GeigerRow]


GeigerGenus<-NewGenusNR

# actually this is 5% cutoff
TenPercentCutoff<-floor(nrow(GeigerGenus)/20)

NonZeroCounts<-list()
for (i in 1:ncol(GeigerGenus)){
  NonZeroCounts[i]<-length(which(GeigerGenus[,i] > 0))
  
}

TenPercentNotZero<-which(NonZeroCounts >= TenPercentCutoff)
GeigerGenusNR<-GeigerGenus[,TenPercentNotZero]

# Ok let's not remove species that are missing in 10% of samples since the sample number are low
GeigerGenusNR<-as.data.frame(t(noise.removal(t(GeigerGenusNR), 0.001)))
#GeigerGenusNR<-as.data.frame(t(noise.removal(t(GeigerGenus), 0.01)))
LowSamples<-which(rowSums(GeigerGenusNR) <= 500000)
# GeigerGenusNR<-GeigerGenusNR[-LowSamples,]
minCount<-min(rowSums(GeigerGenusNR))
GeigerGenusNR<-data.frame(rrarefy(GeigerGenusNR, minCount))
GeigerGenusNR$SampleID<-row.names(GeigerGenusNR)

BackupGeigerGenus<-GeigerGenusNR
#GeigerGenusNR<-BackupGeigerGenus

names(Metadata)[1]<-"SampleID"

GeigerGenusNR<-merge(Metadata, GeigerGenusNR, by="SampleID", all.y = TRUE)
row.names(GeigerGenusNR)<-GeigerGenusNR$SampleID





LongGeigerGenusNR<-melt(GeigerGenusNR, id.vars = c("SampleID", "Recipient", "Background", "Transplant", "ImmuneSystem", "Cage", "Groups", "CageGroups"))
names(LongGeigerGenusNR)<-c("SampleID", "Recipient", "Background", "Transplant", "ImmuneSystem", "Cage", "Groups", "CageGroups", "Genus", "Count")
# LongGeigerGenusNR<-merge(LongGeigerGenusNR, CuratedGenus, by = "SubGenus", all.x = TRUE)
write.csv(LongGeigerGenusNR, file = "GeigerLongGenusTable.csv")


Genus<-as.data.frame(t(GeigerGenusNR[,c(9:ncol(GeigerGenusNR))]))

# Calculate sample diversity using 'diversity' function in Vegan
H <- data.frame(diversity(t(Genus)))
simpson <- data.frame(diversity(t(Genus), "simpson"))
shannon<-data.frame(diversity(t(Genus), "shannon"))
invsimp <- data.frame(diversity(t(Genus), "inv"))
alpha <- data.frame(fisher.alpha(t(Genus)))
## Genus richness (S) and Pielou's evenness (J):
S <- data.frame(specnumber(t(Genus)))
J <- data.frame(H/log(S))
Diversity<-cbind(simpson, shannon, invsimp,  S, J)
Diversity$SampleID<-row.names(GeigerGenusNR)
#Diversity<-Diversity[,c(1,3,5,7,8,9,10)]
names(Diversity)<-c("Simpson", "Shannon",  "InvSimpson",  "GenusNo", "Evenness", "SampleID")

pairs(cbind(simpson, shannon,  J, S), pch="+", col="blue")


DiversityGenus<-merge(Metadata, Diversity, by = "SampleID", all.y = TRUE)

col=tableau.colors


DiversityGenus_Groups<-ggplot(DiversityGenus,
                         aes(x=Groups, y=Shannon, fill=Groups)) + geom_boxplot(lwd=1,aes(color=factor(Groups),fill = NA), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Groups))) + xlab(NULL) +
  ylab("Shannon Diversity Index\n")  +  paramsAngled() + scale_color_tableau() 
ggsave(filename = "DiversityGenus_Groups.pdf", plot = DiversityGenus_Groups, width = 6,  height = 8, limitsize = FALSE)


pairwise.wilcox.test(DiversityGenus$Shannon, DiversityGenus$Groups)


Groups.mrpp<-mrpp(GeigerGenusNR[,3:length(colnames(GeigerGenusNR))], GeigerGenusNR[,2], distance = "bray")
Groups.mrpp

GroupDist<-vegdist(GeigerGenusNR[,3:length(colnames(GeigerGenusNR))], method = "bray")

GeigerGenusNR2<-as.data.frame(t(noise.removal(t(GeigerGenusNR[,9:ncol(GeigerGenusNR)]), 0.1)))

GenusLog<-log(GeigerGenusNR2[,9:ncol(GeigerGenusNR2)])
GenusCutLog<-log2(GeigerGenusNR2[,9:ncol(GeigerGenusNR2)])
GenusCutLog[GenusCutLog==-Inf]<-0

HellingerGenus<- tran(GenusCutLog, method="jaccard")
Cluster <- vegdist(HellingerGenus, method = "euclidean", diag = FALSE, upper = FALSE)

row.clus <- hclust(Cluster, "ward.D2")
cluster2<-vegdist(t(GenusCutLog), method="bray")
col.culst<-hclust(cluster2, "ward.D2")

heatmap.2(as.matrix(GenusCutLog), Rowv = as.dendrogram(row.clus), Colv=as.dendrogram(col.culst),col=metaphlan.colors(100),  margins=c(15,10),trace = "none",
          xlab = NULL, main = NULL, lhei= c(5,20))


GroupMeanDist<-meandist(GroupDist, GeigerGenusNR[,2])
plot(GroupMeanDist)
Cluster<-hclust(GroupDist)
plot(Cluster)

metadata <- GeigerGenusNR[, 1:10]
cts <- as.matrix(GeigerGenusNR[, -(1:10)])
rownames(cts) <- metadata$SampleID

cts_l2 <- glog2(cts)



row.clus <- hclust(Cluster, "ward.D2")
cluster2<-vegdist(t(GenusLog), method="jaccard")
col.culst<-hclust(GroupDist, "ward.D2")

heatmap.2(as.matrix(GenusLog), col=metaphlan.colors(100),  margins=c(15,10),trace = "none", xlab = NULL, main = NULL, lhei= c(5,30))


# col=c("dodgerblue4", "violetred4", "blue", "pink", "orange", "green", "purple")


### Young versus Old HSC : Genus #############

YoungOldHSCGenus<-subset(GeigerGenusNR, GeigerGenusNR$Transplant %in% c("Young_HSCs", "Old_HSCs"))
YoungOldHSCGenus<-YoungOldHSCGenus[, c(1,4,10:ncol(YoungOldHSCGenus))]

Wilcox<-pairwise.wilcox.test(YoungOldHSCGenus[,3], YoungOldHSCGenus[,2])
WilcoxTable<-as.data.frame(Wilcox$p.value)
names(WilcoxTable)<-names(YoungOldHSCGenus)[3]
row.names(WilcoxTable)<-"Wilcox"


for (i in 4:length(colnames(YoungOldHSCGenus))){
  x<-pairwise.wilcox.test(YoungOldHSC[,i], YoungOldHSCGenus[,2])
  y<-as.data.frame(x$p.value)
  names(y)<-names(YoungOldHSCGenus)[i]
  WilcoxTable<-cbind(WilcoxTable, y)
}

WilcoxTable<-as.data.frame(t(WilcoxTable))
names(WilcoxTable)<-c("Unadjusted_p")
WilcoxTable$FDR<-p.adjust(WilcoxTable$Unadjusted_p, method = "fdr")
WilcoxTable$Genus<-row.names(WilcoxTable)
WilcoxTable<-WilcoxTable[order(WilcoxTable$Unadjusted_p),]


GenusCounts<-as.data.frame(colSums(YoungOldHSCGenus[3:length(colnames(YoungOldHSCGenus))]))
names(GenusCounts)<-"GenusTotal"
GenusCountsBigTable<-as.data.frame(colSums(YoungOldHSCGenus[3:length(colnames(YoungOldHSCGenus))]))
names(GenusCountsBigTable)<-"GenusTotal"
TotalCountsAll<-sum(GenusCountsBigTable$GenusTotal)
GenusCounts$Genus<-row.names(GenusCounts)

for (i in 1:length(rownames(GenusCounts))){
  GenusCounts$Fraction[i]<-GenusCounts[i,1] / TotalCountsAll
}

GenusCountMeans<-aggregate(.~YoungOldHSCGenus$Transplant, mean, data=YoungOldHSCGenus[3:ncol(YoungOldHSCGenus)])
row.names(GenusCountMeans)<-GenusCountMeans$`YoungOldHSCGenus$Transplant`
GenusCountMeans<-as.data.frame(t(GenusCountMeans[,2:ncol(GenusCountMeans)]))
GenusCountMeans$Genus<-row.names(GenusCountMeans)
OverallMeans<-as.data.frame(colMeans(YoungOldHSCGenus[3:ncol(YoungOldHSCGenus)]))
OverallMeans$Genus<-row.names(OverallMeans)
names(OverallMeans)<-c("OverallMean", "Genus")
GenusCountMeans<-merge(GenusCountMeans, OverallMeans, by = "Genus", all.x = TRUE)
GenusCountMeans$Young_HSCs<-as.numeric(GenusCountMeans$Young_HSCs)
GenusCountMeans$Old_HSCs<-as.numeric(GenusCountMeans$Old_HSCs)
GenusCountMeans$FoldChange<-(GenusCountMeans$Young_HSCs-GenusCountMeans$Old_HSCs) / GenusCountMeans$OverallMean
GenusCountMeans$Log2<-log2(GenusCountMeans$Young_HSCs/GenusCountMeans$OverallMean)
GenusCountTable<-merge(GenusCountMeans, GenusCounts, by = "Genus", all.x = TRUE)
GenusCountTable<-merge(GenusCountTable, WilcoxTable, by = "Genus", all.x = TRUE)
GenusCountTable<-GenusCountTable[order(GenusCountTable$Unadjusted_p),]
GenusCountTable<-GenusCountTable[,c(1,8,4, 2,3,5,9,10)]
names(GenusCountTable)<-c("Genus", "Abundance", "OverallMean", "Young_HSCs.Mean", "Old_HSCs.Mean", "Fold.Difference", "p_Unadjusted", "FDR")
GenusCountTable$Abundance<-GenusCountTable$Abundance * 100
GenusCountTable$Fold.Difference[is.na(GenusCountTable$Fold.Difference)]<-0
GenusCountTable<-GenusCountTable[order(GenusCountTable$Fold.Difference, decreasing = TRUE),]

SigGenusCountTable<-subset(GenusCountTable, GenusCountTable$p_Unadjusted < 0.05)
SigGenusCountTable$HigherIn<-ifelse(SigGenusCountTable$Fold.Difference > 0, "Young_HSCs", "Old_HSCs")
SigGenusCountTable$HigherIn<-factor(SigGenusCountTable$HigherIn, levels = c("Young_HSCs", "Old_HSCs"))
write.csv(SigGenusCountTable, file = "SignificantGenusYoungvsOldHSC.csv")
SigYoungOldHSCGenus<-SigGenusCountTable$Genus


metadata <- YoungOldHSCGenus[, 1:2]
cts <- as.matrix(YoungOldHSCGenus[, -(1:2)])
rownames(cts) <- metadata$SampleID
grps<-metadata$Transplant

cts_l2 <- glog2(cts)

res <- compute_ef_YoungOldHSCs(d = cts_l2, g = YoungOldHSCGenus$Transplant, min_shrink = 0.3)
names(res)[1]<-"Genus"

res <- arrange(left_join(res, SigGenusCountTable), desc(abs(ef_shrunk))) %>%
  mutate(Genus = as_factor(Genus))

res<-subset(res, ! is.na(res$FDR))

res$HigherIn<-ifelse(res$ef_shrunk > 0, "Young_HSCs", "Old_HSCs")
res$HigherIn<-factor(res$HigherIn, levels = c("Young_HSCs", "Old_HSCs"))

resYoungOldHSCs<-res

# write.table(SigGenusCountTable, file = "SignificantBackgroundGenus.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)
# write.table(GenusCountTable, file = "GroupGeigerGenusNRK2.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)

col=tableau.colors[1:2]

YoungOldHSCEffect<- ggplot(filter(resYoungOldHSCs, resYoungOldHSCs$ef_shrunk > 0.3 | resYoungOldHSCs$ef_shrunk < - 0.3),
                           aes(x = reorder(Genus, ef_shrunk), y = ef_shrunk, fill = factor(HigherIn))) +
  geom_bar(colour="black", stat="identity") +
  #guides(fill=FALSE) + scale_fill_hue(l=40) +
  guides(fill=guide_legend(override.aes=list(colour=col))) +
  coord_flip() +
  labs( x = NULL) +
  ylab(expression(atop("Effect Size"))) +
  guides(fill=guide_legend(title="More abundant in:")) +
  #geom_hline(yintercept = c(-0.2, 0.2), color = "firebrick", lty = 2) +
  theme(axis.text.y = element_text(size= 12, color="black")) +
  theme(legend.title = element_text(size=12)) +
  theme(legend.text = element_text(size = 11)) 


YoungOldHSCEffect<- YoungOldHSCEffect + scale_fill_manual(values = col)

ggsave(filename = "YoungvsOldHSCL_Genus_Effect.pdf", plot = YoungOldHSCEffect, width = 12,
       height = 16, limitsize = FALSE)

# OverlappingGenus<-intersect(SigYoungOldNonTransplantedGenus, SigYoungOldHSCGenus)
DummyGenus<-GeigerGenusNR

Nitrobacter<-ggplot(DummyGenus,
                                aes(x=Groups, y=Nitrobacter, fill=NULL)) + geom_boxplot(lwd=1,aes(color=factor(ImmuneSystem),fill = NULL), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(ImmuneSystem))) + xlab(NULL) + theme_bw() +
  ylab("Nitrobacter\n")  +  paramsAngled() + scale_color_tableau() #+ scale_y_log10()
ggsave(filename = "Nitrobacter.pdf", plot = Nitrobacter, width = 10,  height = 8, limitsize = FALSE)

Alistipes<-ggplot(DummyGenus,
                    aes(x=Groups, y=Alistipes, fill=NULL)) + geom_boxplot(lwd=1,aes(color=factor(ImmuneSystem),fill = NULL), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(ImmuneSystem))) + xlab(NULL) + theme_bw() +
  ylab("Alistipes\n")  +  paramsAngled() + scale_color_tableau() #+ scale_y_log10()
ggsave(filename = "Alistipes.pdf", plot = Alistipes, width = 10,  height = 8, limitsize = FALSE)

Cardiobacterium<-ggplot(DummyGenus,
                  aes(x=Groups, y=Cardiobacterium, fill=NULL)) + geom_boxplot(lwd=1,aes(color=factor(ImmuneSystem),fill = NULL), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(ImmuneSystem))) + xlab(NULL) + theme_bw() +
  ylab("Cardiobacterium\n")  +  paramsAngled() + scale_color_tableau() #+ scale_y_log10()
ggsave(filename = "Cardiobacterium.pdf", plot = Cardiobacterium, width = 10,  height = 8, limitsize = FALSE)

Teredinibacter<-ggplot(DummyGenus,
                        aes(x=Groups, y=Teredinibacter, fill=NULL)) + geom_boxplot(lwd=1,aes(color=factor(ImmuneSystem),fill = NULL), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(ImmuneSystem))) + xlab(NULL) + theme_bw() +
  ylab("Teredinibacter\n")  +  paramsAngled() + scale_color_tableau() #+ scale_y_log10()
ggsave(filename = "Teredinibacter.pdf", plot = Teredinibacter, width = 10,  height = 8, limitsize = FALSE)

Duncaniella<-ggplot(DummyGenus,
                       aes(x=Groups, y=Duncaniella, fill=NULL)) + geom_boxplot(lwd=1,aes(color=factor(ImmuneSystem),fill = NULL), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(ImmuneSystem))) + xlab(NULL) + theme_bw() +
  ylab("Duncaniella\n")  +  paramsAngled() + scale_color_tableau() #+ scale_y_log10()
ggsave(filename = "Duncaniella.pdf", plot = Duncaniella, width = 10,  height = 8, limitsize = FALSE)

Duncaniella<-ggplot(DummyGenus,
                    aes(x=Groups, y=Duncaniella, fill=NULL)) + geom_boxplot(lwd=1,aes(color=factor(ImmuneSystem),fill = NULL), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(ImmuneSystem))) + xlab(NULL) + theme_bw() +
  ylab("Duncaniella\n")  +  paramsAngled() + scale_color_tableau() #+ scale_y_log10()
ggsave(filename = "Duncaniella.pdf", plot = Duncaniella, width = 10,  height = 8, limitsize = FALSE)

Akkermansia<-ggplot(DummyGenus,
                    aes(x=Groups, y=Akkermansia, fill=NULL)) + geom_boxplot(lwd=1,aes(color=factor(ImmuneSystem),fill = NULL), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(ImmuneSystem))) + xlab(NULL) + theme_bw() +
  ylab("Akkermansia\n")  +  paramsAngled() + scale_color_tableau() #+ scale_y_log10()
ggsave(filename = "Akkermansia.pdf", plot = Akkermansia, width = 10,  height = 8, limitsize = FALSE)

save.image(file = "GeigerData20220927")

#######################





Wilcox<-pairwise.wilcox.test(GeigerGenusNR[,3], GeigerGenusNR[,2])
WilcoxTable<-as.data.frame(Wilcox$p.value)
names(WilcoxTable)<-names(GeigerGenusNR)[3]
row.names(WilcoxTable)<-"Wilcox"


for (i in 4:length(colnames(GeigerGenusNR))){
  x<-pairwise.wilcox.test(GeigerGenusNR[,i], GeigerGenusNR[,2])
  y<-as.data.frame(x$p.value)
  names(y)<-names(GeigerGenusNR)[i]
  WilcoxTable<-cbind(WilcoxTable, y)
}

WilcoxTable<-as.data.frame(t(WilcoxTable))
names(WilcoxTable)<-c("Unadjusted_p")
WilcoxTable$FDR<-p.adjust(WilcoxTable$Unadjusted_p, method = "fdr")
WilcoxTable$Genus<-row.names(WilcoxTable)
WilcoxTable<-WilcoxTable[order(WilcoxTable$Unadjusted_p),]

GenusCounts<-as.data.frame(colSums(GeigerGenusNR[4:length(colnames(GeigerGenusNR))]))
names(GenusCounts)<-"GenusTotal"
GenusCountsBigTable<-as.data.frame(colSums(GeigerGenusNR[4:length(colnames(GeigerGenusNR))]))
names(GenusCountsBigTable)<-"GenusTotal"
TotalCountsAll<-sum(GenusCountsBigTable$GenusTotal)
GenusCounts$Genus<-row.names(GenusCounts)

for (i in 1:length(rownames(GenusCounts))){
  GenusCounts$Fraction[i]<-GenusCounts[i,1] / TotalCountsAll
}

GenusGroupMeans<-aggregate(.~GeigerGenusNR$Group, mean, data=GeigerGenusNR[3:ncol(GeigerGenusNR)])
row.names(GenusGroupMeans)<-GenusGroupMeans$`GeigerGenusNR$Group`
GenusGroupMeans<-as.data.frame(t(GenusGroupMeans[,2:ncol(GenusGroupMeans)]))
GenusGroupMeans$Genus<-row.names(GenusGroupMeans)
OverallMeans<-as.data.frame(colMeans(GeigerGenusNR[3:ncol(GeigerGenusNR)]))
OverallMeans$Genus<-row.names(OverallMeans)
names(OverallMeans)<-c("OverallMean", "Genus")
GenusGroupMeans<-merge(GenusGroupMeans, OverallMeans, by = "Genus", all.x = TRUE)
GenusGroupMeans$C57BL<-as.numeric(GenusGroupMeans$C57BL)
GenusGroupMeans$RAG1<-as.numeric(GenusGroupMeans$RAG1)
GenusGroupMeans$FoldChange<-(GenusGroupMeans$C57BL-GenusGroupMeans$RAG1) / GenusGroupMeans$OverallMean
GenusGroupMeans$Log2<-log2(GenusGroupMeans$C57BL/GenusGroupMeans$OverallMean)
GenusGroupTable<-merge(GenusGroupMeans, GenusCounts, by = "Genus", all.x = TRUE)
GenusGroupTable<-merge(GenusGroupTable, WilcoxTable, by = "Genus", all.x = TRUE)
GenusGroupTable<-GenusGroupTable[order(GenusGroupTable$Unadjusted_p),]
GenusGroupTable<-GenusGroupTable[,c(1,8,4, 2,3,5,9,10)]
names(GenusGroupTable)<-c("Genus", "Abundance", "OverallMean", "WT.Mean", "Null.Mean", "Fold.Difference", "p_Unadjusted", "FDR")
GenusGroupTable$Abundance<-GenusGroupTable$Abundance * 100
GenusGroupTable$Fold.Difference[is.na(GenusGroupTable$Fold.Difference)]<-0
GenusGroupTable<-GenusGroupTable[order(GenusGroupTable$Fold.Difference, decreasing = TRUE),]

SigGenusGroupTable<-subset(GenusGroupTable, GenusGroupTable$FDR < 0.1)

SigGenusGroupTable$HigherIn<-ifelse(SigGenusGroupTable$Fold.Difference > 0, "WT", "Null")
SigGenusGroupTable$HigherIn<-factor(SigGenusGroupTable$HigherIn, levels = c("WT", "Null"))

write.csv(SigGenusGroupTable, file = "SignificantGenusGroup.csv")

res <- compute_ef_Group(d = cts_l2, g = GeigerGenusNR$Group, min_shrink = 0.3)
names(res)[1]<-"Genus"

res <- arrange(left_join(res, SigGenusGroupTable), desc(abs(ef_shrunk))) %>%
  mutate(Genus = as_factor(Genus))

res<-subset(res, ! is.na(res$FDR))

res$HigherIn<-ifelse(res$ef_shrunk > 0, "WT", "Null")
res$HigherIn<-factor(res$HigherIn, levels = c("WT", "Null"))

resBackground<-res

write.table(SigGenusGroupTable, file = "SignificantGroupGeigerGenus.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)
# write.table(GenusGroupTable, file = "GroupGeigerGenusNRK2.csv", sep = ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)

col=tableau.colors[1:2]


Background_Effect<- ggplot(filter(resBackground, resBackground$ef_shrunk > 2 | resBackground$ef_shrunk < -1),

                      aes(x = reorder(Genus, ef_shrunk), y = ef_shrunk, fill = factor(HigherIn))) +
  geom_bar(colour="black", stat="identity") +
  #guides(fill=FALSE) + scale_fill_hue(l=40) +
  guides(fill=guide_legend(override.aes=list(colour=col))) +
  coord_flip() +
  labs( x = NULL) +
  ylab(expression(atop("Effect Size"))) +
  guides(fill=guide_legend(title="More abundant in:")) +
  #geom_hline(yintercept = c(-0.2, 0.2), color = "firebrick", lty = 2) +
  theme(axis.text.y = element_text(size= 12, color="black")) +
  theme(legend.title = element_text(size=12)) +
  theme(legend.text = element_text(size = 11)) #+ scale_color_tableau() #+
#labs( caption = "*Unadjusted p-value < 0.05 and SDA Effect Size > 0.2")

Background_Effect<- Background_Effect + scale_fill_manual(values = col)



ggsave(filename = "GroupGeigerGenusNRK2_Effect.pdf", plot = Background_Effect, width = 10,
       height = 12, limitsize = FALSE)


Bacteroides_Groups<-ggplot(GeigerGenusNR,
                           aes(x=Group, y=Bacteroides, fill=Group)) + geom_boxplot(lwd=1,aes(color=factor(Group),fill = NA), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Group))) + xlab(NULL) +
  ylab("Bacteroides\n")  +  paramsBox() + scale_color_tableau() 
ggsave(filename = "Bacteroides_Groups.pdf", plot = Bacteroides_Groups, width = 5,  height = 8, limitsize = FALSE)


Bacteroides_Ovatus<-ggplot(GeigerSpeciesNR,
                              aes(x=Group, y=Bacteroides.ovatus, fill=Group)) + geom_boxplot(lwd=1,aes(color=factor(Group),fill = NA), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Group))) + xlab(NULL) +
  ylab("Bacteroides ovatus\n")  +  paramsBox() + scale_color_tableau() 
ggsave(filename = "Bacteroides_Ovatus.pdf", plot = Bacteroides_Ovatus, width = 5,  height = 8, limitsize = FALSE)


metadata <- GeigerGenusNR[, 1:2]
cts <- as.matrix(GeigerGenusNR[, -(1:2)])
rownames(cts) <- metadata$SampleID

cts_l2 <- glog2(cts)
grps <- GeigerGenusNR$Group

newPCA<-PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "text",
                      pointsize = 2,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Genus Abundance Colored by Group Group",
                      col.ind = grps,
                      addEllipses = TRUE,
                      ellipse.alpha = 0.01,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Group Group",
                      legend.size = 11,
                      mean.point = FALSE,
                      palette = col,
                      axes.linetype = "blank"
)
GroupPCAGenus<- pcaPlot 

pdf("GroupPCAGenus.pdf")
print(GroupPCAGenus)
dev.off()


# # Which species are different in both Group and Transplant group

SigGroupSpecies<-read.csv("SignificantGroupGeigerSpeciesNRK2Copy.csv", header = TRUE, stringsAsFactors = FALSE)
SigGroupSpeciesList<-SigGroupSpecies$Species
SigGroupSpeciesList<-SigSpeciesCountTable$Species

IntersectingSpecies<-intersect(SigGroupSpeciesList, SigGroupSpeciesList)

write.csv(IntersectingSpecies, file = "TransplantGroupInteractingSpecies.csv")

SpeciesVennDiagram<-draw.pairwise.venn(area1 = length(SigGroupSpeciesList), area2 = length(SigGroupSpeciesList),
                                       cross.area = length(IntersectingSpecies), category = c("Transplant",  "Group "), lty = rep("blank", 2),
                                       fill = c("light blue", "pink"),  cat.pos =c(0, 0),
                                       fontfamily = "Helvetica", fontface = "plain", cat.default.pos = "outer", cat.fontfamily = "Helvetica")

ggsave(filename = "SpeciesVennDiagram.pdf", plot = SpeciesVennDiagram, width = 6,
       height = 6, limitsize = FALSE)


save.image(file="GeigerData20200629")


## alright now do heatmaps of microbial species and human transcripts (the human transcripts are at RNASeqAnalysis.R)


# read in non-stratified Pathway table
# the pathway table is named "Murrison_pathabundance_unstratified.tsv"

# I copied this from pCloudDrive
# see the script "Humann3Geiger.sh" to see how the table was Pathwayrated

PathwayTable<-read.csv("GeigerFiles_pathabundance20220910-cpm_unstratified.tsv", stringsAsFactors = FALSE, header = TRUE, sep = "\t")
names(PathwayTable)[1]<-"Pathway"
colnames(PathwayTable)<-gsub(".paired_Abundance.CPM", "", colnames(PathwayTable))
colnames(PathwayTable)<-gsub(".paired_Abundance.RPKs", "", colnames(PathwayTable))
PathwayTable$Pathway<-gsub("UniRef90_", "", PathwayTable$Pathway)
PathwayTable<-PathwayTable[-c(1:2),]


# change this for each sample
# PathwayTable[,c(2:ncol(PathwayTable))]<-PathwayTable[,c(2:ncol(PathwayTable))] * 1000000
PathwayTable$Total<-rowSums(PathwayTable[,c(2:ncol(PathwayTable))])
PathwayTable<-PathwayTable[order(PathwayTable$Total, decreasing = TRUE),]
PathwayTable$Pathway<-gsub(" ", ".", PathwayTable$Pathway)
PathwayTable$Pathway<-gsub(":.", ":", PathwayTable$Pathway)
row.names(PathwayTable)<-PathwayTable$Pathway

# PathwayTableSubset<-subset(PathwayTable, PathwayTable$Total > 500) # this gets us to 2154 Pathways
# PathwayTableSubset2<-subset(PathwayTable, PathwayTable$Total > 25000) # this gets us to 327 Pathways


PathwayCounts<-PathwayTable
PathwayCounts$Total<-NULL
PathwayCounts$Pathway<-NULL

PathwayTableLog<-t(log2(PathwayCounts))
PathwayTableLog2<-PathwayTableLog
PathwayTableLog[PathwayTableLog==-Inf]<-1


data.dist<-vegdist(PathwayTableLog, method="euclidean")
fit<-hclust(data.dist, method="ward.D2")
plot(fit)

Cluster <- vegdist(t(PathwayTableLog), method = "euclidean")
row.clus <- hclust(Cluster, "ward.D2")

colors = c(seq(1,96,length=101))

my_palette <- colorRampPalette(c("red", "black", "green"))(n = 100)
heatmap.2(as.matrix(PathwayTableLog),col=my_palette, margins = c(15,8),
          breaks=colors, density.info="none", trace="none", #Colv = FALSE,
          dendrogram=c("row"), symm=F,symkey=F,symbreaks=T, scale="none")


## Ok let's try PCA with the Pathway abundance table 

GeigerPathwayTable<-PathwayTable
row.names(GeigerPathwayTable)<-GeigerPathwayTable$Pathway
GeigerPathwayTable$Pathway<-NULL
GeigerPathwayTable$Total<-NULL

GeigerPathwayTabledf<-t(GeigerPathwayTable)
GeigerPathwayTabledf<-as.data.frame(GeigerPathwayTabledf)
GeigerPathwayTabledf$SampleID<-row.names(GeigerPathwayTabledf)

GeigerPathways<-merge(Metadata, GeigerPathwayTabledf, by = "SampleID", all = TRUE)
# names(GeigerPathways)<-gsub(".x", "", names(GeigerPathways))
# GeigerPathways<-GeigerPathways[,c(1, 7, 9:ncol(GeigerPathways))]
GeigerPathways<-subset(GeigerPathways, ! is.na(GeigerPathways$Groups))
colnames(GeigerPathways)<-gsub(": ", ".", colnames(GeigerPathways), fixed = TRUE)
colnames(GeigerPathways)<-gsub(" ", ".", colnames(GeigerPathways), fixed = TRUE)
colnames(GeigerPathways)<-gsub("-", ".", colnames(GeigerPathways), fixed = TRUE)
colnames(GeigerPathways)<-gsub(",", ".", colnames(GeigerPathways), fixed = TRUE)

# GeigerPathways<-GeigerPathways[-22,]

metadata <- GeigerPathways[, 1:10]
cts <- as.matrix(GeigerPathways[, -(1:10)])
rownames(cts) <- metadata$SampleID

cts_l2 <- glog2(cts)
grps <- factor(GeigerPathways$Groups)

# col=c(rep(tableau.colors[1], 5),rep(tableau.colors[2],5), rep(tableau.colors[3], 5), rep(tableau.colors[4], 5)) 
col=tableau.colors
newPCA<-PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
pcaPlot<-fviz_pca_ind(newPCA, 
                      geom.ind = "point",
                      pointsize = 0.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Pathway Abundance Data Colored by Mouse Group",
                      col.ind = grps,
                      addEllipses = TRUE,
                      ellipse.alpha = 0.01,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Sample",
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank"
)
PathwaysPCA<- pcaPlot 


pdf("PathwaysPCA.pdf")
print(PathwaysPCA)
dev.off()

TransplantedPathwayInd<-GeigerPathways$SampleID[GeigerPathways$Transplant %in% c("Young_HSCs", "Old_HSCs")]

metadata <- GeigerPathways[, 1:10]
cts <- as.matrix(GeigerPathways[, -(1:10)])
rownames(cts) <- metadata$SampleID

cts_l2 <- glog2(cts)
grps <- factor(GeigerPathways$GroupsExp)

pcaPlot<-fviz_pca_ind(newPCA, 
                      geom.ind = "point",
                      pointsize = 0.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Pathway Abundance Data Colored by Mouse Group",
                      col.ind = grps,
                      addEllipses = TRUE,
                      ellipse.alpha = 0.01,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Sample",
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank",
                      select.ind = list(name = TransplantedPathwayInd)
)
TransplantedPathwaysPCA<- pcaPlot + scale_color_tableau()


pdf("TransplantedPathwaysPCA.pdf")
print(TransplantedPathwaysPCA)
dev.off()



UntransplantedPathwayInd<-GeigerPathways$SampleID[GeigerPathways$Transplant %in% c("None") & ! GeigerPathways$Background == "RAG1"]

metadata <- GeigerPathways[, 1:10]
cts <- as.matrix(GeigerPathways[, -(1:10)])
rownames(cts) <- metadata$SampleID

cts_l2 <- glog2(cts)
grps <- factor(GeigerPathways$GroupsExp)

pcaPlot<-fviz_pca_ind(newPCA, 
                      geom.ind = "point",
                      pointsize = 0.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Pathway Abundance Data Colored by Mouse Group",
                      col.ind = grps,
                      addEllipses = TRUE,
                      ellipse.alpha = 0.01,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Sample",
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank",
                      select.ind = list(name = UntransplantedPathwayInd)
)
UntransplantedPathwaysPCA<- pcaPlot + scale_color_tableau()


pdf("UntransplantedPathwaysPCA.pdf")
print(UntransplantedPathwaysPCA)
dev.off()

# Now Pathway differences

# first do Rag1-Young vs Rag1-Old
GeigerPathwaysYO<-subset(GeigerPathways, GeigerPathways$Groups %in% c("RAG1-Young_HSCs", "RAG1-Old_HSCs"))
GeigerPathwaysYO<-GeigerPathwaysYO[,c(1,4,11:ncol(GeigerPathwaysYO))]

Wilcox<-pairwise.wilcox.test(GeigerPathwaysYO[,3], GeigerPathwaysYO[,2])
WilcoxTable<-as.data.frame(Wilcox$p.value)
names(WilcoxTable)<-names(GeigerPathwaysYO)[3]
row.names(WilcoxTable)<-"Wilcox"


for (i in 4:length(colnames(GeigerPathwaysYO))){
  x<-pairwise.wilcox.test(GeigerPathwaysYO[,i], GeigerPathwaysYO[,2])
  y<-as.data.frame(x$p.value)
  names(y)<-names(GeigerPathwaysYO)[i]
  WilcoxTable<-cbind(WilcoxTable, y)
}

WilcoxTable<-as.data.frame(t(WilcoxTable))
names(WilcoxTable)<-c("Unadjusted_p")
WilcoxTable$FDR<-p.adjust(WilcoxTable$Unadjusted_p, method = "fdr")
WilcoxTable$Species<-row.names(WilcoxTable)
WilcoxTable<-WilcoxTable[order(WilcoxTable$Unadjusted_p),]
names(WilcoxTable)<-c("Unadjusted_p", "FDR", "Pathway")

write.csv(WilcoxTable, file = "GeigerMicrobialPathwayWilcoxTable.csv")
sigPathwaysYORag<-subset(WilcoxTable, WilcoxTable$Unadjusted_p < 0.05)
write.csv(sigPathwaysYORag, file = "SignificantGeigerPathwaysYO.csv")

GeigerRAG_YvsOHSC<-sigPathwaysYORag

GeigerPathwaysYO.mrpp<-mrpp(GeigerPathwaysYO[,3:length(colnames(GeigerPathwaysYO))], GeigerPathwaysYO[,2], distance = "euclidean", na.rm = TRUE)
GeigerPathwaysYO.mrpp

# Now Rag vs C57BL


GeigerPathwaysRC<-subset(GeigerPathways, GeigerPathways$Groups %in% c("Young_RAG1-Untransplanted", "Young_C57BL-Untransplanted"))
GeigerPathwaysRC<-GeigerPathwaysRC[,c(1,3,11:ncol(GeigerPathwaysRC))]


Wilcox<-pairwise.wilcox.test(GeigerPathwaysRC[,3], GeigerPathwaysRC[,2])
WilcoxTable<-as.data.frame(Wilcox$p.value)
names(WilcoxTable)<-names(GeigerPathwaysRC)[3]
row.names(WilcoxTable)<-"Wilcox"


for (i in 4:length(colnames(GeigerPathwaysRC))){
  x<-pairwise.wilcox.test(GeigerPathwaysRC[,i], GeigerPathwaysRC[,2])
  y<-as.data.frame(x$p.value)
  names(y)<-names(GeigerPathwaysRC)[i]
  WilcoxTable<-cbind(WilcoxTable, y)
}

WilcoxTable<-as.data.frame(t(WilcoxTable))
names(WilcoxTable)<-c("Unadjusted_p")
WilcoxTable$FDR<-p.adjust(WilcoxTable$Unadjusted_p, method = "fdr")
WilcoxTable$Species<-row.names(WilcoxTable)
WilcoxTable<-WilcoxTable[order(WilcoxTable$Unadjusted_p),]
names(WilcoxTable)<-c("Unadjusted_p", "FDR", "Pathway")

write.csv(WilcoxTable, file = "GeigerPathwayWilcoxTableRagC57.csv")
sigPathways<-subset(WilcoxTable, WilcoxTable$Unadjusted_p < 0.05)
write.csv(sigPathways, file = "SignificantGeigerPathwaysRagC57.csv")

GeigerPathwaysRC.mrpp<-mrpp(GeigerPathwaysRC[,3:length(colnames(GeigerPathwaysRC))], GeigerPathwaysRC[,2], distance = "euclidean")
GeigerPathwaysRC.mrpp #  p = 0.025

# make a heatmap of logratio:

# let's calculate fold change and also significance of difference:
## this is messed up because I changed the columns included. Check later.

save.image(file = "GeigerData20221006")

GeigerNewPathwayTable<-as.data.frame(t(aggregate( . ~ Groups, data =GeigerPathwaysYO[,2:ncol(GeigerPathwaysYO)], FUN = mean)))
GeigerNewPathwayTable<-GeigerNewPathwayTable[-1,]
GeigerNewPathwayTable<-as.data.frame(data.matrix(GeigerNewPathwayTable))
colnames(GeigerNewPathwayTable)<-c("RAG1_Young_HSCs", "RAG1_Old_HSCs")

GeigerNewPathwayTable$Difference<-(GeigerNewPathwayTable$RAG1_Young_HSCs -GeigerNewPathwayTable$RAG1_Old_HSCs)
GeigerNewPathwayTable$FC<-foldchange(GeigerNewPathwayTable$RAG1_Young_HSCs, GeigerNewPathwayTable$RAG1_Old_HSCs)
GeigerNewPathwayTable$logratio<-foldchange2logratio(GeigerNewPathwayTable$FC)
# GeigerNewPathwayTable$log2Delta<-ifelse(sign(GeigerNewPathwayTable$Difference) == -1, -GeigerNewPathwayTable$log2Delta,
#                                       ifelse(GeigerNewPathwayTable$Difference == 0, 0, GeigerNewPathwayTable$log2Delta))
GeigerNewPathwayTable<-GeigerNewPathwayTable[order(GeigerNewPathwayTable$logratio, decreasing = FALSE),]
GeigerNewPathwayTable$Pathway<-row.names(GeigerNewPathwayTable)

# Ok that's cool. Now figure out which are statistically signficant from the GeigerNewGeneNR table
Kruskal<-kruskal.test(GeigerPathwaysYO[,3], factor(GeigerPathwaysYO[,2]))
KruskalTable<-as.data.frame(Kruskal$p.value)
names(KruskalTable)<-names(GeigerPathwaysYO)[3]
row.names(KruskalTable)<-"Kruskal"

for (i in 4:length(colnames(GeigerPathwaysYO))){
  x<-kruskal.test(GeigerPathwaysYO[,i], factor(GeigerPathwaysYO[,2]))
  y<-as.data.frame(x$p.value)
  names(y)<-names(GeigerPathwaysYO)[i]
  KruskalTable<-cbind(KruskalTable, y)
}

KruskalTable<-as.data.frame(t(KruskalTable))
names(KruskalTable)<-c("Unadjusted_p")
KruskalTable$FDR<-p.adjust(KruskalTable$Unadjusted_p, method = "fdr")
KruskalTable$Pathway<-row.names(KruskalTable)
KruskalTable<-KruskalTable[order(KruskalTable$Unadjusted_p),]

GeigerNewPathwayTable<-merge(GeigerNewPathwayTable, KruskalTable, by = "Pathway", all.x = TRUE)

GeigerNewPathwayTable2<-subset(GeigerNewPathwayTable, GeigerNewPathwayTable$Unadjusted_p < 0.05)
GeigerNewPathwayTable2<-GeigerNewPathwayTable2[order(GeigerNewPathwayTable2$logratio, decreasing = FALSE),]
GeigerNewPathwayTable3<-subset(GeigerNewPathwayTable, GeigerNewPathwayTable$FDR < 0.2)
GeigerNewPathwayTable<-GeigerNewPathwayTable[order(GeigerNewPathwayTable$logratio, decreasing = FALSE),]

GeigerNewPathwayTable$Pathway<-gsub("\t", " ", GeigerNewPathwayTable$Pathway)
GeigerNewPathwayTable$Pathway<-gsub("   ", " ", GeigerNewPathwayTable$Pathway)
GeigerNewPathwayTable$Pathway<-gsub("  ", " ", GeigerNewPathwayTable$Pathway)
GeigerNewPathwayTable$Pathway<-gsub(",", "", GeigerNewPathwayTable$Pathway)

SigRAGYOPathway<-subset(GeigerNewPathwayTable, GeigerNewPathwayTable$Unadjusted_p < 0.05)

matrixtable<-GeigerNewPathwayGeigerNewPathwayTableTable2[,6:7]
row.names(matrixtable)<-GeigerNewPathwayTable2$Pathway
matrixtable<-as.data.frame(matrixtable)
matrixtable<-subset(matrixtable, matrixtable$logratio < -1 | matrixtable$logratio > 1)
matrixtable<-as.matrix(matrixtable)
heatmap.2(matrixtable, Rowv = NULL, Colv = FALSE, col = diverging.colors(100), trace = "none",
          tracecol= "black",hline = NULL, vline = NULL, margins=c(2,25),
          xlab = NULL, main = NULL)


heatmap.2(t(matrixtable), Rowv = NULL, Colv = FALSE, col = diverging.colors(100), trace = "none",
          tracecol= "black",hline = NULL, vline = NULL, margins=c(15,2),
          xlab = NULL, main = NULL)



# Redo for untransplanted young versus old:

GeigerRCPathwayTable<-as.data.frame(t(aggregate( . ~ Groups, data =GeigerPathwaysRC[,2:ncol(GeigerPathwaysRC)], FUN = mean)))
GeigerRCPathwayTable<-GeigerRCPathwayTable[-1,]
GeigerRCPathwayTable<-as.data.frame(data.matrix(GeigerRCPathwayTable))
colnames(GeigerRCPathwayTable)<-c("Young_C57BL_Untransplanted", "Young_RAG1_Untransplanted")

GeigerRCPathwayTable$Difference<-(GeigerRCPathwayTable$Young_C57BL_Untransplanted -GeigerRCPathwayTable$Young_RAG1_Untransplanted)
GeigerRCPathwayTable$FC<-foldchange(GeigerRCPathwayTable$Young_C57BL_Untransplanted, GeigerRCPathwayTable$Young_RAG1_Untransplanted)
GeigerRCPathwayTable$logratio<-foldchange2logratio(GeigerRCPathwayTable$FC)
# GeigerRCPathwayTable$log2Delta<-ifelse(sign(GeigerRCPathwayTable$Difference) == -1, -GeigerRCPathwayTable$log2Delta,
#                                       ifelse(GeigerRCPathwayTable$Difference == 0, 0, GeigerRCPathwayTable$log2Delta))
GeigerRCPathwayTable<-GeigerRCPathwayTable[order(GeigerRCPathwayTable$logratio, decreasing = FALSE),]
GeigerRCPathwayTable$Pathway<-row.names(GeigerRCPathwayTable)

# Ok that's cool. Now figure out which are statistically signficant from the GeigerNewGeneNR table
Kruskal<-kruskal.test(GeigerPathwaysRC[,3], factor(GeigerPathwaysRC[,2]))
KruskalTable<-as.data.frame(Kruskal$p.value)
names(KruskalTable)<-names(GeigerPathwaysRC)[3]
row.names(KruskalTable)<-"Kruskal"

for (i in 4:length(colnames(GeigerPathwaysRC))){
  x<-kruskal.test(GeigerPathwaysRC[,i], factor(GeigerPathwaysRC[,2]))
  y<-as.data.frame(x$p.value)
  names(y)<-names(GeigerPathwaysRC)[i]
  KruskalTable<-cbind(KruskalTable, y)
}

KruskalTable<-as.data.frame(t(KruskalTable))
names(KruskalTable)<-c("Unadjusted_p")
KruskalTable$FDR<-p.adjust(KruskalTable$Unadjusted_p, method = "fdr")
KruskalTable$Pathway<-row.names(KruskalTable)
KruskalTable<-KruskalTable[order(KruskalTable$Unadjusted_p),]

GeigerRCPathwayTable<-merge(GeigerRCPathwayTable, KruskalTable, by = "Pathway", all.x = TRUE)

GeigerRCPathwayTable2<-subset(GeigerRCPathwayTable, GeigerRCPathwayTable$Unadjusted_p < 0.05)
GeigerRCPathwayTable2<-GeigerRCPathwayTable2[order(GeigerRCPathwayTable2$logratio, decreasing = FALSE),]
GeigerRCPathwayTable3<-subset(GeigerRCPathwayTable, GeigerRCPathwayTable$FDR < 0.1)
GeigerRCPathwayTable<-GeigerRCPathwayTable[order(GeigerRCPathwayTable$logratio, decreasing = FALSE),]

GeigerRCPathwayTable$Pathway<-gsub("\t", " ", GeigerRCPathwayTable$Pathway)
GeigerRCPathwayTable$Pathway<-gsub("   ", " ", GeigerRCPathwayTable$Pathway)
GeigerRCPathwayTable$Pathway<-gsub("  ", " ", GeigerRCPathwayTable$Pathway)
GeigerRCPathwayTable$Pathway<-gsub(",", "", GeigerRCPathwayTable$Pathway)

GeigerRCPathwayTable3$Pathway<-gsub("\t", " ", GeigerRCPathwayTable3$Pathway)
GeigerRCPathwayTable3$Pathway<-gsub("   ", " ", GeigerRCPathwayTable3$Pathway)
GeigerRCPathwayTable3$Pathway<-gsub("  ", " ", GeigerRCPathwayTable3$Pathway)
GeigerRCPathwayTable3$Pathway<-gsub(",", "", GeigerRCPathwayTable3$Pathway)



matrixtable<-GeigerRCPathwayTable3[,6:7]
row.names(matrixtable)<-GeigerRCPathwayTable3$Pathway
matrixtable<-as.data.frame(matrixtable)
matrixtable<-subset(matrixtable, matrixtable$logratio < -1 | matrixtable$logratio > 1)
matrixtable<-as.matrix(matrixtable)
heatmap.2(matrixtable, Rowv = NULL, Colv = FALSE, col = diverging.colors(100), trace = "none",
          tracecol= "black",hline = NULL, vline = NULL, margins=c(2,25),
          xlab = NULL, main = NULL)


heatmap.2(t(matrixtable), Rowv = NULL, Colv = FALSE, col = diverging.colors(100), trace = "none",
          tracecol= "black",hline = NULL, vline = NULL, margins=c(15,2),
          xlab = NULL, main = NULL)

# Young vs old C57BL

GeigerPathwaysYOWT<-subset(GeigerPathways, GeigerPathways$Groups %in% c ("Young_C57BL-Untransplanted", "Old_C57BL-Untransplanted"))

GeigerPathwayWTYO.mrpp<-mrpp(GeigerPathwaysYOWT[,3:length(colnames(GeigerPathwaysYOWT))], GeigerPathwaysYOWT[,2], distance = "euclidean")
GeigerPathwayWTYO.mrpp



GeigerPathwayTableWTYO<-as.data.frame(t(aggregate( . ~ Groups, data =GeigerPathwaysYOWT[,2:ncol(GeigerPathwaysYOWT)], FUN = mean)))
GeigerPathwayTableWTYO<-GeigerPathwayTableWTYO[-1,]
GeigerPathwayTableWTYO<-as.data.frame(data.matrix(GeigerPathwayTableWTYO))
colnames(GeigerPathwayTableWTYO)<-c("Young_C57BL_Untransplanted", "Old_C57BL_Untransplanted")

GeigerPathwayTableWTYO$Difference<-(GeigerPathwayTableWTYO$Young_C57BL_Untransplanted -GeigerPathwayTableWTYO$Old_C57BL_Untransplanted)
GeigerPathwayTableWTYO$FC<-foldchange(GeigerPathwayTableWTYO$Young_C57BL_Untransplanted, GeigerPathwayTableWTYO$Old_C57BL_Untransplanted)
GeigerPathwayTableWTYO$logratio<-foldchange2logratio(GeigerPathwayTableWTYO$FC)
# GeigerPathwayTableWTYO$log2Delta<-ifelse(sign(GeigerPathwayTableWTYO$Difference) == -1, -GeigerPathwayTableWTYO$log2Delta,
#                                       ifelse(GeigerPathwayTableWTYO$Difference == 0, 0, GeigerPathwayTableWTYO$log2Delta))
GeigerPathwayTableWTYO<-GeigerPathwayTableWTYO[order(GeigerPathwayTableWTYO$logratio, decreasing = FALSE),]
GeigerPathwayTableWTYO$Pathway<-row.names(GeigerPathwayTableWTYO)

# Ok that's cool. Now figure out which are statistically signficant from the GeigerNewGeneNR table
Kruskal<-kruskal.test(GeigerPathwaysYOWT[,3], factor(GeigerPathwaysYOWT[,2]))
KruskalTable<-as.data.frame(Kruskal$p.value)
names(KruskalTable)<-names(GeigerPathwaysYOWT)[3]
row.names(KruskalTable)<-"Kruskal"

for (i in 4:length(colnames(GeigerPathwaysYOWT))){
  x<-kruskal.test(GeigerPathwaysYOWT[,i], factor(GeigerPathwaysYOWT[,2]))
  y<-as.data.frame(x$p.value)
  names(y)<-names(GeigerPathwaysYOWT)[i]
  KruskalTable<-cbind(KruskalTable, y)
}

KruskalTable<-as.data.frame(t(KruskalTable))
names(KruskalTable)<-c("Unadjusted_p")
KruskalTable$FDR<-p.adjust(KruskalTable$Unadjusted_p, method = "fdr")
KruskalTable$Pathway<-row.names(KruskalTable)
KruskalTable<-KruskalTable[order(KruskalTable$Unadjusted_p),]

GeigerPathwayTableWTYO<-merge(GeigerPathwayTableWTYO, KruskalTable, by = "Pathway", all.x = TRUE)

GeigerPathwayTableWTYO2<-subset(GeigerPathwayTableWTYO, GeigerPathwayTableWTYO$Unadjusted_p < 0.05)
GeigerPathwayTableWTYO2<-GeigerPathwayTableWTYO2[order(GeigerPathwayTableWTYO2$logratio, decreasing = FALSE),]
GeigerPathwayTableWTYO3<-subset(GeigerPathwayTableWTYO, GeigerPathwayTableWTYO$FDR < 0.2)
GeigerPathwayTableWTYO<-GeigerPathwayTableWTYO[order(GeigerPathwayTableWTYO$logratio, decreasing = FALSE),]

GeigerPathwayTableWTYO2$Pathway<-gsub("\t", " ", GeigerPathwayTableWTYO2$Pathway)
GeigerPathwayTableWTYO2$Pathway<-gsub("   ", " ", GeigerPathwayTableWTYO2$Pathway)
GeigerPathwayTableWTYO2$Pathway<-gsub("  ", " ", GeigerPathwayTableWTYO2$Pathway)
GeigerPathwayTableWTYO2$Pathway<-gsub(",", "", GeigerPathwayTableWTYO2$Pathway)



matrixtable<-GeigerPathwayTableWTYO2[,6:7]
row.names(matrixtable)<-GeigerPathwayTableWTYO2$Pathway
matrixtable<-as.data.frame(matrixtable)
matrixtable<-subset(matrixtable, matrixtable$logratio < -1 | matrixtable$logratio > 1)
matrixtable<-as.matrix(matrixtable)
heatmap.2(matrixtable, Rowv = NULL, Colv = FALSE, col = diverging.colors(100), trace = "none",
          tracecol= "black",hline = NULL, vline = NULL, margins=c(2,25),
          xlab = NULL, main = NULL)


heatmap.2(t(matrixtable), Rowv = NULL, Colv = FALSE, col = diverging.colors(100), trace = "none",
          tracecol= "black",hline = NULL, vline = NULL, margins=c(15,2),
          xlab = NULL, main = NULL)

## 20211018: Look at overlapping genes/pathways between aged and young C57 and RAG 
## do intersecting graphs

OverlappingYOPathways<-intersect(GeigerPathwayTableWTYO2$Pathway, SigRAGYOPathway$Pathway)


SpeciesVennDiagram<-draw.pairwise.venn(area1 = length(GeigerPathwayTableWTYO2$Pathway), area2 = length(SigRAGYOPathway$Pathway),
                                       cross.area = length(OverlappingYOPathways), category = c("Young.Old.Mice",  "Young.Old.HSC"), lty = rep("blank", 2),
                                       fill = c("light blue", "pink"),  cat.pos =c(0, 0),
                                       fontfamily = "Helvetica", fontface = "plain", cat.default.pos = "outer", cat.fontfamily = "Helvetica")

ggsave(filename = "PathwayVennDiagram.pdf", plot = SpeciesVennDiagram, width = 6,
       height = 6, limitsize = FALSE)

YOMicePathwayList<-subset(GeigerPathwayTableWTYO2, GeigerPathwayTableWTYO2$Pathway %in% OverlappingYOPathways)
YORAGPathwayList<-subset(SigRAGYOPathway, SigRAGYOPathway$Pathway  %in% OverlappingYOPathways)

ComparisonTable<-merge(YOMicePathwayList, YORAGPathwayList, by = "Pathway", all = TRUE)
write.csv(ComparisonTable, file = "IntersectingPathways.csv")


COBALSYN.PWY:adenosylcobalamin.salvage.from.cobinamide.I

col<-tableau.colors
COBALSYN<-ggplot(GeigerPathways,
                           aes(x=Groups, y=`COBALSYN.PWY:adenosylcobalamin.salvage.from.cobinamide.I`, fill=Groups)) + geom_boxplot(lwd=1,aes(color=factor(Groups),fill = NA), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Groups))) + xlab(NULL) +
  ylab("COBALSYN.PWY\n")  +  paramsAngled() + scale_color_tableau() 
ggsave(filename = "COBALSYN.pdf", plot = COBALSYN, width = 5,  height = 8, limitsize = FALSE)


save.image(file="GeigerData20220927")


### Gene families


GeneTable<-read.csv("GeigerFiles_genefamilies20211018-cpm_unstratified.tsv", stringsAsFactors = FALSE, header = TRUE, sep = "\t")
names(GeneTable)[1]<-"Gene"
colnames(GeneTable)<-gsub(".paired_Abundance.CPM", "", colnames(GeneTable))
colnames(GeneTable)<-gsub(".paired_Abundance.RPKs", "", colnames(GeneTable))
GeneTable$Gene<-gsub("UniRef90_", "", GeneTable$Gene)
GeneTable<-GeneTable[-c(1:2),]


# change this for each sample
# GeneTable[,c(2:ncol(GeneTable))]<-GeneTable[,c(2:ncol(GeneTable))] * 1000000
GeneTable$Total<-rowSums(GeneTable[,c(2:ncol(GeneTable))])
GeneTable<-GeneTable[order(GeneTable$Total, decreasing = TRUE),]
GeneTable$Gene<-gsub(" ", ".", GeneTable$Gene)
GeneTable$Gene<-gsub(":.", ":", GeneTable$Gene)
row.names(GeneTable)<-GeneTable$Gene

GeneTableSubset<-subset(GeneTable, GeneTable$Total > 500) # this gets us to 5355 Genes
GeneTableSubset2<-subset(GeneTable, GeneTable$Total > 1000) # this gets us to 966 Genes


GeneCounts<-GeneTableSubset
GeneCounts$Total<-NULL
GeneCounts$Gene<-NULL

GeneTableLog<-t(log2(GeneCounts))
GeneTableLog2<-GeneTableLog
GeneTableLog[GeneTableLog==-Inf]<-1


data.dist<-vegdist(GeneTableLog, method="euclidean")
fit<-hclust(data.dist, method="ward.D2")
plot(fit)

Cluster <- vegdist(t(GeneTableLog), method = "euclidean")
row.clus <- hclust(Cluster, "ward.D2")

colors = c(seq(1,96,length=101))

my_palette <- colorRampPalette(c("red", "black", "green"))(n = 100)
heatmap.2(as.matrix(GeneTableLog),col=my_palette, margins = c(15,8),
          breaks=colors, density.info="none", trace="none", #Colv = FALSE,
          dendrogram=c("row"), symm=F,symkey=F,symbreaks=T, scale="none")


## Ok let's try PCA with the Gene abundance table 

GeigerGeneTable<-GeneTableSubset
row.names(GeigerGeneTable)<-GeigerGeneTable$Gene
GeigerGeneTable$Gene<-NULL
GeigerGeneTable$Total<-NULL

GeigerGeneTable<-t(GeigerGeneTable)
GeigerGeneTabledf<-as.data.frame(GeigerGeneTable)
GeigerGeneTabledf$SampleID<-row.names(GeigerGeneTabledf)

GeigerGenes<-merge(Metadata, GeigerGeneTabledf, by = "SampleID", all.y = TRUE)
GeigerGenes<-GeigerGenes[,c(1, 7, 9:ncol(GeigerGenes))]
GeigerGenes<-subset(GeigerGenes, ! is.na(GeigerGenes$Groups))
colnames(GeigerGenes)<-gsub(": ", ".", colnames(GeigerGenes), fixed = TRUE)
colnames(GeigerGenes)<-gsub(" ", ".", colnames(GeigerGenes), fixed = TRUE)
colnames(GeigerGenes)<-gsub("-", ".", colnames(GeigerGenes), fixed = TRUE)
colnames(GeigerGenes)<-gsub(",", ".", colnames(GeigerGenes), fixed = TRUE)

# GeigerGenes<-GeigerGenes[-22,]

metadata <- GeigerGenes[, 1:2]
cts <- as.matrix(GeigerGenes[, -(1:2)])
rownames(cts) <- metadata$SampleID

cts_l2 <- glog2(cts)
grps <- factor(GeigerGenes$Groups)

# col=c(rep(tableau.colors[1], 5),rep(tableau.colors[2],5), rep(tableau.colors[3], 5), rep(tableau.colors[4], 5)) 
col=tableau.colors
newPCA<-PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
pcaPlot<-fviz_pca_ind(newPCA, 
                      geom.ind = "point",
                      pointsize = 0.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Gene Abundance Data Colored by Mouse Group",
                      col.ind = grps,
                      addEllipses = TRUE,
                      ellipse.alpha = 0.01,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Sample",
                      legend.size = 11,
                      mean.point = TRUE,
                      palette = col,
                      axes.linetype = "blank"
)
GenesPCA<- pcaPlot 


pdf("GenesPCA.pdf")
print(GenesPCA)
dev.off()


# look at significant genes and make heatmap

GeigerGenesRC<-subset(GeigerGenes, GeigerGenes$Groups %in% c("Young_C57BL-Untransplanted", "Young_RAG1-Untransplanted"))

GeigerRCGeneTable<-as.data.frame(t(aggregate( . ~ Groups, data =GeigerGenesRC[,2:ncol(GeigerGenesRC)], FUN = mean)))
GeigerRCGeneTable<-GeigerRCGeneTable[-1,]
GeigerRCGeneTable<-as.data.frame(data.matrix(GeigerRCGeneTable))
colnames(GeigerRCGeneTable)<-c("Young_C57BL_Untransplanted", "Young_RAG1_Untransplanted")

GeigerRCGeneTable$Difference<-(GeigerRCGeneTable$Young_C57BL_Untransplanted -GeigerRCGeneTable$Young_RAG1_Untransplanted)
GeigerRCGeneTable$FC<-foldchange(GeigerRCGeneTable$Young_C57BL_Untransplanted, GeigerRCGeneTable$Young_RAG1_Untransplanted)
GeigerRCGeneTable$logratio<-foldchange2logratio(GeigerRCGeneTable$FC)
# GeigerRCGeneTable$log2Delta<-ifelse(sign(GeigerRCGeneTable$Difference) == -1, -GeigerRCGeneTable$log2Delta,
#                                       ifelse(GeigerRCGeneTable$Difference == 0, 0, GeigerRCGeneTable$log2Delta))
GeigerRCGeneTable<-GeigerRCGeneTable[order(GeigerRCGeneTable$logratio, decreasing = FALSE),]
GeigerRCGeneTable$Gene<-row.names(GeigerRCGeneTable)

# Ok that's cool. Now figure out which are statistically signficant from the GeigerNewGeneNR table
Kruskal<-kruskal.test(GeigerGenesRC[,3], factor(GeigerGenesRC[,2]))
KruskalTable<-as.data.frame(Kruskal$p.value)
names(KruskalTable)<-names(GeigerGenesRC)[3]
row.names(KruskalTable)<-"Kruskal"

for (i in 4:length(colnames(GeigerGenesRC))){
  x<-kruskal.test(GeigerGenesRC[,i], factor(GeigerGenesRC[,2]))
  y<-as.data.frame(x$p.value)
  names(y)<-names(GeigerGenesRC)[i]
  KruskalTable<-cbind(KruskalTable, y)
}

KruskalTable<-as.data.frame(t(KruskalTable))
names(KruskalTable)<-c("Unadjusted_p")
KruskalTable$FDR<-p.adjust(KruskalTable$Unadjusted_p, method = "fdr")
KruskalTable$Gene<-row.names(KruskalTable)
KruskalTable<-KruskalTable[order(KruskalTable$Unadjusted_p),]

GeigerRCGeneTable<-merge(GeigerRCGeneTable, KruskalTable, by = "Gene", all.x = TRUE)

GeigerRCGeneTable2<-subset(GeigerRCGeneTable, GeigerRCGeneTable$Unadjusted_p < 0.05)
GeigerRCGeneTable2<-GeigerRCGeneTable2[order(GeigerRCGeneTable2$logratio, decreasing = FALSE),]
GeigerRCGeneTable3<-subset(GeigerRCGeneTable, GeigerRCGeneTable$FDR < 0.01)
GeigerRCGeneTable<-GeigerRCGeneTable[order(GeigerRCGeneTable$logratio, decreasing = FALSE),]

GeigerRCGeneTable$Gene<-gsub("\t", " ", GeigerRCGeneTable$Gene)
GeigerRCGeneTable$Gene<-gsub("   ", " ", GeigerRCGeneTable$Gene)
GeigerRCGeneTable$Gene<-gsub("  ", " ", GeigerRCGeneTable$Gene)
GeigerRCGeneTable$Gene<-gsub(",", "", GeigerRCGeneTable$Gene)

GeigerRCGeneTable3$Gene<-gsub("\t", " ", GeigerRCGeneTable3$Gene)
GeigerRCGeneTable3$Gene<-gsub("   ", " ", GeigerRCGeneTable3$Gene)
GeigerRCGeneTable3$Gene<-gsub("  ", " ", GeigerRCGeneTable3$Gene)
GeigerRCGeneTable3$Gene<-gsub(",", "", GeigerRCGeneTable3$Gene)



matrixtable<-GeigerRCGeneTable3[,6:7]
row.names(matrixtable)<-GeigerRCGeneTable3$Gene
matrixtable<-as.data.frame(matrixtable)
matrixtable<-subset(matrixtable, matrixtable$logratio < -1 | matrixtable$logratio > 1)
matrixtable<-as.matrix(matrixtable)
heatmap.2(matrixtable, Rowv = NULL, Colv = FALSE, col = diverging.colors(100), trace = "none",
          tracecol= "black",hline = NULL, vline = NULL, margins=c(2,25),
          xlab = NULL, main = NULL)


heatmap.2(t(matrixtable), Rowv = NULL, Colv = FALSE, col = diverging.colors(100), trace = "none",
          tracecol= "black",hline = NULL, vline = NULL, margins=c(15,2),
          xlab = NULL, main = NULL)

GeigerGenesRC.mrpp<-mrpp(GeigerGenesRC[,3:length(colnames(GeigerGenesRC))], GeigerGenesRC[,2], distance = "euclidean")
GeigerGenesRC.mrpp

GeigerGenesYO.mrpp<-mrpp(GeigerGenesYO[,3:length(colnames(GeigerGenesRC))], GeigerGenesYO[,2], distance = "euclidean")
GeigerGenesYO.mrpp




### Genes different between YoungHSC and OldHSC

GeigerGenesYO<-subset(GeigerGenes, GeigerGenes$Groups %in% c("RAG1-Young_HSCs", "RAG1-Old_HSCs"))

GeigerYOGeneTable<-as.data.frame(t(aggregate( . ~ Groups, data =GeigerGenesYO[,2:ncol(GeigerGenesYO)], FUN = mean)))
GeigerYOGeneTable<-GeigerYOGeneTable[-1,]
GeigerYOGeneTable<-as.data.frame(data.matrix(GeigerYOGeneTable))
colnames(GeigerYOGeneTable)<-c("RAG1_Young_HSCs", "RAG1_Old_HSCs")

GeigerYOGeneTable$Difference<-(GeigerYOGeneTable$RAG1_Young_HSCs -GeigerYOGeneTable$RAG1_Old_HSCs)
GeigerYOGeneTable$FC<-foldchange(GeigerYOGeneTable$RAG1_Young_HSCs, GeigerYOGeneTable$RAG1_Old_HSCs)
GeigerYOGeneTable$logratio<-foldchange2logratio(GeigerYOGeneTable$FC)
# GeigerYOGeneTable$log2Delta<-ifelse(sign(GeigerYOGeneTable$Difference) == -1, -GeigerYOGeneTable$log2Delta,
#                                       ifelse(GeigerYOGeneTable$Difference == 0, 0, GeigerYOGeneTable$log2Delta))
GeigerYOGeneTable<-GeigerYOGeneTable[order(GeigerYOGeneTable$logratio, decreasing = FALSE),]
GeigerYOGeneTable$Gene<-row.names(GeigerYOGeneTable)

# Ok that's cool. Now figure out which are statistically signficant from the GeigerNewGeneNR table
Kruskal<-kruskal.test(GeigerGenesYO[,3], factor(GeigerGenesYO[,2]))
KruskalTable<-as.data.frame(Kruskal$p.value)
names(KruskalTable)<-names(GeigerGenesYO)[3]
row.names(KruskalTable)<-"Kruskal"

for (i in 4:length(colnames(GeigerGenesYO))){
  x<-kruskal.test(GeigerGenesYO[,i], factor(GeigerGenesYO[,2]))
  y<-as.data.frame(x$p.value)
  names(y)<-names(GeigerGenesYO)[i]
  KruskalTable<-cbind(KruskalTable, y)
}

KruskalTable<-as.data.frame(t(KruskalTable))
names(KruskalTable)<-c("Unadjusted_p")
KruskalTable$FDR<-p.adjust(KruskalTable$Unadjusted_p, method = "fdr")
KruskalTable$Gene<-row.names(KruskalTable)
KruskalTable<-KruskalTable[order(KruskalTable$Unadjusted_p),]

GeigerYOGeneTable<-merge(GeigerYOGeneTable, KruskalTable, by = "Gene", all.x = TRUE)

GeigerYOGeneTable2<-subset(GeigerYOGeneTable, GeigerYOGeneTable$Unadjusted_p < 0.05)
GeigerYOGeneTable2<-GeigerYOGeneTable2[order(GeigerYOGeneTable2$logratio, decreasing = FALSE),]
GeigerYOGeneTable3<-subset(GeigerYOGeneTable, GeigerYOGeneTable$FDR < 0.2)
GeigerYOGeneTable<-GeigerYOGeneTable[order(GeigerYOGeneTable$logratio, decreasing = FALSE),]

GeigerYOGeneTable$Gene<-gsub("\t", " ", GeigerYOGeneTable$Gene)
GeigerYOGeneTable$Gene<-gsub("   ", " ", GeigerYOGeneTable$Gene)
GeigerYOGeneTable$Gene<-gsub("  ", " ", GeigerYOGeneTable$Gene)
GeigerYOGeneTable$Gene<-gsub(",", "", GeigerYOGeneTable$Gene)

GeigerYOGeneTable2$Gene<-gsub("\t", " ", GeigerYOGeneTable2$Gene)
GeigerYOGeneTable2$Gene<-gsub("   ", " ", GeigerYOGeneTable2$Gene)
GeigerYOGeneTable2$Gene<-gsub("  ", " ", GeigerYOGeneTable2$Gene)
GeigerYOGeneTable2$Gene<-gsub(",", "", GeigerYOGeneTable2$Gene)



matrixtable<-GeigerYOGeneTable2[,6:7]
row.names(matrixtable)<-GeigerYOGeneTable2$Gene
matrixtable<-as.data.frame(matrixtable)
matrixtable<-subset(matrixtable, matrixtable$logratio < -1 | matrixtable$logratio > 1)
matrixtable<-as.matrix(matrixtable)
heatmap.2(matrixtable, Rowv = NULL, Colv = FALSE, col = diverging.colors(100), trace = "none",
          tracecol= "black",hline = NULL, vline = NULL, margins=c(2,25),
          xlab = NULL, main = NULL)


heatmap.2(t(matrixtable), Rowv = NULL, Colv = FALSE, col = diverging.colors(100), trace = "none",
          tracecol= "black",hline = NULL, vline = NULL, margins=c(15,2),
          xlab = NULL, main = NULL)

### Genes in YO WT

GeigerGeneYOWT<-subset(GeigerGenes, GeigerGenes$Groups %in% c ("Young_C57BL-Untransplanted", "Old_C57BL-Untransplanted"))

GeigerGeneWTYO.mrpp<-mrpp(GeigerGeneYOWT[,3:length(colnames(GeigerGeneYOWT))], GeigerGeneYOWT[,2], distance = "euclidean")
GeigerGeneWTYO.mrpp



GeigerGeneTableWTYO<-as.data.frame(t(aggregate( . ~ Groups, data =GeigerGeneYOWT[,2:ncol(GeigerGeneYOWT)], FUN = mean)))
GeigerGeneTableWTYO<-GeigerGeneTableWTYO[-1,]
GeigerGeneTableWTYO<-as.data.frame(data.matrix(GeigerGeneTableWTYO))
colnames(GeigerGeneTableWTYO)<-c("Young_C57BL_Untransplanted", "Old_C57BL_Untransplanted")

GeigerGeneTableWTYO$Difference<-(GeigerGeneTableWTYO$Young_C57BL_Untransplanted -GeigerGeneTableWTYO$Old_C57BL_Untransplanted)
GeigerGeneTableWTYO$FC<-foldchange(GeigerGeneTableWTYO$Young_C57BL_Untransplanted, GeigerGeneTableWTYO$Old_C57BL_Untransplanted)
GeigerGeneTableWTYO$logratio<-foldchange2logratio(GeigerGeneTableWTYO$FC)
# GeigerGeneTableWTYO$log2Delta<-ifelse(sign(GeigerGeneTableWTYO$Difference) == -1, -GeigerGeneTableWTYO$log2Delta,
#                                       ifelse(GeigerGeneTableWTYO$Difference == 0, 0, GeigerGeneTableWTYO$log2Delta))
GeigerGeneTableWTYO<-GeigerGeneTableWTYO[order(GeigerGeneTableWTYO$logratio, decreasing = FALSE),]
GeigerGeneTableWTYO$Gene<-row.names(GeigerGeneTableWTYO)

# Ok that's cool. Now figure out which are statistically signficant from the GeigerNewGeneNR table
Kruskal<-kruskal.test(GeigerGeneYOWT[,3], factor(GeigerGeneYOWT[,2]))
KruskalTable<-as.data.frame(Kruskal$p.value)
names(KruskalTable)<-names(GeigerGeneYOWT)[3]
row.names(KruskalTable)<-"Kruskal"

for (i in 4:length(colnames(GeigerGeneYOWT))){
  x<-kruskal.test(GeigerGeneYOWT[,i], factor(GeigerGeneYOWT[,2]))
  y<-as.data.frame(x$p.value)
  names(y)<-names(GeigerGeneYOWT)[i]
  KruskalTable<-cbind(KruskalTable, y)
}

KruskalTable<-as.data.frame(t(KruskalTable))
names(KruskalTable)<-c("Unadjusted_p")
KruskalTable$FDR<-p.adjust(KruskalTable$Unadjusted_p, method = "fdr")
KruskalTable$Gene<-row.names(KruskalTable)
KruskalTable<-KruskalTable[order(KruskalTable$Unadjusted_p),]

GeigerGeneTableWTYO<-merge(GeigerGeneTableWTYO, KruskalTable, by = "Gene", all.x = TRUE)

GeigerGeneTableWTYO2<-subset(GeigerGeneTableWTYO, GeigerGeneTableWTYO$Unadjusted_p < 0.05)
GeigerGeneTableWTYO2<-GeigerGeneTableWTYO2[order(GeigerGeneTableWTYO2$logratio, decreasing = FALSE),]
GeigerGeneTableWTYO3<-subset(GeigerGeneTableWTYO, GeigerGeneTableWTYO$FDR < 0.2)
GeigerGeneTableWTYO<-GeigerGeneTableWTYO[order(GeigerGeneTableWTYO$logratio, decreasing = FALSE),]

GeigerGeneTableWTYO2$Gene<-gsub("\t", " ", GeigerGeneTableWTYO2$Gene)
GeigerGeneTableWTYO2$Gene<-gsub("   ", " ", GeigerGeneTableWTYO2$Gene)
GeigerGeneTableWTYO2$Gene<-gsub("  ", " ", GeigerGeneTableWTYO2$Gene)
GeigerGeneTableWTYO2$Gene<-gsub(",", "", GeigerGeneTableWTYO2$Gene)



matrixtable<-GeigerGeneTableWTYO2[,6:7]
row.names(matrixtable)<-GeigerGeneTableWTYO2$Gene
matrixtable<-as.data.frame(matrixtable)
matrixtable<-subset(matrixtable, matrixtable$logratio < -1 | matrixtable$logratio > 1)
matrixtable<-as.matrix(matrixtable)
heatmap.2(matrixtable, Rowv = NULL, Colv = FALSE, col = diverging.colors(100), trace = "none",
          tracecol= "black",hline = NULL, vline = NULL, margins=c(2,25),
          xlab = NULL, main = NULL)


heatmap.2(t(matrixtable), Rowv = NULL, Colv = FALSE, col = diverging.colors(100), trace = "none",
          tracecol= "black",hline = NULL, vline = NULL, margins=c(15,2),
          xlab = NULL, main = NULL)


save.image(file="GeigerData20220927")

# Now Gene differences (skip this and use Kruskall instead)

# Wilcox<-pairwise.wilcox.test(GeigerGenes[,3], GeigerGenes[,2])
# WilcoxTable<-as.data.frame(Wilcox$p.value)
# names(WilcoxTable)<-names(GeigerGenes)[3]
# row.names(WilcoxTable)<-"Wilcox"
# 
# 
# for (i in 4:length(colnames(GeigerGenes))){
#   x<-pairwise.wilcox.test(GeigerGenes[,i], GeigerGenes[,2])
#   y<-as.data.frame(x$p.value)
#   names(y)<-names(GeigerGenes)[i]
#   WilcoxTable<-cbind(WilcoxTable, y)
# }
# 
# WilcoxTable<-as.data.frame(t(WilcoxTable))
# names(WilcoxTable)<-c("Unadjusted_p")
# WilcoxTable$FDR<-p.adjust(WilcoxTable$Unadjusted_p, method = "fdr")
# WilcoxTable$Species<-row.names(WilcoxTable)
# WilcoxTable<-WilcoxTable[order(WilcoxTable$Unadjusted_p),]
# names(WilcoxTable)<-c("Unadjusted_p", "FDR", "Gene")
# 
# write.csv(WilcoxTable, file = "GeigerMicrobialGeneWilcoxTable.csv")
# sigGenes<-subset(WilcoxTable, WilcoxTable$Unadjusted_p < 0.05)
# write.csv(sigGenes, file = "SignificantGeigerGenes.csv")
# 
# GeigerGenes.mrpp<-mrpp(GeigerGenes[,3:length(colnames(GeigerGenes))], GeigerGenes[,2], distance = "euclidean")
# GeigerGenes.mrpp

save.image(file = "GeigerData20220927")


