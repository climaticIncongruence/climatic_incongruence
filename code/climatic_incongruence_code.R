
#### Load required packages
library(ggplot2)
library(plyr)
library(dplyr)
library(vegan)
library(dendextend)
library(cowplot)
library(parallel)

####%%%%% Cluster analysis
clim <- read.csv("../data/clim_inc_dat.csv")
climVar <- clim[, 5:17] 
rownames(climVar) <- clim$Population
climDist <- dist(climVar)
climClust <- hclust(climDist)

## Call the two partitions...
CC <- cutree(climClust, 7)
SZ <- clim$SZ
Pop <- clim$Population

## Generate colour palettes for the two partitions (following paper)
ccPalette <- c("#8DD3C7", "#E988E9", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69")
sz.pal <- c("#E41A1C", "#000099", "#00CC33", "#984EA3", "#FF0000", "#000000", "#FFCC00")

#### Make the dendrogram (Fig 1...) - code to produce map and the 'interaction plot not included here - needed?

cl <- as.dendrogram(climClust)
prettydend <- cl %>% set("branches_k_color", k = 7) %>% 
  set("branches_lwd", 2) %>% set("labels_cex", 0.65)
d1 <- color_branches(prettydend, 7, col = c("#8DD3C7", "#FDB462","#E988E9","#BEBADA", "#B3DE69","#FB8072", "#80B1D3"))

par(bg = NA)
plot(d1, horiz = TRUE, axes = FALSE)

####%%%%% Make the map

Points <- cbind(select(clim, Eastings, Northings), CC)
seed_zone_lines <- read.csv("../data/SP_biochemical_regions.csv")

ggplot(seed_zone_lines, aes(long, lat, group = group))+
  geom_path()+
  coord_fixed(ratio = 1) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  theme(panel.background = element_blank()) +
  geom_point(data = Points, aes(Eastings, Northings, group = factor(CC), fill = factor(CC)),  size = 4.5, pch = 21)+
  scale_fill_manual(values = ccPalette) +
  guides(fill = FALSE)

####%%%%%  Draw the lines connecting dendrogram to seed zone partition
## This bit was put together with the dendrogram (produced above) for final figure using Inkscape (https://inkscape.org/)

rnks <- read.csv("../data/SZ_CC_ranks.csv")

ggplot(rnks, aes(Rank, Score, group = Name_CC, fill = SZ))+
  geom_point(size = 4, pch = 21)+
  geom_line()+
  guides(fill = FALSE) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  theme(panel.background = element_blank())

####%%%%% Comparison of classifiers.
#### The idea here is to compare two classifiers and determine whether they agree more than is to be to expected
#### by chance
#### We have N observations (individuals) divided into K classes
#### Since the classes in the classifiers are arbitrary, the assignments can't be compared directly.
#### Instead, we will transform each vector with the class assignments into a membership matrix, M
#### The element M[i, j] is 1 when i and j are in the same class and 0 otherwise
#### We then proceed to compare the membership matrices for the two classifiers.
#### Here, the function ComputeStat() will compute an agreement measure between the two membership matrices
#### Currently, it counts how entries are in agreement. (Note: the minimum is N and the maximum is N^2)
#### Next step is to compute the agreement "statistic" between one of the classifiers and 'Nsim' simulated
#### classifiers, with the same number of individuals, classes and class frequencies.
#### This re-sampling step allows us to assess the "significance" of the "statistic"  we computed
#### for our two classifiers, stored in the 'realStat' object.
#### Our code essentially implements the agreement index devised by Rand, 1971.
############################

## Begin Functions
create.bin.matrix <- function(cl){
  ## takes one classifier and returns a binary agreement matrix
  n <- length(cl)
  Grid <- expand.grid(1:n, 1:n)
  res <- apply(Grid, 1, function(pair){
    ifelse(cl[pair[1]] == cl[pair[2]], 1, 0)
  })
  return(data.frame(Grid, ind = res))
}
# ComputeStat <- function(agree1, agree2){ ## alternative (dis)agreement measure
#   sum(agree1 - agree2)
# }
ComputeStat <- function(agree1, agree2){
  ## computes the agreement "statistic" between membership matrices
  sum(agree1 == agree2)
}
nclasses <- function(x) length(unique(x))
classProbs <- function(x) table(x)/sum(table(x))
Simulate <- function(i, N, K, Pp) create.bin.matrix(sample(1:K, N, prob = Pp, replace = TRUE))
## End Functions
################
################
### Simulate some data
set.seed(666)
N <- 84
K <- 7

dt <- data.frame(Pop, c1 = CC, c2 = SZ)
MatC1 <- create.bin.matrix(dt$c1)
MatC2 <- create.bin.matrix(dt$c2)
ProbsC2 <- classProbs(dt$c2) ## notice we compute frequencies for the SECOND classifier

#### Re-sampling step. Generating 'Nsim' random classifications compatible with the SECOND classifier
Nsim <- 10000
Ncores <- 10
RandomStats <- unlist(parallel::mclapply(1:Nsim, function(i){
  ComputeStat(agree1 = MatC1$ind, agree2 = Simulate(1, N = N, K = K, Pp = ProbsC2)$ind)
}, mc.cores = Ncores))

realStat <- ComputeStat(agree1 = MatC1$ind, agree2 =  MatC2$ind)
## Proportion of randomly generated data sets that result in a greater "statistic" than the observed one
mean(RandomStats > realStat)

###Make the picture...

jPerm <- data.frame(RandomStats)

Perm <- ggplot(jPerm, aes(RandomStats/N^2)) +
  geom_histogram(colour = "black") +
  geom_vline(xintercept = realStat/N^2, lty = "dashed") +
  theme_bw(base_size = 6) +
  theme(panel.grid = element_blank()) +
  theme(panel.background = element_blank()) +
  labs(x = "Rand index", y = "Count")

Perm
# # Optionally, run with clusterCrit
# install.packages("clusterCrit")
# library(clusterCrit)
# extIdx <- extCriteria(part1 = as.integer(dt$c1), part2 = as.integer(dt$c2),"rand")
# extIdx
# # if you're feeling adventurous you can also do
# extCriteria(part1 = as.integer(dt$c1), part2 = as.integer(dt$c2), "all")
# #and have fun!

####%%%%% Do the PCA
pc <- princomp(climVar, cor = TRUE)
pcScores <- data.frame(scores(pc))
pcScores <- cbind(pcScores[, 1:3], SZ, factor(CC))
names(pcScores) = c("Comp.1", "Comp.2", "Comp.3", "SZ", "CC")

pc1 <- pcScores$Comp.1
pc2 <- pcScores$Comp.2
pc3 <- pcScores$Comp.3

## Generate hulls for PCA biplot for the two partitions

chulls <- ddply(pcScores, .(CC), function(df) df[chull(df$Comp.1, df$Comp.2), ])
chulls_SZ <- ddply(pcScores, .(SZ), function(df) df[chull(df$Comp.1, df$Comp.2), ])

## Make the PCA biplots
CC <- ggplot(pcScores, aes(Comp.1, Comp.2, fill = CC)) +
  geom_polygon(data = chulls, aes(Comp.1, Comp.2, group = CC), colour = "black", alpha = 0.45)+
  geom_point(size = 2, pch = 21) +
  scale_fill_manual(values = ccPalette, guide = guide_legend(title = "Climatic clusters")) +
  theme_bw(base_size = 6) +
  theme(panel.grid = element_blank()) +
  ggtitle("b. Climatic clusters") +
  theme(plot.title = element_text(hjust = 0)) +
  labs(x = "PC1 (61%)", y = "PC2 (24%)")

SZ <- ggplot(pcScores, aes(Comp.1, Comp.2, fill = SZ)) +
  geom_polygon(data = chulls_SZ, aes(Comp.1, Comp.2, group = SZ), colour = "black", alpha = 0.45) +
  geom_point(size = 2, pch = 21) +
  scale_fill_manual(values = sz.pal) +
  guides(fill = guide_legend(title = "Seed zones")) +
  theme_bw(base_size = 6) +
  theme(panel.grid = element_blank()) +
  ggtitle("a. Seed zones") +
  theme(plot.title = element_text(hjust = 0)) +
  labs(x = "PC1 (61%)", y = "PC2 (24%)")

## Plot biplots side by side (requires package:cowplot) ::: Figure 2
plot_grid(SZ, CC, align = 'h')

####%%%%% Make figure 3 (PC3 boxplots)
ggplot(pcScores, aes(CC, Comp.3, fill = CC))+
  geom_boxplot()+ scale_fill_manual(values = ccPalette) + 
  guides(fill = FALSE)+
  labs(x = "Climatic cluster", y = "PC3 (7%)")+
  theme_bw(base_size = 6)+
  theme(panel.grid = element_blank())

####%%%%%  Do the ANOVAs
lm(pc1 ~ CC, data = pcScores) %>% summary() 
lm(pc2 ~ CC, data = pcScores) %>% summary()
lm(pc1 ~ SZ, data = pcScores) %>% summary()
lm(pc2 ~ SZ, data = pcScores) %>% summary()
