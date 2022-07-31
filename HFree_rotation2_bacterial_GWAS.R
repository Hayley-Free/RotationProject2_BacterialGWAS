## Rotation Project 2: Bacterial Genome Wide Association Study (GWAS) ##
# 31/07/2022
# Hayley Free - BBSRC Oxford Interdisiplinary Biosicence DTP
# Supervisor: Dr Daniel Wilson
# Big Data Institute, Old Road Campus

# This is the R code used for the multiple parts of analysis for bacterial GWAS
# on Campylobacter coli host specificity to humans, chickens and pig sources.

# Files used:
# metadata_4329_ccoli.xlsx - the original metadata as downloaded from PubMLST
# metadata_fastaname_new.csv - the metadata with edited 'fasta name' column
# pca_source.csv - source metadata for multidimensional scaling plots
# tree_extra_data.xlsx - tidied metadata for labelling tree plots
# sender20_nucleotide31.patternmerge.patternIndex.txt.gz
# sender20_pvals/ppattern.txt - all p_LRT values for sender 20
# sender20_nucleotide31.kmermerge.txt.gz - kmer file
# sender225_nucleotide31.kmermerge.txt.gz - kmer file
# sender225_nucleotide31.patternmerge.patternIndex.txt.gz 
# sender225_pvals/ppattern.txt - all p_LRT values for sender 225

# Index
# Part 1 - Metadata
# Part 2 - Mash Data
# Part 3 - Multidimensional Scaling
# Part 4a - NJ Trees All Data
# Part 4b - NJ Trees Chicken & Pig Data
# Part 4c - NJ Trees Human & Chicken Data
# Part 5 - Commonality

# Packages required:
library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(readxl)
library(wesanderson)
#install.packages("remotes")
#remotes::install_github("YuLab-SMU/ggtreeExtra")
library(ggtreeExtra)
library(ggnewscale)
#install.packages("PCMBase")
library(PCMBase)
#install.packages("ggpubr")
library("ggpubr")
#install.packages("cowplot")
library(cowplot)

# The working directory is set to where the files are sourced and saved
setwd("~/OneDrive - Oxford Brookes University/Hayley_2022/Rotation_2_GWAS/c_coli_data")


#### Part 1 - Metadata ####
# The meta data downloaded from PubMLST for Campylobacter coli (see write up for
# filter details). This was uploaded to use to annotate plots. 

# The original metadata:
orig_metadata_4329_ccoli <- read_excel("metadata_4329_ccoli.xlsx")
# Creating new metadata file
metadata_fasta_name <- orig_metadata_4329_ccoli
# Creating new column of ID_isolate - which is the name of the fasta files
metadata_fasta_name$fasta_name <- paste(metadata_fasta_name$id, metadata_fasta_name$isolate,sep="_")
# The new column had to be manually inspected in Excel as they were not correct
# Updated fasta names - new column called fasta_name_correct
metadata_fastaname_new <- read.csv("~/OneDrive - Oxford Brookes University/Hayley_2022/Rotation_2_GWAS/c_coli_data/metadata_fastaname_new.csv")
# Shortening the dataframe to columns needed 
metatake2 <- metadata_fastaname_new[,c(30,15)]
# Filtering for the sources
metatake2<- dplyr::filter(metatake2, source == "human stool" | source == "chicken" | source == "pig")
# Checking the source are human, chicken and pig only
table(metatake2$source)


#### Part 2 - Mash Data ####
# The Mash algorithm was used to create a distance matrix of relatedness.

# Upload a copy of the mash output with no row names but column names as headers
distance_fix5_copy <- read.delim("~/OneDrive - Oxford Brookes University/Hayley_2022/Rotation_2_GWAS/c_coli_data/c_coli_contigs/another4_correct/distance_fix5_copy.tsv")
# Original mash output but no headers set
distance_fix5 <- read.delim("~/OneDrive - Oxford Brookes University/Hayley_2022/Rotation_2_GWAS/c_coli_data/c_coli_contigs/another4_correct/distance_fix5.tsv", header=FALSE, comment.char="#")
# Getting the row names to create the full mash output
distance_fix5[,1]
row_names <- as.vector(distance_fix5[,1])
row.names(distance_fix5_copy)<-row_names
# Remove these variables
rm(row_names, distance_fix5)
# Cleaning the names
rownames(distance_fix5_copy)<-gsub("folder1/","",rownames(distance_fix5_copy))
rownames(distance_fix5_copy)<-gsub(".fas","",rownames(distance_fix5_copy))
colnames(distance_fix5_copy)<-gsub(".fas","",colnames(distance_fix5_copy))
colnames(distance_fix5_copy)<-gsub("folder1.","",colnames(distance_fix5_copy))


#### Part 3 - Multidimensional Scaling ####
# Multidimensional scaling (MDS) was used to identify outliers in the contigs
# and to identify the clustering of the C. coli contigs

# Using cmdscale for MDS: returns the best-fitting dimensional representation,
# takes n by n matrix and returns n by p configuration matrix. p is the dimension
# of the smallest space in which the n points whose interpoint distances are given
# by the matrix can be embedded.
pca_all_outliers <- cmdscale(distance_fix5_copy)
# Source metadata for multidimensional scaling plots
pca_source_csv <- read.csv("~/OneDrive - Oxford Brookes University/Hayley_2022/Rotation_2_GWAS/c_coli_data/pca_source.csv", row.names=1)
pca_source_csv$source = as.factor(pca_source_csv$source)
# Plot first MDS
plot(pca_all_outliers,col=pca_source_csv$source, pch=16, cex=0.75, ylab="V1",
     xlab="V2")
legend("topright",legend=levels(factor(pca_source_csv$source)), pch=19,
       col = seq_along(levels(factor(pca_source_csv$source))))

# Finding the means of the matrix to get an average relatedness
means <- rowMeans(distance_fix5_copy)
# Below plot is means of matrix - the average relatedness of each contig
barplot(means,
        xlab = "C. coli contigs",
        ylab = "Average genetic distance",
        las=2,
        cex.axis=1,
        xaxt='n',
        cex.lab=1,
        ylim=c(0,1.1))
text(898, 0.35, "*")
text(3200, 1.05, "*")
# This graphs show outliers - those with a larger distance between them.
# From the graph, creat a cut off for outlier range seen
sum(means>0.15) 
outliers_first <- which(means>0.15)
# True and false of all values - to use as a sum to see how many are included
gd_first <- means<0.15
# MDS of data with outliers over 0.15 cut off removed
pca_one_outlier <- cmdscale(distance_fix5_copy[-outliers_first,-outliers_first])
# MDS with outliers over 0.15 cut off removed
plot(pca_one_outlier,col=pca_source_csv$source, pch=16, cex=0.75, ylab="V1",
     xlab="V2")
legend("topleft",legend=levels(factor(pca_source_csv$source)), pch=19,
       col=seq_along(levels(factor(pca_source_csv$source)))) 
# From the plot there is still one outlier left in top right
# Creat bar plot to find cut off for outlier
final_outlier_plot <- barplot(means[gd_first],
                              ylab = "Average genetic distance",
                              xlab = "C. coli contigs",
                              las=2,
                              cex.axis=1,
                              xaxt='n',
                              cex.lab=1,
                              ylim=c(0,0.1029))
text(309, 0.102, "*")
# New cut off
sum(means>0.08)
outliers <- which(means>0.08)
gd<- means<0.08
# Creating test3 to add colour with plot and cmdscale
test3 <- cmdscale(distance_fix5_copy[-outliers,-outliers])
# Final MDS
plot(test3,col=pca_source_csv$source, pch=16, cex=0.75, ylab="V1", xlab="V2")
legend("topleft",legend=levels(factor(pca_source_csv$source)), pch=19,
       col=seq_along(levels(factor(pca_source_csv$source))))





#### Part 4a - NJ Trees All Data ####
# The distance matrix (mash output) is used to create neighbour-joining (NJ) trees.
# These consist of the total data and are then filtered for the successive phenotype
# and corresponding metadata. 

# Removing outliers from distance matrix found in MDS
no_outliers_distance <- distance_fix5_copy[-outliers,-outliers]
# Making a matrix and then NJ using ape
no_out_m_distance_fix5 <- as.matrix(no_outliers_distance)
no_out_nj_distance_fix5 <- ape::nj(as.dist(no_out_m_distance_fix5))

# The metadata was cleaned in Excel for easy labelling of trees
tree_extra_data <- read_excel("tree_extra_data.xlsx")

# The values for the branches (edge length) is taken as an absolute value:
abs_adjst <- abs(no_out_nj_distance_fix5$edge.length)
abs_no_out_nj_dist_fix5 <- no_out_nj_distance_fix5
abs_no_out_nj_dist_fix5$edge.length <-abs_adjst
# and square rooted to minimise the distances between points and negate the
# negative branch lengths:
sqt_abs <- sqrt(abs_no_out_nj_dist_fix5$edge.length)
abs_no_out_nj_dist_fix5_sqt <- abs_no_out_nj_dist_fix5
abs_no_out_nj_dist_fix5_sqt$edge.length <- sqt_abs

# The original full chicken, human and pig sourced C. coli contigs
# geom_fruit() is used to plot additional metadata in rings around the tree
layer1 <- ggtree(abs_no_out_nj_dist_fix5_sqt,layout="fan")%<+% metatake2 + 
  geom_tippoint(aes(color=source), size=0.5)+  scale_color_manual(values=c(wes_palette("GrandBudapest1", n = 3)))
layer2 <- layer1 + new_scale_fill() +
  geom_fruit(
    data=tree_extra_data,
    geom=geom_tile,
    mapping=aes(y=tip_label, fill=country),
    width=0.1,
    offset=0.1,
    pwidth=0.2)+
  scale_fill_manual(
    name="Country",
    values=wes_palette("Darjeeling2", 31, type = "continuous"))
layer3 <- layer2 + 
  new_scale_fill() +
  geom_fruit(
    data=tree_extra_data,
    geom=geom_tile,
    mapping=aes(y=tip_label, fill=as.character(sender)),
    width=0.1,
    offset=0.1,
    pwidth=0.2)+
  scale_fill_manual(
    name="Sender",
    values=wes_palette("Moonrise3", n = 44,type = "continuous"))
layer3


#### Part 4b - NJ Trees Chicken & Pig Data ####
# According to the phenotypes being investigated in GWAS, the human-stool
# sourced contigs are removed from the tree. 
# Create a tree variable of initial data
treeo <- ggtree(abs_no_out_nj_dist_fix5_sqt,layout="fan")
# Get human stool tip names 
human_tip <- tree_extra_data$tip_label[tree_extra_data$source=="human stool"]
# Use drop.tip with human stool tip names to drop them from the tree variable
treeo_1 <- drop.tip(as.phylo(treeo), tip=human_tip)

# Inspecting the country and sender of the chicken and pig sourced contigs
table(tree_extra_data$country[tree_extra_data$source!="human stool"])
table(tree_extra_data$sender[tree_extra_data$source!="human stool"])

# Adding rings of metadata to chicken and pig sourced contigs in NJ tree
layerA <- ggtree(treeo_1,layout="fan")%<+% metatake2 + 
  geom_tippoint(aes(color=source), size=0.5)+
  scale_color_manual(values=c(wes_palette("GrandBudapest1", n = 2)))
layerB <- layerA +
  new_scale_fill() +
  geom_fruit(
    data=tree_extra_data,
    geom=geom_tile,
    mapping=aes(y=tip_label, fill=country),
    width=0.1,
    offset=0.1,
    pwidth=0.2)+
  scale_fill_manual(
    name="Country",
    values=wes_palette("Darjeeling2", 22, type = "continuous"))
layerC <- layerB +
  new_scale_fill() +
  geom_fruit(
    data=tree_extra_data,
    geom=geom_tile,
    mapping=aes(y=tip_label, fill=as.character(sender)),
    width=0.1,
    offset=0.1,
    pwidth=0.2)+
  scale_fill_manual(
    name="Sender",
    values=wes_palette("Moonrise3", n = 31,type = "continuous"))
layerC

# The data was now filtered by sender - the contributor of the contigs in PubMLST
# This was to remove confounding from batch effect (e.g. differences in study)

# The table of sender in tree without human i.e. sender of chicken and pig contigs
# shows which sender the majority of contigs are from
table(tree_extra_data$sender[tree_extra_data$source!="human stool"])
# There are four senders with the majority of contigs: 20, 216, 225 and 438
# Data frame of just chicken and pig
cvp <- tree_extra_data[tree_extra_data$source == "chicken"| 
                         tree_extra_data$source == "pig",]
# There should be 1252 chicken and pig sequences

# Creating bar plot of country to show that the majority are from a select few
cvptablecountry <- as.data.frame(table(cvp$country))
ggplot(cvptablecountry, aes(x=Var1, y=Freq, fill=Var1))+geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(name="Country",values=wes_palette("Darjeeling2", 22,
                                                      type = "continuous")) + 
  xlab("Country") + ylab("Frequency")+
  theme(legend.position="none")+
  ggtitle("Country of C. coli contigs of chicken and pig source")


# Dropping senders which are not 20, 216, 225 or 438
# Create tree variable
treeo_2 <- ggtree(treeo_1,layout="fan")
# Create list of senders hat are not 20, 216, 225, 438 tip names
no_senders <- cvp$tip_label[cvp$sender == "12" | cvp$sender == "89" | 
                              cvp$sender == "99" | cvp$sender == "119" | 
                              cvp$sender == "126"| cvp$sender == "167"| 
                              cvp$sender == "174"| cvp$sender == "197"| 
                              cvp$sender == "213"| cvp$sender == "255"| 
                              cvp$sender == "269"| cvp$sender == "319"| 
                              cvp$sender == "328"| cvp$sender == "338"| 
                              cvp$sender == "346"| cvp$sender == "385"| 
                              cvp$sender == "412"| cvp$sender == "423"| 
                              cvp$sender == "444"| cvp$sender == "460"| 
                              cvp$sender == "508"| cvp$sender == "512"| 
                              cvp$sender == "514"| cvp$sender == "564"| 
                              cvp$sender == "565"| cvp$sender == "568"| 
                              cvp$sender == "571"] 
# Use drop.tip with tip names
treeo_3 <- drop.tip(as.phylo(treeo_2), tip=no_senders)

ggtree(treeo_3,layout="fan")%<+% metatake2 + 
  geom_tippoint(aes(color=source), size=1)+ 
  scale_color_manual(values=c("#FCC214", "#FF7F50"),name="Source") + 
  new_scale_fill() +
  geom_fruit(
    data=tree_extra_data,
    geom=geom_tile,
    mapping=aes(y=tip_label, fill=as.character(sender)),
    width=0.1,
    offset=0.1,
    pwidth=0.2)+
  scale_fill_manual(
    name="Sender",
    values=c(wes_palette("Moonrise3", n = 4)))
# The sender 20 and sender 225 seem to cover both chicken and pig sourced contigs
# Contigs from these two sender will be used in separate GWAS



#### Part 4c - NJ Trees Human & Chicken Data ####
# Another phenotype studies was human and chicken sourced specificity of C. coli
# The pig sourced contigs needed removing from the tree, including any chicken 
# and human sourced contigs which clustered with the clade (to reduce population
# stratification)

# Tree with all data and with node labels (of number)
all <- ggtree(abs_no_out_nj_dist_fix5_sqt,layout="fan")%<+% metatake2 + 
  geom_tippoint(aes(color=source), size=0.5)+ 
  scale_color_manual(values=c(wes_palette("GrandBudapest1", n = 3))) + 
  geom_text(aes(label=node), hjust=-.3, size=2.5)
# Highlighting the node which needs to be removed
ggtree(abs_no_out_nj_dist_fix5_sqt,layout="fan")%<+% metatake2 + 
  geom_tippoint(aes(color=source), size=0.5)+ 
  scale_color_manual(values=c(wes_palette("GrandBudapest1", n = 3)))+
  geom_hilight(node=3412, fill="yellow")
# Upwards from the node want removed is 3412

# Extracting the pig clade to attain which contig IDs are present using PCMBase
# Making PCM tree format of all of the data
pmctree <- PCMTreePlot(as.phylo(all), 
                       palette=c(a = "red", b = "green", c = "blue")) +
  ggtree::geom_nodelab(angle = 45) + ggtree::geom_tiplab(angle = 45)
# Selecting the pig clade
blueTree <- PCMTreeExtractClade(as.phylo(pmctree), 3412)
PCMTreeGetPartRegimes(blueTree)
# Forming a tree variable with the pig clade
pmctree_pigbranch <- PCMTreePlot(blueTree,
                                 palette=c(a = "red", b = "green", c = "blue")) +
  ggtree::geom_nodelab(angle = 45) + ggtree::geom_tiplab(angle = 45)
# Filtering the pig clade tree for the data
pmctree_pigbranch_data <- pmctree_pigbranch$data
# Writing the pig clade tree data out to get the names
write.csv(pmctree_pigbranch_data,"pmctree_pigbranch_data.csv")
# got these names through filtering in excel but could use pmctree_pigbranch_data$label
hctoremove <- read_excel("hctoremove.xlsx")
# Could also attain the names through 
# hctoremove <- pmctree_pigbranch_data$label
# Getting data to correct format to use in drop.tip
rmhc<- as.list(hctoremove)
huch_inpig <- rmhc$fasta_name_correct
# Using drop.tip to remove clade
# Tree variable
treeo_4 <- ggtree(abs_no_out_nj_dist_fix5_sqt,layout="fan")
# Usig drop.tip with tip names
treeo_5 <- drop.tip(as.phylo(treeo_4), tip=huch_inpig)
# Current tree still needs pig sourced contigs not in main clade removed 
# Using drop.tip again
# Create tree variable
treeo_6 <- ggtree(treeo_5,layout="fan")
# Attain pig tip names
pig_tip <- tree_extra_data$tip_label[tree_extra_data$source=="pig"]
# Use drop.tip with tip names
treenp <- drop.tip(as.phylo(treeo_6), tip=pig_tip)
# Plotting tree - no pigs and no humans and chickens which were in the pig clade:
nopigtree <- ggtree(treenp,layout="fan")%<+% metatake2 + 
  geom_tippoint(aes(color=source), size=1) + scale_color_manual(values=c(wes_palette("GrandBudapest1", n = 2)),name="Source")
nopigtree
# Attaining the IDs of the humans and chickens in the tree with pig clade which
# were removed
tip_label <- nopigtree$data
tip_label <- tip_label$label
tip_label <- as.data.frame(tip_label)
tip_label <- na.omit(tip_label)
# Combine cvp_tree_extra_data with tip_label to get the metadata
hc_treedata <- merge(tip_label,tree_extra_data, by="tip_label")
# The GWAS is then altered as change what is included in the id file
write.csv(tip_label, "HvP_tips.csv")

# Adding metadata in rings
nopig_sourcetree <- nopigtree +
  new_scale_fill() +
  geom_fruit(
    data=tree_extra_data,
    geom=geom_tile,
    mapping=aes(y=tip_label, fill=as.character(sender)),
    width=0.1,
    offset=0.1,
    pwidth=0.2)+
  scale_fill_manual(
    name="Sender",
    values=(wes_palette("Moonrise3", 36, type="continuous")))
nopig_sourcetree + 
  new_scale_fill() +
  geom_fruit(
    data=tree_extra_data,
    geom=geom_tile,
    mapping=aes(y=tip_label, fill=as.character(country)),
    width=0.1,
    offset=0.1,
    pwidth=0.2)+
  scale_fill_manual(
    name="Country",
    values=(wes_palette("Darjeeling2", 17, type="continuous")))

# Showing the trees before and after the pig clade (and other pig contigs)
# were removed, with the clade highlighted
ggarrange((ggtree(abs_no_out_nj_dist_fix5_sqt,layout="fan")%<+% metatake2 + 
             geom_tippoint(aes(color=source), size=0.5)+ scale_color_manual(values=c(wes_palette("GrandBudapest1", n = 3)))+geom_hilight(node=3412, fill="yellow")), 
          nopigtree,labels = c("All Phenotypes", "Pig Clade Removed"),
          ncol = 2, nrow = 1)



#### Part 5 - Commonality ####
# The commonality across the GWAS of sender 20 and sender 225 contigs was assessed
# by plotting the pairs of p-values for top kmers for both studies
# Files are used from the bugwas pipeline:
# The *.patternmerge.patternIndex.txt.gz file describes the presence or absence
# of each pattern (unique patterns of kmer presence (1) or absence (0) across the
# sample)
# The *.kmermerge.txt.gz file contains all kmers present at least once in the
# samples in the ID file
# The *.assoc.txt.gz files in the GEMMA output directory containing the pattern
# number, p-values and log likelihood. These files were used to attain the
# likelihood ratio test p-value (p_LRT) per kmer, which are initially described
# per pattern in the *.assoc.txt.gz files.



# Sender 20
# Presence/absence of kmer
index20 <- read.table("~/OneDrive - Oxford Brookes University/Hayley_2022/Rotation_2_GWAS/c_coli_data/pipeline_files/sender20_pvals/sender20_nucleotide31.patternmerge.patternIndex.txt.gz", quote="\"", comment.char="")
# Empty list same length as above (number of kmers)
pkmer20 <- seq(1,6689727,1)
# The p_LRT was taken from each *.assoc.txt.gz and combined in Excel
ppattern20 <- read.table("~/OneDrive - Oxford Brookes University/Hayley_2022/Rotation_2_GWAS/c_coli_data/pipeline_files/sender20_pvals/ppattern.txt", quote="\"", comment.char="")
# Whenever 'index' is zero, set pkmer to NA replace values in pkmer to NA if 
# value in index is 0
pkmer20[index20$V1==0] = NA
# Changing pkmer20: if index is not 0, pkmer = the ppattern only if index is not 0
pkmer20[index20$V1!=0] = ppattern20$V1[index20$V1[index20$V1!=0]]
# The individual kmers are attached to them
sender20_nucleotide31.kmermerge.txt <- read.table("~/OneDrive - Oxford Brookes University/Hayley_2022/Rotation_2_GWAS/c_coli_data/pipeline_files/sender20_pvals/sender20_nucleotide31.kmermerge.txt.gz", quote="\"", comment.char="")
sender20_pvals <- cbind(sender20_nucleotide31.kmermerge.txt,pkmer20)

# Sender 225
# Presence/absence of kmer
index225 <- read.table("~/OneDrive - Oxford Brookes University/Hayley_2022/Rotation_2_GWAS/c_coli_data/pipeline_files/sender225_pvals/sender225_nucleotide31.patternmerge.patternIndex.txt.gz", quote="\"", comment.char="")
# Empty list same length as above (number of kmers)
pkmer225 <- seq(1,8291751,1)
# The p_LRT was taken from each *.assoc.txt.gz and combined in Excel
ppattern225 <- read.table("~/OneDrive - Oxford Brookes University/Hayley_2022/Rotation_2_GWAS/c_coli_data/pipeline_files/sender225_pvals/ppattern.txt", quote="\"", comment.char="")
# Whenever 'index' is zero, set pkmer to NA replace values in pkmer to NA if 
# value in index is 0
pkmer225[index225$V1==0] = NA
# Changing pkmer20: if index is not 0, pkmer = the ppattern only if index is not 0
pkmer225[index225$V1!=0] = ppattern225$V1[index225$V1[index225$V1!=0]]
# The individual kmers are attached to them
sender225_nucleotide31.kmermerge.txt <- read.table("~/OneDrive - Oxford Brookes University/Hayley_2022/Rotation_2_GWAS/c_coli_data/pipeline_files/sender225_pvals/sender225_nucleotide31.kmermerge.txt.gz", quote="\"", comment.char="")
sender225_pvals <- cbind(sender225_nucleotide31.kmermerge.txt,pkmer225)

# Obtaining which kmers are shared in sender20_pvals and sender225_pvals
# Combining the pvals together to get shared kmers and pvals
mergepvals <- merge(sender20_pvals,sender225_pvals)
# The harmonic mean p-values were taken which detects statistically significant
# of two groups which are individually non-significant
hmp <- 1/(1/mergepvals$pkmer20+1/mergepvals$pkmer225)
# data was ordered to obtain the most significant p_LRT values.
od <- order(hmp,decreasing = F)
# 1 million values which are of interest as ordered by harmonic p-value
subsetpvals <- mergepvals[od[1:1000000],]
# Plotting the scatter plot of kmer p-values
# Added significance threshold lines: Bonferroni threshold obtained from the GWAS
# and nominal significance arbitrarily at 1e-3 (the horizontal lines are for
# sender 225 GWAS and vertical for sender 20)
sp<- ggplot(subsetpvals, aes(x=pkmer20, y=pkmer225)) + geom_point()+
  geom_hline(yintercept =10^-6.551499, linetype="dashed",color = "red")+
  geom_vline(xintercept =10^-6.489483, linetype="dashed",color = "red")+ 
  geom_hline(yintercept =1e-3, linetype="dashed",color = "blue")+
  geom_vline(xintercept =1e-3, linetype="dashed",color = "blue") 
# Transforming to log10
sp + scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10')
# Retrieving the kmers which are significant for both by Bonferroni thresholds
gd3 <- mergepvals$pkmer20<10^-6.489483 & mergepvals$pkmer225 <10^-6.551499
sigkmer_senders <- mergepvals[gd3,]
sigkmer_senders <- na.omit(sigkmer_senders)
# The GWAS of sender 225 is the initial and study 20 is the verifying
# So within Bonferroni of sender 225 and within nominal of sender 20
gd2 <- mergepvals$pkmer20<1e-3 & mergepvals$pkmer225 <10^-6.551499
sigkmer_senders_nom <- mergepvals[gd2,]
sigkmer_senders_nom <- na.omit(sigkmer_senders_nom)
# Obtaining files
write.csv(sigkmer_senders,"signif_kmers_bon_bothstudies.csv")
write.csv(sigkmer_senders_nom,"signif_kmers_nom_bothstudies.csv")


