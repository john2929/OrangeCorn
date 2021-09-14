library(car)
library(qiime2R)
library(ape)
library(dplyr)
library(ggplot2)
library(gplots)
library(lme4)
library(phangorn)
library(phyloseq)
library(plotly)
library(tidyr)
library(vegan)
library(VennDiagram)
library(naniar)
library(zoo)
library(tidyverse)
library(ggpubr)
library(rstatix)

##############################################
#Set UP
#
#These are the things that  we need from Qiime:
#
#sample-metadata.tsv
#taxonomy.tsv
#core-metrics-results-sd-6181/bray_curtis_pcoa_results.qza
#core-metrics-results-sd-6181/weighted_unifrac_pcoa_results.qza
#core-metrics-results-sd-6181/rarefied_table.qza
#rooted-tree.qza
#taxonomy.qza
#core-metrics-results-sd-6181/evenness_vector.qza
#core-metrics-results-sd-6181/faith_pd_vector.qza
#core-metrics-results-sd-6181/observed_otus_vector.qza
#core-metrics-results-sd-6181/shannon_vector.qza
#core-metrics-results-sd-6181/chao1_vector.qza
##############################################

setwd("/Volumes/ag_ansc/Users/john2185/Students/Isaac/data_depot/OrangeCornAgro/DK18_02/WITHOUTMOCK/by_pen/")

qq.line = function(x) {
  # following four lines from base R's qqline()
  y <- quantile(x[!is.na(x)], c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  slope <- diff(y)/diff(x)
  int <- y[1L] - slope * x[1L]
  return(c(int = int, slope = slope))
}

if(!dir.exists("output"))
  dir.create("output")


metadata <- read.delim("sample-metadata-bypen-grouped.txt", sep = "\t", header = T, quote = "", stringsAsFactors = F)
str(metadata)
metadata$DietTreatment <- factor(metadata$DietTreatment, levels = c("WHITE", "YELLOW", "ORANGE"))
#metadata$sample.id <- factor(metadata$sample.id)
row.names(metadata) <- metadata[,1]
#metadata <- metadata[,-1]
metadata$sample.id

#Taxonomy of each OTU
#tax = read.delim("taxonomy.tsv", header=TRUE, sep="\t")

#tax2 = separate(tax, Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")

#This warning means that some cells are empty and that R is replacing empty cells with NA. Now there are other cells that are unclassified that  say, for example `s__` or `g__`. All these need to be removed and replaced with `NA`. 

#All this is OK except that in future use of the taxonomy table, these ASVs will be ignored because they are not classified. Why are ASVs not classified? Its because there is not a close enough match in the database. Just because there is not a good match in the database does not mean they don’t exist, so I wanted to make sure this data was not lost. So in my new code, from lines 200 – 224 I make it so that ASVs that are unclassified at any level are classified as the lowest taxonomic level for which there is a classification.
#Next, all these `NA` classifications with the last level that was classified

##Qiime2r method of reading in the taxonomy files
taxonomy<-read_qza("../taxonomy.qza")
head(taxonomy$data)

tax.clean<-parse_taxonomy(taxonomy$data)
head(tax.clean)

tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("unclassified ", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("unclassified ", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("unclassified ", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("unclassified ", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("unclassified ", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("unclassified ",tax.clean$Genus[i], sep = "_")
  }
}

bc_PCoA<-read_qza("core-metrics-result-sd-6181/bray_curtis_pcoa_results.qza")

corn_colors <- c("Gray", "Yellow", "Orange")
#litter_colors <- c("Black", "Gray")

bc_meta <- bc_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "sample.id"))

#my_column <- "LitterTreatment"
my_column <- "DietTreatment"
centroids <- aggregate(cbind(PC1,PC2)~get(my_column),bc_meta,mean)

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=corn_colors, name = my_column)
ggsave(paste0("output/BC-ellipse_", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches


##SAME thing but with weighted UniFrac

Wuni_PCoA<-read_qza("core-metrics-result-sd-6181//weighted_unifrac_pcoa_results.qza")

corn_colors <- c("Gray", "Yellow", "Orange")
#litter_colors <- c("Black", "Gray")

Wuni_meta <- Wuni_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "sample.id"))

#my_column <- "LitterTreatment"
my_column <- "DietTreatment"
centroids <- aggregate(cbind(PC1,PC2)~get(my_column),Wuni_meta,mean)

ggplot(Wuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Wuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Wuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=corn_colors, name = "Diet Treatment")
ggsave(paste0("output/Wuni-ellipse_", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

##SAME thing but with unweighted UniFrac

uWuni_PCoA<-read_qza("core-metrics-result-sd-6181/unweighted_unifrac_pcoa_results.qza")

uWuni_meta <- uWuni_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "sample.id"))

#my_column <- "LitterTreatment"
my_column <- "DietTreatment"
centroids <- aggregate(cbind(PC1,PC2)~get(my_column),Wuni_meta,mean)

ggplot(uWuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Wuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Wuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=corn_colors, name = "Diet Treatment")
ggsave(paste0("output/uWuni-ellipse_", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches


#################################################################
##Taxa barplot
#################################################################

ASVs <- read_qza("core-metrics-result-sd-6181/rarefied_table.qza")
ASVs$data[1:5,1:5]

physeq <- qza_to_phyloseq(
  features="core-metrics-result-sd-6181/rarefied_table.qza",
  tree="../rooted-tree.qza",
  "../taxonomy.qza",
  metadata = "sample-metadata-bypen-grouped.txt"
)


#First get the OTU table from physeq
physeq_otu_table <- data.frame(otu_table(physeq), check.names = F)
tax.clean = tax.clean[row.names(tax.clean) %in% rownames(physeq_otu_table),]
metadata.filtered = metadata[row.names(metadata) %in% colnames(physeq_otu_table),]

#Assign as variables to be feed into phyloseq
OTU.physeq = otu_table(as.matrix(physeq_otu_table), taxa_are_rows=TRUE)

#our edited and formatted taxonomy table from the top of this script
tax.physeq = tax_table(as.matrix(tax.clean))    
meta.physeq = sample_data(metadata.filtered)

#We then merge these into an object of class phyloseq.

physeq_bar_plot = phyloseq(OTU.physeq, tax.physeq, meta.physeq)



# Set colors for plotting
my_colors <- c(
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black"
)

#If you want different taxonomic level, find and replace the taxonomic level listed here
my_level <- c("Phylum", "Family", "Genus")
my_column <- "DietTreatment"

rm(taxa.summary)

abund_filter <- 0.01
ml ="Family"

for(ml in my_level){
  print(ml)
  
  taxa.summary <- physeq_bar_plot %>%
    tax_glom(taxrank = ml, NArm = FALSE) %>%  # agglomerate at `ml` level
    transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
    psmelt()  %>%                               # Melt to long format
    group_by(get(my_column), get(ml)) %>%
    summarise(Abundance.average=mean(Abundance)) 
  taxa.summary <- as.data.frame(taxa.summary)
  colnames(taxa.summary)[1] <- my_column
  colnames(taxa.summary)[2] <- ml
  
  physeq.taxa.average <- taxa.summary %>% 
    group_by(get(ml)) %>%
    summarise(overall.max=max(Abundance.average))
  
  physeq.taxa.average <- as.data.frame(physeq.taxa.average)
  colnames(physeq.taxa.average)[1] <- ml
  
  # merging the phyla means with the metadata #
  physeq_meta <- merge(taxa.summary, physeq.taxa.average)
  

  physeq_meta_filtered <- filter(physeq_meta, overall.max>abund_filter)
  #str(physeq_meta_filtered)

  #physeq_meta_filtered$body.site.ord = factor(physeq_meta_filtered$body.site, c("left palm", "right palm", "gut", "tongue"))
  
  # Plot 
  ggplot(physeq_meta_filtered, aes(x = get(my_column), y = Abundance.average, fill = get(ml))) + 
    #facet_grid(.~LitterTreatment) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = my_colors) +
    # Remove x axis title
    #theme(axis.title.x = element_blank()) + 
    ylim(c(0,1)) +
    guides(fill = guide_legend(reverse = F, keywidth = .5, keyheight = .5, ncol = 1)) +
    theme(legend.text=element_text(size=8)) +
    #theme(legend.position="bottom") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    theme(legend.title = element_blank()) +
    ylab("Relative Abundance") +
    xlab(my_column) 
    #ggtitle(paste0(ml, " (>", abund_filter * 100,"%) in at least 1 sample")) 
  ggsave(paste0("output/20210716/", ml, "BarPlot_DietOnly.png"), height = 3.5, width = 3.3)
#phylum width = 2.8
#Family width = 3.3
#genus width = 3.5
}


#################################################################
###Differential Abundance with DESeq2
#################################################################


#Adapted from https://joey711.github.io/phyloseq-extensions/DESeq2.html

#First load DESeq2.


library("DESeq2")


#To use DESeq, we need no zeros in our OTU table. So we will edit the table by multiplying by 2 and + 1

#First get the OTU table from physeq
physeq_otu_table <- data.frame(otu_table(physeq), check.names = FALSE)

OTU.clean2 <- physeq_otu_table + 1


#Now make the phyloseq object:
  
  
OTU.physeq = otu_table(as.matrix(OTU.clean2), taxa_are_rows=TRUE)
tax.physeq = tax_table(as.matrix(tax.clean))
meta.physeq = sample_data(metadata.filtered)


#We then merge these into an object of class phyloseq.


physeq_deseq = phyloseq(OTU.physeq, tax.physeq, meta.physeq)


#The following two lines actually do all the complicated DESeq2 work. The function phyloseq_to_deseq2 converts your phyloseq-format microbiome data into a DESeqDataSet with dispersions estimated, using the experimental design formula, also shown (the ~body.site term). The DESeq function does the rest of the testing, in this case with default testing framework, but you can actually use alternatives.


diagdds = phyloseq_to_deseq2(physeq_deseq, ~ DietTreatment)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#the test type of "Wald" tests for significance of coefficients in a Negative Binomial GLM. This is generally a pretty good assumption for sequencing experiments. This was designed with RNA-seq in mind, but also pretty good for 16S sequencing.


###Investigate test results table

#The following results function call creates a table of the results of the tests. Very fast. The hard work was already stored with the rest of the DESeq2-related data in our latest version of the diagdds object (see above). I then order by the adjusted p-value, removing the entries with an NA value. The rest of this example is just formatting the results table with taxonomic information for nice(ish) display in the HTML output.

#Contrast: this argument specifies what comparison to extract from the object to build a results table. There are exactly three elements:
  
#  1. the name of a factor in the design formula, 
#  2. the name of the numerator level for the fold change, and 
#  3. the name of the denominator level for the fold change (simplest case)

alpha = 0.05
my_contrast = c("DietTreatment", "WHITE", "YELLOW")
my_contrast = c("DietTreatment", "WHITE", "ORANGE")
my_contrast = c("DietTreatment", "YELLOW", "ORANGE")
res = results(diagdds, contrast = my_contrast, cooksCutoff = FALSE)

sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq_deseq)[rownames(sigtab), ], "matrix"))
#If you get the following error message after this step, it means you have 
#no significant differences between the groups at the given "alpha" level
###   Error in dimnames(x) <- dn : 
###       length of 'dimnames' [1] not equal to array extent

head(sigtab)


###Volcano Plot

#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", xlim=c(-15,15)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


#Let's look at the OTUs that were significantly different between the two treatment groups. The following makes a nice ggplot2 summary of the results.


theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
scale_fill_brewer(palette = palname, ...)
}
# Phylum order
#x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
#x = sort(x, TRUE)
#sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
DESeq_fig = ggplot(sigtab, aes(x=Genus, y = log2FoldChange, color=Phylum)) + 
  geom_point(size=3) + 
  ylab(paste0("(", my_contrast[2], "/", my_contrast[3], ")\n", "log2FoldChange")) +
  scale_color_manual(values = my_colors[c(4,6)]) +
  #ylim(0,8) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

ggsave(paste0("output/20210716/DESeq2-", my_contrast[2], "-", my_contrast[3], ".png"), DESeq_fig, height = 3, width = 2.25)




#################################################################
### Alpha Diversity
#################################################################


evenness = read_qza("core-metrics-result-sd-6181/evenness_vector.qza")$data %>% rownames_to_column("SampleID")
faith_pd = read_qza("core-metrics-result-sd-6181/faith_pd_vector.qza")$data %>% rownames_to_column("SampleID")
observed_otus = read_qza("core-metrics-result-sd-6181/observed_otus_vector.qza")$data %>% rownames_to_column("SampleID")
shannon = read_qza("core-metrics-result-sd-6181/shannon_vector.qza")$data %>% rownames_to_column("SampleID")

##Merge all alpha diversity measures with metadata

alpha_diversity = merge(x=faith_pd, y=evenness, by.x = "SampleID", by.y = "SampleID")
alpha_diversity = merge(alpha_diversity, observed_otus, by.x = "SampleID", by.y = "SampleID")
alpha_diversity = merge(alpha_diversity, shannon, by.x = "SampleID", by.y = "SampleID")
meta = merge(metadata, alpha_diversity, by.x = 0, by.y = "SampleID")
row.names(meta) = meta$Row.names
meta = meta[,-1]
str(meta)

OTU.clean <- physeq_otu_table 
OTU.clean = OTU.clean[order(row.names(OTU.clean)),]
meta = meta[order(row.names(meta)),]

#Explore the normality of the alpha diversity metrics


#Plots
hist(meta$shannon, main="Shannon diversity", xlab="", breaks=10)
ggqqplot(meta$shannon)
shapiro.test(meta$shannon)
#p>0.05 means normal

hist(meta$faith_pd, main="Faith phylogenetic diversity", xlab="", breaks=10)
ggqqplot(meta$faith_pd)
shapiro.test(meta$faith_pd)
#p>0.05 means normal

hist(meta$pielou_e, main="Evenness", xlab="", breaks=15)
ggqqplot(meta$pielou_e)
shapiro.test(meta$pielou_e)
#p>0.05 means normal

hist(meta$observed_otus, main="Observed OTUs", xlab="", breaks=15)
ggqqplot(meta$observed_otus)
shapiro.test(meta$observed_otus)
#p>0.05 means normal

#To test for normalcy statistically, we ran the Shapiro-Wilk test of normality.
#If the p value from Shapiro-Wilk is > 0.05 we assume normal distribution.

alpha_measures <- c("shannon", "pielou_e", "faith_pd", "observed_otus")
#Run the ANOVA and save it as an object
i <- 3

for(i in 1:length(alpha_measures)){
  print(alpha_measures[i])
  aov.alpha_measures = aov(get(alpha_measures[i]) ~ DietTreatment, data=meta)
  #Call for the summary of that ANOVA, which will include P-values
  print(summary(aov.alpha_measures))
  
  #To do all the pairwise comparisons 
  #between groups and correct for multiple comparisons, 
  #we run Tukey's honest significance test of our ANOVA.
  
  TukeyHSD(aov.alpha_measures)
  
  #We clearly see that the evenness between wet and dry litter are different. 
  
  alpha_plot <- ggplot(meta, aes(DietTreatment, get(alpha_measures[i]))) + 
    geom_boxplot(aes(DietTreatment)) + 
    #facet_grid(.~DietTreatment) +
    #ylim(c(0.7,0.9)) +
    stat_summary(fun=mean, geom="point", shape=3, size=2, color="black", fill="black") +
    ylab(alpha_measures[i]) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  ggsave(paste0("output/", alpha_measures[i], ".png"), alpha_plot, height = 3, width = 3)
}
