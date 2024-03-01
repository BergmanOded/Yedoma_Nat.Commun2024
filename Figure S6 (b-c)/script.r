#  relative stacked area chart

library(ggplot2)
library(dplyr)
library(phyloseq)
library(tidyverse)
library(RColorBrewer)
library(plyr)
# we start with absolute abundance and later convert to relative

feature_table <- read.table(file = "/home/oded/Documents/orit_sivan/raw_data/170-170/cluster100/exclude/alaska/4runs/final/export/absolute-ASV-table/alaska-feature-table.tsv", sep = ",", header = T, row.names = 1, comment.char = "", check.names = FALSE)
str(feature_table)

metadata <- read.table(file = "/home/oded/Documents/orit_sivan/raw_data/170-170/cluster100/exclude/alaska/4runs/final/sample-metadata.tsv", sep = "\t", header = T, 
                       comment.char = "", row.names = 1, check.names = FALSE)

taxonomy <- read.table(file = "/home/oded/Documents/orit_sivan/raw_data/170-170/cluster100/exclude/alaska/4runs/final/export/taxonomy-original-for_r.tsv", sep = "\t", header = T, 
                       comment.char = "", row.names = 1, check.names = FALSE)

#add Kingdom column from the taxonomy, so we can filter to bacteria and archaea

# first change taxonomy rows order to fit feature table rows 
taxonomy_ord <- taxonomy[match(rownames(feature_table), rownames(taxonomy)), ]

# then add Taxon
feature_table_Taxon <- cbind(feature_table, Taxon = taxonomy_ord$Taxon, Kingdom = taxonomy_ord$Kingdom)



# separate to Bacteria and Archaea feature tables
str(feature_table_ord_Taxon)

feature_table_arch<- dplyr::filter(feature_table_ord_Taxon, grepl('Archaea', Kingdom))
feature_table_bact <- dplyr::filter(feature_table_ord_Taxon, grepl('Bacteria', Kingdom))

# now omit the Taxon and Kingdom columns (so we can import to phyloseq)

feature_table_arch_final <- feature_table_arch[c(1:21)]
feature_table_bact_final <- feature_table_bact[c(1:21)]
feature_table_final <- feature_table_ord_Taxon[c(1:21)]

# separate taxonomy to Bacteria and Archaea  
taxonomy_arch<- dplyr::filter(taxonomy_ord, grepl('Archaea', Kingdom))
taxonomy_bact <- dplyr::filter(taxonomy_ord, grepl('Bacteria', Kingdom))


# separate the taxa column into individual columns

#archeaea
taxonomy_arch[c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')] <- str_split_fixed(taxonomy_arch$Taxon, ';', 7)

taxonomy_arch_final <- taxonomy_arch[c(3:9)]

# remove the d__, p__...

taxonomy_arch_final$kingdom<-gsub("d__","",as.character(taxonomy_arch_final$kingdom))
taxonomy_arch_final$phylum<-gsub("p__","",as.character(taxonomy_arch_final$phylum))
taxonomy_arch_final$class<-gsub("c__","",as.character(taxonomy_arch_final$class))
taxonomy_arch_final$order<-gsub("o__","",as.character(taxonomy_arch_final$order))
taxonomy_arch_final$family<-gsub("f__","",as.character(taxonomy_arch_final$family))
taxonomy_arch_final$genus<-gsub("g__","",as.character(taxonomy_arch_final$genus))
taxonomy_arch_final$species<-gsub("s__","",as.character(taxonomy_arch_final$species))

#bacteria
taxonomy_bact[c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')] <- str_split_fixed(taxonomy_bact$Taxon, ';', 7)


taxonomy_bact_final <- taxonomy_bact[c(3:9)]

# remove the d__, p__...

taxonomy_bact_final$kingdom<-gsub("d__","",as.character(taxonomy_bact_final$kingdom))
taxonomy_bact_final$phylum<-gsub("p__","",as.character(taxonomy_bact_final$phylum))
taxonomy_bact_final$class<-gsub("c__","",as.character(taxonomy_bact_final$class))
taxonomy_bact_final$order<-gsub("o__","",as.character(taxonomy_bact_final$order))
taxonomy_bact_final$family<-gsub("f__","",as.character(taxonomy_bact_final$family))
taxonomy_bact_final$genus<-gsub("g__","",as.character(taxonomy_bact_final$genus))
taxonomy_bact_final$species<-gsub("s__","",as.character(taxonomy_bact_final$species))

# bacteria+archaea

taxonomy_ord[c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')] <- str_split_fixed(taxonomy_ord$Taxon, ';', 7)

taxonomy_final <- taxonomy_ord[c(3:9)]

# remove the d__, p__...

taxonomy_final$kingdom<-gsub("d__","",as.character(taxonomy_final$kingdom))
taxonomy_final$phylum<-gsub("p__","",as.character(taxonomy_final$phylum))
taxonomy_final$class<-gsub("c__","",as.character(taxonomy_final$class))
taxonomy_final$order<-gsub("o__","",as.character(taxonomy_final$order))
taxonomy_final$family<-gsub("f__","",as.character(taxonomy_final$family))
taxonomy_final$genus<-gsub("g__","",as.character(taxonomy_final$genus))
taxonomy_final$species<-gsub("s__","",as.character(taxonomy_final$species))


# now create phyloseq objects for bacteria and archaea

ASVs_arch = otu_table(as.matrix(feature_table_arch_final), taxa_are_rows = TRUE)
SAMPLE <- sample_data(metadata)
taxasums_arch <- tax_table(as.matrix(taxonomy_arch_final))
phylo_arch <- phyloseq(ASVs_arch, SAMPLE, taxasums_arch)

ASVs_bact = otu_table(as.matrix(feature_table_bact_final), taxa_are_rows = TRUE)
SAMPLE <- sample_data(metadata)
taxasums_bact <- tax_table(as.matrix(taxonomy_bact_final))
phylo_bact <- phyloseq(ASVs_bact, SAMPLE, taxasums_bact)


ASVs = otu_table(as.matrix(feature_table_final), taxa_are_rows = TRUE)
SAMPLE <- sample_data(metadata)
taxasums_all <- tax_table(as.matrix(taxonomy_final))
phylo_all <- phyloseq(ASVs, SAMPLE, taxasums_all)



#first I will create a table with ASV ID (QIIME2 md5 hashes), counts and taxonomy - bacteria and archeae.
# we will also do this for relative ASVs and both at the genus (qiime2 collapse) level to report

# extract md5 hashes
ASVs_md5_hashes_arch <- tibble::rownames_to_column(feature_table_arch_final, "row_names")
ASVs_md5_hashes_arch <- ASVs_md5_hashes_arch[1]

ASVs_md5_hashes_bact <- tibble::rownames_to_column(feature_table_bact_final, "row_names")
ASVs_md5_hashes_bact <- ASVs_md5_hashes_bact[1]

#extract taxonomy for all ASVs (bact+arch)

ASVs_taxonomy_count_arch <-cbind(ASVs_md5_hashes_arch, feature_table_arch_final, taxonomy_arch_final)
rownames(ASVs_taxonomy_count_arch) <- NULL
colnames(ASVs_taxonomy_count_arch) <- c('ASVs_ID', 'BH1_23cm', 'BH1_56cm', 'BH1_68cm', 'BH1_86cm', 'BH1_106cm', 'BH1_166cm', 'BH1_186cm', 'BH1_259cm', 'BH1_295cm', 'BH1_305cm','BH6_50cm', 'BH6_162cm', 'BH6_198cm', 'BH6_100cm', 'BH6_285cm',	'BH6_345cm', 'BH6_535cm', 'BH6_564cm', 'BH6_655cm', 'BH6_695cm', 'BH6_711cm', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
write.csv(ASVs_taxonomy_count_arch, "/home/oded/Documents/orit_sivan/raw_data/170-170/cluster100/exclude/alaska/4runs/final/export/ASVs_taxonomy_count_arch.csv", row.names=FALSE)

ASVs_taxonomy_count_bact <-cbind(ASVs_md5_hashes_bact, feature_table_bact_final, taxonomy_bact_final)
rownames(ASVs_taxonomy_count_bact) <- NULL
colnames(ASVs_taxonomy_count_bact) <- c('ASVs_ID', 'BH1_23cm', 'BH1_56cm', 'BH1_68cm', 'BH1_86cm', 'BH1_106cm', 'BH1_166cm', 'BH1_186cm', 'BH1_259cm', 'BH1_295cm', 'BH1_305cm','BH6_50cm', 'BH6_162cm', 'BH6_198cm', 'BH6_100cm', 'BH6_285cm',	'BH6_345cm', 'BH6_535cm', 'BH6_564cm', 'BH6_655cm', 'BH6_695cm', 'BH6_711cm', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
write.csv(ASVs_taxonomy_count_bact, "/home/oded/Documents/orit_sivan/raw_data/170-170/cluster100/exclude/alaska/4runs/final/export/ASVs_taxonomy_count_bact.csv", row.names=FALSE)


# now the same table with relative abundance

feature_table_rel <- read.table(file = "/home/oded/Documents/orit_sivan/raw_data/170-170/cluster100/exclude/alaska/4runs/final/export/relative-ASV-table/alaska-relative-table.tsv", sep = "\t", header = T, row.names = 1, comment.char = "", check.names = FALSE)
taxonomy_rel <- taxonomy[match(rownames(feature_table_rel), rownames(taxonomy)), ]
feature_table_rel_Taxon <- cbind(feature_table_rel, Taxon = taxonomy_rel$Taxon, Kingdom = taxonomy_ord_rel$Kingdom)
feature_table_arch_rel <- dplyr::filter(feature_table_rel_Taxon, grepl('Archaea', Kingdom))
feature_table_bact_rel <- dplyr::filter(feature_table_rel_Taxon, grepl('Bacteria', Kingdom))

# now omit the Taxon column 

feature_table_arch_rel_final <- feature_table_arch_rel[c(1:21)]
feature_table_bact_rel_final <- feature_table_bact_rel[c(1:21)]


ASVs_taxonomy_rel_arch <-cbind(ASVs_md5_hashes_arch, feature_table_arch_rel_final, taxonomy_arch_final)
rownames(ASVs_taxonomy_rel_arch) <- NULL
colnames(ASVs_taxonomy_rel_arch) <- c('ASVs_ID', 'BH1_23cm', 'BH1_56cm', 'BH1_68cm', 'BH1_86cm', 'BH1_106cm', 'BH1_166cm', 'BH1_186cm', 'BH1_259cm', 'BH1_295cm', 'BH1_305cm','BH6_50cm', 'BH6_162cm', 'BH6_198cm', 'BH6_100cm', 'BH6_285cm',	'BH6_345cm', 'BH6_535cm', 'BH6_564cm', 'BH6_655cm', 'BH6_695cm', 'BH6_711cm', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
write.csv(ASVs_taxonomy_rel_arch, "/home/oded/Documents/orit_sivan/raw_data/170-170/cluster100/exclude/alaska/4runs/final/export/ASVs_taxonomy_rel_arch.csv", row.names=FALSE)

ASVs_taxonomy_rel_bact <-cbind(ASVs_md5_hashes_bact, feature_table_bact_rel_final, taxonomy_bact_final)
rownames(ASVs_taxonomy_rel_bact) <- NULL
colnames(ASVs_taxonomy_rel_bact) <- c('ASVs_ID', 'BH1_23cm', 'BH1_56cm', 'BH1_68cm', 'BH1_86cm', 'BH1_106cm', 'BH1_166cm', 'BH1_186cm', 'BH1_259cm', 'BH1_295cm', 'BH1_305cm','BH6_50cm', 'BH6_162cm', 'BH6_198cm', 'BH6_100cm', 'BH6_285cm',	'BH6_345cm', 'BH6_535cm', 'BH6_564cm', 'BH6_655cm', 'BH6_695cm', 'BH6_711cm', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
write.csv(ASVs_taxonomy_rel_bact, "/home/oded/Documents/orit_sivan/raw_data/170-170/cluster100/exclude/alaska/4runs/final/export/ASVs_taxonomy_rel_bact.csv", row.names=FALSE)

#Turn all (not filtered to methanotrophs) ASVs into "groupde" genus counts

phylo_bact_genus <- tax_glom(phylo_bact, taxrank = 'genus')
phylo_arch_genus <- tax_glom(phylo_arch, taxrank = 'genus')
phylo_all_genus <- tax_glom(phylo_all, taxrank = 'genus')

#transform to relative abundance
phylo_arch_genus_rel = transform_sample_counts(phylo_arch_genus, function(x) x / sum(x) )
phylo_bact_genus_rel = transform_sample_counts(phylo_bact_genus, function(x) x / sum(x) )
phylo_all_genus_rel = transform_sample_counts(phylo_all_genus, function(x) x / sum(x) )


# filter to a mean threshold 1%

genus_arch_filter_rel = filter_taxa(phylo_arch_genus_rel, function(x) mean(x) > 0.01, TRUE)
genus_bact_filter_rel = filter_taxa(phylo_bact_genus_rel, function(x) mean(x) > 0.01, TRUE)
genus_all_filter_rel = filter_taxa(phylo_all_genus_rel, function(x) mean(x) > 0.01, TRUE)

genus_arch_filter_rel
genus_bact_filter_rel
genus_all_filter_rel


# create dataframe from phyloseq object 
genus_bact_filter_rel_df <- psmelt(genus_bact_filter_rel)
genus_arch_filter_rel_df <- psmelt(genus_arch_filter_rel)
genus_all_filter_rel_df <- psmelt(genus_all_filter_rel)

# omit the ASV column, as it only show the first one. While the genus is multiple ASVs

genus_all_filter_final <- genus_all_filter_rel_df[, c(-1)]

write.csv(genus_all_filter_final, "/home/oded/Documents/orit_sivan/raw_data/170-170/cluster100/exclude/alaska/4runs/final/export/genus_filter_bact_arch.csv", row.names=FALSE)


#Turn all ASVs into "groupde" order counts

phylo_bact_order <- tax_glom(phylo_bact, taxrank = 'order')
phylo_arch_order <- tax_glom(phylo_arch, taxrank = 'order')
phylo_all_order <- tax_glom(phylo_all, taxrank = 'order')

#transform to relative abundance
phylo_arch_order_rel = transform_sample_counts(phylo_arch_order, function(x) x / sum(x) )
phylo_bact_order_rel = transform_sample_counts(phylo_bact_order, function(x) x / sum(x) )
phylo_all_order_rel = transform_sample_counts(phylo_all_order, function(x) x / sum(x) )


# filter to a mean threshold 1%

order_arch_filter_rel = filter_taxa(phylo_arch_order_rel, function(x) mean(x) > 0.01, TRUE)
order_bact_filter_rel = filter_taxa(phylo_bact_order_rel, function(x) mean(x) > 0.01, TRUE)
order_all_filter_rel = filter_taxa(phylo_all_order_rel, function(x) mean(x) > 0.01, TRUE)

order_arch_filter_rel
order_bact_filter_rel
order_all_filter_rel


# create dataframe from phyloseq object 
order_bact_filter_rel_df <- psmelt(order_bact_filter_rel)
order_arch_filter_rel_df <- psmelt(order_arch_filter_rel)
order_all_filter_rel_df <- psmelt(order_all_filter_rel)

# omit the ASV column, as it only shows the first one. While the order is multiple ASVs

order_all_filter_final <- order_all_filter_rel_df[, c(-1)]

write.csv(order_all_filter_final, "/home/oded/Documents/orit_sivan/raw_data/170-170/cluster100/exclude/alaska/4runs/final/export/order_filter_bact_arch.csv", row.names=FALSE)


# make sure genus is character

str(genus_bact_filter_rel_df)
str(genus_arch_filter_rel_df)
str(genus_all_filter_rel_df)
# if not: df$phylum <- as.character(df$phylum) #convert to character

str(order_bact_filter_rel_df)
str(order_arch_filter_rel_df)
str(order_all_filter_rel_df)

# creat new column that combines phylum and order

order_arch_filter_rel_df$phylum_order <- paste(order_arch_filter_rel_df$phylum, order_arch_filter_rel_df$order)
order_bact_filter_rel_df$phylum_order <- paste(order_bact_filter_rel_df$phylum, order_bact_filter_rel_df$order)

# stack plot

#first Archaea
# create color pallet

nb.cols_arch_filter <- 17
mycolors_arch_filter <- colorRampPalette(brewer.pal(12, 'Paired'))(nb.cols_arch_filter)

# we will sort levels by value for depth (Sample)
order_arch_filter_rel_df$depth <- factor(order_arch_filter_rel_df$depth,levels = c("305", "295", "259", "186", "166", "106", "86", "68", "56", "23", "711", "695", "655", "564", "535", "345", "285", "198", "162", "100", "50")) 

str(order_arch_filter_rel_df)



#to plot subset the df by core 1 and 2
order_arch_filter_BH1 <- order_arch_filter_rel_df[order_arch_filter_rel_df$core == "BH1", ]
order_arch_filter_BH6 <- order_arch_filter_rel_df[order_arch_filter_rel_df$core == "BH6", ]


p_arch_BH1 <- ggplot(order_arch_filter_BH1, aes(x=depth, y=Abundance, fill=phylum_order, group=phylum_order))
p1_arch_BH1 <- p_arch_BH1 + geom_area(position = "fill", colour = "black", linewidth = .5, alpha = .7) +
  scale_fill_manual(values = mycolors_arch_filter) + 
    labs(y=expression(paste("Relative abundance (proportion)")), x=expression(paste("Depth (cm)"))) +
  theme(legend.position='right') + theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank())+
  theme(legend.text=element_text(size=10))+
  theme(axis.title.y = element_text(size=12, color = "black"), axis.text.y = element_text(size=12, color = "black")) + 
  theme(axis.title.x = element_text(size=12, color = "black"), axis.text.x = element_text(size=12, color = "black")) +  
  theme(legend.title = element_text(size = 10))+theme(legend.background = element_rect(fill  = "transparent"),  legend.box.background = element_rect(linetype=2, colour="gray"))+
  theme(legend.position = "right")
  p1_arch_BH1+coord_flip() 

  
  p_arch_BH6 <- ggplot(order_arch_filter_BH6, aes(x=depth, y=Abundance, fill=phylum_order, group=phylum_order))
  p1_arch_BH6 <- p_arch_BH6 + geom_area(position = "fill", colour = "black", linewidth = .5, alpha = .7) +
    scale_fill_manual(values = mycolors_arch_filter) + 
    labs(y=expression(paste("Relative abundance (proportion)")), x=expression(paste("Depth (cm)"))) +
    theme(legend.position='right') + theme_classic() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(panel.border = element_blank())+
    theme(legend.text=element_text(size=10))+
    theme(axis.title.y = element_text(size=12, color = "black"), axis.text.y = element_text(size=12, color = "black")) + 
    theme(axis.title.x = element_text(size=12, color = "black"), axis.text.x = element_text(size=12, color = "black")) +  
    theme(legend.title = element_text(size = 10))+theme(legend.background = element_rect(fill  = "transparent"),  legend.box.background = element_rect(linetype=2, colour="gray"))+
    theme(legend.position = "right")
  p1_arch_BH6+coord_flip()  
  
  
  pdf("/home/oded/Documents/orit_sivan/raw_data/170-170/cluster100/exclude/alaska/4runs/final/export/plot_area/23.08.16_alaska_stack_archea_plot_BH1.pdf",  width = 10, height = 6)
  p1_arch_BH1+coord_flip() 
  dev.off()
  
  pdf("/home/oded/Documents/orit_sivan/raw_data/170-170/cluster100/exclude/alaska/4runs/final/export/plot_area/23.08.16_alaska_stack_archea_plot_BH6.pdf",  width = 10, height = 6)
  p1_arch_BH6+coord_flip() 
  dev.off()
  

  

#Now Bacteria
# we have many phyla for Bacteria, so we will convert to others (<1% abundance)

order_bact_filter_rel_df$phylum_order[order_bact_filter_rel_df$Abundance < 0.01] <- "others"

# so after conversion we have several "others" phyla for the same depth, these need to be summed

order_bact_filter_rel_df_sumed <- order_bact_filter_rel_df %>% group_by(Sample, depth, core, phylum_order) %>% summarize_at("Abundance", sum) 

# many Phyla do not appear in all Samples, so they must be add with Abundance 0

area_plot <- order_bact_filter_rel_df_sumed %>%
  # Ungroup or `complete` won't work as expected
  ungroup() %>%
  # Add all (Type, year) combinations, filling in with 0s where `n()` is not observed
  complete(phylum_order, Sample, fill = list(`n()` = 0))
#now replace NA with 0
area_plot$Abundance[is.na(area_plot$Abundance)] <- 0

# now the depth and core columns should also be completed 

area_plot[area_plot$Sample == '68_BH1'& is.na(area_plot$depth), "depth"] <- 68
area_plot[area_plot$Sample == '86_BH1'& is.na(area_plot$depth), "depth"] <- 86
area_plot[area_plot$Sample == '106_BH1'& is.na(area_plot$depth), "depth"] <- 106
area_plot[area_plot$Sample == '166_BH1'& is.na(area_plot$depth), "depth"] <- 166
area_plot[area_plot$Sample == '186_BH1'& is.na(area_plot$depth), "depth"] <- 186
area_plot[area_plot$Sample == '259_BH1'& is.na(area_plot$depth), "depth"] <- 259
area_plot[area_plot$Sample == '23_BH1'& is.na(area_plot$depth), "depth"] <- 23
area_plot[area_plot$Sample == '56_BH1'& is.na(area_plot$depth), "depth"] <- 56
area_plot[area_plot$Sample == '295_BH1'& is.na(area_plot$depth), "depth"] <- 295
area_plot[area_plot$Sample == '305_BH1'& is.na(area_plot$depth), "depth"] <- 305
area_plot[area_plot$Sample == '50_BH6'& is.na(area_plot$depth), "depth"] <- 50
area_plot[area_plot$Sample == '162_BH6'& is.na(area_plot$depth), "depth"] <- 162
area_plot[area_plot$Sample == '198_BH6'& is.na(area_plot$depth), "depth"] <- 198
area_plot[area_plot$Sample == '100_BH6'& is.na(area_plot$depth), "depth"] <- 100
area_plot[area_plot$Sample == '285_BH6'& is.na(area_plot$depth), "depth"] <- 285
area_plot[area_plot$Sample == '345_BH6'& is.na(area_plot$depth), "depth"] <- 345
area_plot[area_plot$Sample == '535_BH6'& is.na(area_plot$depth), "depth"] <- 535
area_plot[area_plot$Sample == '564_BH6'& is.na(area_plot$depth), "depth"] <- 564
area_plot[area_plot$Sample == '655_BH6'& is.na(area_plot$depth), "depth"] <- 655
area_plot[area_plot$Sample == '695_BH6'& is.na(area_plot$depth), "depth"] <- 695
area_plot[area_plot$Sample == '711_BH6'& is.na(area_plot$depth), "depth"] <- 711


area_plot[area_plot$Sample == '68_BH1'& is.na(area_plot$core), "core"] <- 'BH1'
area_plot[area_plot$Sample == '86_BH1'& is.na(area_plot$core), "core"] <- 'BH1'
area_plot[area_plot$Sample == '106_BH1'& is.na(area_plot$core), "core"] <- 'BH1'
area_plot[area_plot$Sample == '166_BH1'& is.na(area_plot$core), "core"] <- 'BH1'
area_plot[area_plot$Sample == '186_BH1'& is.na(area_plot$core), "core"] <- 'BH1'
area_plot[area_plot$Sample == '259_BH1'& is.na(area_plot$core), "core"] <- 'BH1'
area_plot[area_plot$Sample == '23_BH1'& is.na(area_plot$core), "core"] <- 'BH1'
area_plot[area_plot$Sample == '56_BH1'& is.na(area_plot$core), "core"] <- 'BH1'
area_plot[area_plot$Sample == '295_BH1'& is.na(area_plot$core), "core"] <- 'BH1'
area_plot[area_plot$Sample == '305_BH1'& is.na(area_plot$core), "core"] <- 'BH1'
area_plot[area_plot$Sample == '50_BH6'& is.na(area_plot$core), "core"] <- 'BH6'
area_plot[area_plot$Sample == '162_BH6'& is.na(area_plot$core), "core"] <- 'BH6'
area_plot[area_plot$Sample == '198_BH6'& is.na(area_plot$core), "core"] <- 'BH6'
area_plot[area_plot$Sample == '100_BH6'& is.na(area_plot$core), "core"] <- 'BH6'
area_plot[area_plot$Sample == '285_BH6'& is.na(area_plot$core), "core"] <- 'BH6'
area_plot[area_plot$Sample == '345_BH6'& is.na(area_plot$core), "core"] <- 'BH6'
area_plot[area_plot$Sample == '535_BH6'& is.na(area_plot$core), "core"] <- 'BH6'
area_plot[area_plot$Sample == '564_BH6'& is.na(area_plot$core), "core"] <- 'BH6'
area_plot[area_plot$Sample == '655_BH6'& is.na(area_plot$core), "core"] <- 'BH6'
area_plot[area_plot$Sample == '695_BH6'& is.na(area_plot$core), "core"] <- 'BH6'
area_plot[area_plot$Sample == '711_BH6'& is.na(area_plot$core), "core"] <- 'BH6'
# assign colors

nb.cols_bact_filter <- 24
mycolors_bact_filter <- colorRampPalette(brewer.pal(12, 'Paired'))(nb.cols_bact_filter)

# we will sort levels by value for depth (Sample)
area_plot$depth <- factor(area_plot$depth,levels = c("305", "295", "259", "186", "166", "106", "86", "68", "56", "23", "711", "695", "655", "564", "535", "345", "285", "198", "162", "100", "50")) 

str(area_plot)


#to plot subset the df by core 1 and 2
order_bact_filter_BH1 <- area_plot[area_plot$core == "BH1", ]
order_bact_filter_BH6 <- area_plot[area_plot$core == "BH6", ]



p_bact_BH1 <- ggplot(order_bact_filter_BH1, aes(x=depth, y=Abundance, fill=phylum_order, group=phylum_order))
p1_bact_BH1 <- p_bact_BH1 + geom_area(position = "fill", colour = "black", linewidth = .5, alpha = .7) +
  scale_fill_manual(values = mycolors_bact_filter) + 
  labs(y=expression(paste("Relative abundance (proportion)")), x=expression(paste("Depth (cm)"))) +
  theme(legend.position='right') + theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank())+
  theme(legend.text=element_text(size=10))+
  theme(axis.title.y = element_text(size=12, color = "black"), axis.text.y = element_text(size=12, color = "black")) + 
  theme(axis.title.x = element_text(size=12, color = "black"), axis.text.x = element_text(size=12, color = "black")) +  
  theme(legend.title = element_text(size = 10))+theme(legend.background = element_rect(fill  = "transparent"),  legend.box.background = element_rect(linetype=2, colour="gray"))+
  theme(legend.position = "right")
p1_bact_BH1+coord_flip() 

pdf("/home/oded/Documents/orit_sivan/raw_data/170-170/cluster100/exclude/alaska/4runs/final/export/plot_area/23.08.16_alaska_stack_bacteria_plot_BH1.pdf",  width = 10, height = 6, family="ArialMT")
p1_bact_BH1+coord_flip()
dev.off() 


pdf("/home/oded/Documents/orit_sivan/raw_data/170-170/cluster100/exclude/alaska/4runs/final/export/plot_area/23.08.16_alaska_stack_bacteria_plot_BH2.pdf",  width = 10, height = 6, family="ArialMT")
p1_bact_BH2+coord_flip()
dev.off() 



p_bact_BH6 <- ggplot(order_bact_filter_BH6, aes(x=depth, y=Abundance, fill=phylum_order, group=phylum_order))
p1_bact_BH6 <- p_bact_BH6 + geom_area(position = "fill", colour = "black", size = .5, alpha = .7) +
  scale_fill_manual(values = mycolors_bact_filter) + 
  labs(y=expression(paste("Relative abundance (proportion)")), x=expression(paste("Depth (cm)"))) +
  theme(legend.position='right') + theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank())+
  theme(legend.text=element_text(size=10))+
  theme(axis.title.y = element_text(size=12, color = "black"), axis.text.y = element_text(size=12, color = "black")) + 
  theme(axis.title.x = element_text(size=12, color = "black"), axis.text.x = element_text(size=12, color = "black")) +  
  theme(legend.title = element_text(size = 10))+theme(legend.background = element_rect(fill  = "transparent"),  legend.box.background = element_rect(linetype=2, colour="gray"))+
  theme(legend.position = "right")
p1_bact_BH6+coord_flip() 

pdf("/home/oded/Documents/orit_sivan/raw_data/170-170/cluster100/exclude/alaska/4runs/final/export/plot_area/23.08.16_alaska_stack_bacteria_plot_BH6.pdf",  width = 10, height = 6, family="ArialMT")
p1_bact_BH6+coord_flip()
dev.off() 
