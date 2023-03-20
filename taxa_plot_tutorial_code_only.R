library(tidyverse)
library(phyloseq)

toy_ps <- read_rds("/Users/lgschaer/Desktop/toy_phyloseq.rds")  
toy_ps

head(sample_data(toy_ps))

genusabundance <- toy_ps %>%
  tax_glom(taxrank = "Genus") %>%                        # Set to smallest taxonomic level you are interested in
  transform_sample_counts(function(x) {x/sum(x)} ) %>%   # Transform to rel. abundance
  psmelt()                                               # Melt to long format
head(genusabundance)

all <- genusabundance %>%
  select(Phylum, Class, Family, Genus, Sample, Abundance, Substrate, Substrate_Label, Replicate) %>%
  filter(Abundance != 0) %>%
  mutate(
    Phylum = as.character(Phylum),
    Class = as.character(Class),
    Family = as.character(Family),
    Genus = as.character(Genus))
head(all)

phylum <- all %>%
  select(Substrate, Substrate_Label, Replicate, Phylum, Abundance) %>%  #choose variables to work with
  group_by(Substrate, Substrate_Label, Replicate) %>%                   #group by variables used to plot NOT taxonomic level
  mutate(totalSum = sum(Abundance)) %>%                                 #calculate total abundance of each Phylum
  ungroup() %>%                                                         #remove grouping variables
  group_by(Substrate, Substrate_Label, Replicate, Phylum) %>%           #group by same variables PLUS taxonomic level
  summarise(                                                            
    Abundance = sum(Abundance),                                         #sum abundance in each phylum for each  group
    totalSum,
    RelAb = Abundance/totalSum) %>%                                     #calculate relative abundance
  unique()

head(phylum)

max(phylum$RelAb)
mean(phylum$RelAb)
min(phylum$RelAb)

length(unique(phylum$Phylum))

phylum_colors <- c(
  "grey22", "darkcyan", "orchid1", "green", "orange", "blue", "tomato2", "olivedrab", "grey47",
  "cyan", "coral3", "darkgreen", "magenta", "palegoldenrod", "dodgerblue", "firebrick", "yellow", "purple4",
  "lightblue", "grey77", "mediumpurple1", "tan4", "red", "darkblue", "yellowgreen")

length(phylum_colors)

ggplot(phylum)+
  geom_col(mapping = aes(x = Replicate, y = RelAb, fill = Phylum), color = "black", position = "stack", show.legend = TRUE)+
  facet_grid(cols = vars(Substrate_Label))+
  ylab("Proportion of Community") +
  xlab(NULL)+
  scale_fill_manual(values = phylum_colors) +
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 1, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 13),
        legend.position = "bottom",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text = element_text(size = 18, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill=guide_legend(ncol=3,byrow=TRUE))


genus <- all %>%
  select(Substrate, Substrate_Label, Replicate, Genus, Abundance) %>%
  group_by(Substrate, Substrate_Label, Replicate) %>%
  mutate(totalSum = sum(Abundance)) %>%
  ungroup() %>%
  group_by(Substrate, Substrate_Label, Replicate, Genus) %>%
  summarise(
    Abundance = sum(Abundance),
    totalSum,
    RelAb = Abundance/totalSum) %>%
  unique()

head(genus)

max(genus$RelAb)
mean(genus$RelAb)
min(genus$RelAb)


length(unique(genus$Genus))

genus <- all %>%
  select(Substrate, Substrate_Label, Replicate, Genus, Abundance) %>%
  group_by(Substrate, Substrate_Label, Replicate) %>%
  mutate(totalSum = sum(Abundance)) %>%
  ungroup() %>%
  group_by(Substrate, Substrate_Label, Replicate, Genus, totalSum) %>%
  summarise(
    Abundance = sum(Abundance),
    Genus = ifelse(Abundance < 0.05, "< 5 %", Genus)) %>%               #change Genus label to group low abundance taxa together
  group_by(Substrate, Substrate_Label, Replicate, Genus, totalSum) %>%  #now group and summarize again to group newly labeled low abundance taxa together
  summarise(
    Abundance = sum(Abundance),
    RelAb = Abundance/totalSum) %>%
  unique()

head(genus)
length(unique(genus$Genus))

colFunc <- colorRampPalette(c("purple4", "lightpink", "firebrick", "orange", "lemonchiffon", "olivedrab4", "darkcyan", "lightblue", "darkblue"))
color_list <- colFunc(length(unique(genus$Genus)))
genus_colors <- c("black", color_list)
length(genus_colors)

ggplot(genus)+
  geom_col(mapping = aes(x = Replicate, y = RelAb, fill = Genus), color = "black", position = "stack", show.legend = TRUE)+
  facet_grid(cols = vars(Substrate_Label))+
  ylab("Proportion of Community") +
  xlab(NULL)+
  scale_fill_manual(values = genus_colors) +
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 1, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 13),
        legend.position = "bottom",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text = element_text(size = 18, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill=guide_legend(ncol=3,byrow=TRUE))








