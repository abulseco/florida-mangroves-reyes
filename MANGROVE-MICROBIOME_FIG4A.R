# Flow-through Plot 1
# PCoA and statistics
# Florida mangrove collab

# SETUP----
# This code assumes you have run phyloseq successfully ('MANGROVE-MICROBIOME-PHYLO.R')
# Uses the Phyloseq object 'phy.f.nolow.feb'


## Load libraries----
library(phyloseq); library(data.table); library(ggplot2); library(dplyr); library(viridis);library(vegan)

## Graph formatting function----
pretty.theme <- function(){
  theme_bw() +
    theme(axis.text.x=element_text(size = 18, color="black"),
          axis.text.y=element_text(size = 18, color="black"),
          # panel.grid.major = element_blank(),
          # panel.grid.minor = element_blank(),
          text=element_text(size=18),
          panel.background = element_blank(),
          panel.border = element_rect(fill=NA),
          axis.line = element_line(colour = "black"))
}

# PLOTTING ORDINATION----
# Directions for plotting
# "#DDCC77" for the "Mixed-Basin"
# "#CC6677" for the "Red-Fringe"
# Add the ecotype (i.e., basin and fringe)

## Ordination--------
BC_distance <- phyloseq::distance(phy.f.nolow.feb, "bray")
bcOrd <- ordinate(phy.f.nolow.feb, "PCoA", BC_distance)
plot_scree(bcOrd)

p1 <- plot_ordination(phy.f.nolow.feb, bcOrd, color = "Ecotype", shape = "Treatment") +
  geom_point(size = 6, alpha = 0.9) +
  stat_ellipse(aes(group = Ecotype), type = "t", level = 0.95) +
  # geom_text(aes(label = Sample_ID), vjust = 1, size = 5, nudge_x = 0.025) + 
  scale_color_manual(values = c("#DDCC77", "#CC6677")) +
  pretty.theme() +
  labs(color = "Stand Location", shape = "Treatment") +
  theme(legend.position = "left")
p1

## Permanova----
set.seed(1)

feb_bray <- phyloseq::distance(phy.f.nolow.feb, method = "bray")
sampledf <- data.frame(sample_data(phy.f.nolow.feb))
adonis2(feb_bray ~ Ecotype, data = sampledf)
adonis2(feb_bray ~ Treatment, data = sampledf)

# Homogeneity of dispersion test
beta <- betadisper(feb_bray, sampledf$Ecotype, type = "median")
permutest(beta)
