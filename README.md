# ADMIXTURE 
<img width="1966" height="372" alt="image" src="https://github.com/user-attachments/assets/92b6c700-f9eb-4d9b-a07a-71a4386c05f7" />

## Introduction
This page is a work in progress!
This repo explains how to run ADMIXTURE on Linux and R from a .vcf file. This pipeline runs in Linux and R and relies on BCFtools (version BCFtools/1.21-GCC-13.3.0), PLINK (version PLINK/2.00a3.7-gfbf-2023a), and VCFtools (version VCFtools/0.1.16-GCC-13.2.0). For this entire pipeline, it may be necessary to split each step up into smaller sample groups. All code should be run in scratch.

Contact: Camille Block (camilleblock@vt.edu)

## Remember to create a SLURM Script
```bash
#!/bin/bash
#SBATCH --nodes=2
#SBATCH --cpus-per-task=20
#SBATCH --time=12:00:00
#SBATCH --job-name STACKS
#SBATCH --account=bedbug
#SBATCH --partition=normal_q
```
## Load needed directories 
```bash
module load BCFtools/1.21-GCC-13.3.0
module load PLINK/2.00a3.7-gfbf-2023a
module load VCFtools/0.1.16-GCC-13.2.0
```
## Subset vcf file so that it runs faster 
Break up samples into groups of 30-50.
```bash
bcftools view -S keep1.txt snps_merged.integer.maf001.geno75.final.vcf.gz -Oz -o subset1.vcf.gz
```
## Index files so that Plink can read them 
This must be done separately for each group of samples.
```bash
tabix -p vcf subset1.vcf.gz
```
## Convert each vcf to bed files (tab delimited file that stores chromosomes as coordinates) and filter 
```bash
plink --bfile snps1 \
     --maf 0.001 \
     --geno 0.75 \
     --allow-extra-chr \
     --make-bed \
     --out snps_subset1.int75
```
## Load ADMIXTURE 
```bash
export PATH=/home/camilleblock/.local/easybuild/software/ADMIXTURE/1.3.0-x86_64:$PATH
echo 'export PATH=/home/camilleblock/.local/easybuild/software/ADMIXTURE/1.3.0-x86_64:$PATH' >> ~/.bashrc
source ~/.bashrc
```
## Run ADMIXTURE for each bed file and for K=1 through K=8 (or more if needed)
```bash
admixture --cv snps_subset1.int75.bed 1 > logsubset1_1.out
```
# Graph ADMIXTURE plots using R 
The K with the lowest CV error is probably correct but you must look at the graphs too. 

## Load necessary libraries 
```bash
library(tidyverse)
library(ggplot2)
```
## Read ADMIXTURE Q file
```bash
Q <- read.table("snps_filtered.int75.5.Q")
K <- ncol(Q)
```
## Add sample names from metadata
```bash
meta <- read.table("pop_map_nodupes.txt", header = TRUE)
````
## Make sure order matches Q file row order
```bash
Q$sample <- meta$sample
colnames(Q)[1:K] <- paste0("Cluster", 1:K)
```
## Merge metadata
```bash
Q <- left_join(Q, meta, by = "sample")
```
## Convert to long format
```bash
Q_long <- Q %>%
  pivot_longer(cols = starts_with("Cluster"),
               names_to = "Cluster",
               values_to = "Ancestry")
```
## Determine major cluster per individual
```bash
Q_long <- Q_long %>%
  group_by(sample) %>%
  mutate(
    MajorCluster = Cluster[which.max(Ancestry)],
    MajorProp = max(Ancestry)
  ) %>%
  ungroup()
```
## Force geographic city order 
```bash
Q_long$city <- factor(
  Q_long$city,
  
  levels = c(
    "Paterson",
    "Hackensack",
    "Newark",
    "Irvington",
    "Jersey_City",
    "Bayonne",
    "Linden",
    "Aberdeen",
    "New_Brunswick",
    "Trenton",
    "Lakewood",
    "Indianpolis"
  )
)
```
## Create hierarchical order
```bash
sample_order <- Q_long %>%
  distinct(sample, city, MajorCluster, MajorProp) %>%
  arrange(city, MajorCluster, desc(MajorProp)) %>%
  pull(sample)

Q_long$sample <- factor(Q_long$sample, levels = sample_order)
```
## Calculate city block sizes
```bash
city_positions <- Q_long %>%
  distinct(sample, city) %>%
  arrange(sample) %>%
  group_by(city) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    start = lag(cumsum(n), default = 0) + 1,
    end = cumsum(n),
    midpoint = (start + end) / 2
  )
```
## Final ADMIXTURE plot
```bash
admixture_plot <- ggplot(Q_long,
                         aes(x = sample,
                             y = Ancestry,
                             fill = Cluster)) +
  scale_fill_manual(values = c("Cluster1" = "#3F858CFF",
                               "Cluster2" = "#6D62AFFF",
                               "Cluster3" = "#F2D43DFF",
                               "Cluster4" = "#D9814EFF",
                               "Cluster5" = "#731A12FF"))+

  geom_bar(stat = "identity", width = 1) +
  
  scale_y_continuous(expand = c(0, 0)) +
  
  scale_x_discrete(expand = c(0, 0)) +
  
  # Vertical separator lines between cities
  geom_vline(data = city_positions[-nrow(city_positions), ],
             aes(xintercept = end + 0.5),
             inherit.aes = FALSE,
             linewidth = 1.0) +
  
  
  # City labels INSIDE the plot, just under bars
  geom_text(data = city_positions,
            aes(x = midpoint, y = -0.001, label = city),
            inherit.aes = FALSE,
            angle = 45,
            hjust = 1,
            size = 3.5,
            fontface = "bold") +
  
  labs(y = "Ancestry proportion",
       title = paste("ADMIXTURE K =", K)) +
  
  coord_cartesian(clip = "off") +  # allow labels outside plot
  
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    plot.margin = margin(10, 10, 60, 10)  # extra bottom space
  )

admixture_plot
```
