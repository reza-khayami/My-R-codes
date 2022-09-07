# Load necessary packages
library(purrr)
library(broom)
library(dplyr)
library(tidyr)

#create example data
## I created this to look like the data you have a screenshot of.
genes.table <- data.frame(gene = paste0("gene", 1:10),
                           mut.cont = runif(10, 0, 10) %>% round(),
                           WT.cont = runif(10, 0, 10) %>% round(),
                           mut.dis = runif(10, 0, 10) %>% round(),
                           WT.dis = runif(10, 0, 10) %>% round())

# Tidy the data
## this makes it easier to convert to 2x2 tables later, which are needed for fisher.test
genes.tidy <- genes.table %>%
  gather(-gene, key = treatment, value = expression) %>% #gathers treatments into one column
  separate(treatment, into = c("genotype", "treatment")) #splits into a genotype and treatment column


# Split by gene
## this creates a list of dataframes
gene.list <- genes.tidy %>%
  split(.$gene)

# Create tables
## This makes 2x2 tables for each gene
gene.tables <- map(gene.list, ~xtabs(expression ~ genotype + treatment, data = .))


# Do fisher test and get p.value
output <- map_dfr(gene.tables,
                  ~fisher.test(.) %>% tidy(), .id = "gene") %>%
  select(gene, p.value) %>% 
  mutate(p.adj = p.adjust(p.value, "fdr"))

output
