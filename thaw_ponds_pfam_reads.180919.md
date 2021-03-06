---
title: "Thaw ponds: differential abundance analysis for genes involved in carbon degradation"
author: "Domenico Simone"
date: "November 1st, 2018"
#output: html_document
output:
  html_document:
    keep_md: true
    code_folding: hide
    toc: true
#    toc_float:
#      collapsed: false
#      smooth_scroll: true
---



## Summary

Differential abundance analysis of genes involved in carbon degradation in samples from thaw ponds, according to:

- Ages in epilimnion samples
- Layers in old ponds
- Oxygen content in old ponds.

The analyses were performed with two different strategies:

- Quantifying reads directly mapping on PFAM entries of interest;
- Quantifying reads mapping on genes annotated on contigs from a coassembly of all read datasets.

## Methods in brief

Raw reads were trimmed with **sickle**v1.210 (default parameters) and mapped on Pfam HMMs of interest with **metadomain**. Read datasets were splitted in subsets of 500K reads for this analysis and the results combined in this protocol. 
Contigs co-assembled from all read datasets were annotated against PFAM using **hmmer**. Read datasets were mapped against the contigs with **bowtie2**. Reads mapping on annotated PFAM entries were counted with **htseq-count**. 
Two separate data exploration and DA analysis were performed, by using read mappings on PFAM entries of interest and read mappings on contigs assembled from metagenomic datasets.
Data exploration was performed with the **metaMDS** function available through the R **vegan** package, to perform NMDS to compare the distribution of PFAM entries of interest in the analysed samples. Differential abundance (DA) analysis was performed with R package **DESeq2**. DA analyses were carried out by performing pairwise comparisons of tested multiple conditions, namely *old*, *emerge*, *medium* within epilimnion samples, and *epilimnion*, *metalimnion* and *hypolimnion* within old samples. An additional DA analysis was performed considering altogether *epilimnion* and *metalimnion* as *oxic* samples, compared to *hypolimnion* samples which were considered as *anoxic* samples. For each DA analysis, two outputs are shown: 1) MA-plot (described in the next section) and 2) a searchable/filterable table with DA results reported for each PFAM entry. Each table is also exported as csv file.

### Note on MA plots

The MA-plot shows the log2 fold changes over the mean of normalized counts, i.e. the average of counts normalized by size factor. The x axis is the average abundance over all samples, the y axis the log2 fold change between the first condition and the second condition *eg* if an MA-plot shows epilimnion vs metalimnion, PFAM entries with a positive log2 fold change are more abundant in epilimnion samples, when compared to metalimnion samples. For each PFAM entry an adjusted p-value is calculated by the DA pipeline. PFAM entries with an adjusted p-value below a threshold (here 0.1, the suggested default) are shown in red.

### Functions used in this report (click "Code" on the right to expand)


```r
# code to read all metadomain results
getBigTable <- function(i){
  small.table <- read.table(i, header = FALSE, stringsAsFactors = FALSE)
  colnames(small.table) <- c("PFAM_ID", "PFAM_length", "n_reads", "PFAM_cov", "PFAM_status")
  # i is eg ../metadomain_1M/A-1_S21_L005.0.1.fasta.metadomain
  sample_name <- strsplit(basename(i), "_")[[1]][1]
  p <- data.frame(sample=c(sample_name), stringsAsFactors = FALSE)
  small.table <- cbind(p, small.table)
}

getBigTable.general <- function(i, ColNames=NULL){
  small.table <- read.table(i, header = FALSE, stringsAsFactors = FALSE)
  colnames(small.table) <- ColNames
  sample_name <- strsplit(basename(i), "_")[[1]][1]
  p <- data.frame(sample=c(sample_name), stringsAsFactors = FALSE)
  small.table <- cbind(p, small.table)
}

# Calculate TPM for each PFAM in each sample
#calculate.RPK <- function(table.counts, sample, n_reads, feature.length){
calculate.RPK <- function(table.counts, sample, feature.name, pfam.data){
  RPK <- table.counts %>%
            base::subset(sample == sample)
  return(RPK)
}

# this is for read counts from mapping on contigs which were annotated against PFAM
getBigTable.tpm.pfam <- function(i, pfam.lengths.nt){
  small.table <- read.table(i, header = FALSE, stringsAsFactors = FALSE, col.names = c("PFAM_ID", "n_reads")) %>%
    filter(!base::grepl("^_", PFAM_ID)) %>%
    mutate_at("PFAM_ID", str_trunc, width=7, ellipsis="") %>%
#strsplit(pfam.lengths$PFAM_ID, split = "\\.") %>% unlist %>% matrix(ncol = 2, byrow = TRUE) %>% as.data.frame %>% .$V1
    left_join(pfam.lengths, by = c("PFAM_ID" = "PFAM_ID.noversion")) %>%
    select(-c(PFAM_ID.y, PFAM_length)) %>%
    mutate(tpm=n_reads/PFAM_length_nt/sum(n_reads/PFAM_length_nt)*1e06)
    #big.table.pfam[base::grep("^_", big.table.pfam$PFAM_ID, invert=TRUE),]
  #colnames(small.table) <- c("PFAM_ID", "n_reads")
  # i is eg ../metadomain_1M/A-1_S21_L005.0.1.fasta.metadomain
  sample_name <- strsplit(basename(i), "_")[[1]][1]
  p <- data.frame(sample=c(sample_name), stringsAsFactors = FALSE)
  small.table <- cbind(p, small.table)
}

plotMA.significant <- function(df, t=0.1, cex=0.8, main="DE analysis - MA plot"){
  plotMA(df, main=main)
  xx <- df$baseMean
  yy <- df$log2FoldChange
  # Replace NAs in padj with values that you'll never want to consider :-)
  signSubset <- replace(df$padj, is.na(df$padj), 1000) < t
  if(count(signSubset, TRUE) > 0)
    text( xx[signSubset],
          yy[signSubset],
          labels=rownames(df)[signSubset],
          pos=1, cex=cex)
    text( xx[c(TRUE)],
          yy[c(TRUE)],
          labels=c(""))
}
# # pointLabel(x, y, as.character(round(x,5)), offset = 0, cex = .7)
# 
# plotMA.significant.pointlabel <- function(df, t=0.1, cex=0.8, main="DE analysis - MA plot"){
#   plotMA(df, main=main)
#   xx <- df$baseMean
#   yy <- df$log2FoldChange
#   # Replace NAs in padj with values that you'll never want to consider :-)
#   signSubset <- replace(df$padj, is.na(df$padj), 1000) < t
#   if(count(signSubset, TRUE) > 0)
#     pointLabel(xx[signSubset],
#                yy[signSubset],
#                labels=rownames(df)[signSubset],
#                pos=1, cex=cex)
#     # text( xx[signSubset],
#     #       yy[signSubset],
#     #       labels=rownames(df)[signSubset],
#     #       pos=1, cex=cex)
#     text( xx[c(TRUE)],
#           yy[c(TRUE)],
#           labels=c(""))
# }
# 
# plotMA.significant.thigmo <- function(df, t=0.1, cex=0.8, main="DE analysis - MA plot"){
#   plotMA(df, main=main)
#   xx <- df$baseMean
#   yy <- df$log2FoldChange
#   # Replace NAs in padj with values that you'll never want to consider :-)
#   signSubset <- replace(df$padj, is.na(df$padj), 1000) < t
#   if(count(signSubset, TRUE) > 0)
#     thigmophobe.labels(xx[signSubset],
#                yy[signSubset],
#                labels=rownames(df)[signSubset],
#                cex=cex)
#     # pointLabel(xx[signSubset],
#     #            yy[signSubset],
#     #            labels=rownames(df)[signSubset],
#     #            pos=1, cex=cex)
#     # text( xx[signSubset],
#     #       yy[signSubset],
#     #       labels=rownames(df)[signSubset],
#     #       pos=1, cex=cex)
#     text( xx[c(TRUE)],
#           yy[c(TRUE)],
#           labels=c(""))
# }

# since the DESeq2 table has PFAM ids as rownames,
# we need to extract the PFAM ids and add them to the dataframe
# as the first column
DESeq2datatable <- function(DESeq2.data){
  DESeq2.data <- DESeq2.data[which(DESeq2.data$padj < 0.1),]
  datatable(data.frame(cbind(data.frame(PFAM_ID=row.names(DESeq2.data)), DESeq2.data)), rownames = FALSE, filter = "top")
}

# function to write out table DESeq2 table in csv format
DESeq2csv <- function(DESeq2.data, file="DESeq2_data.csv"){
  DESeq2.data <- DESeq2.data[which(DESeq2.data$padj < 0.1),]
  write.table(data.frame(cbind(data.frame(PFAM_ID=row.names(DESeq2.data)), DESeq2.data)), row.names = FALSE, file=file, sep = ",")
}
```

## Mapping of reads on HMM profiles

### Data preparation


```r

## Create results dir
dir.create("results", showWarnings = FALSE)

# Read sample metadata and join it with table with dataset size:

metadata <- read_csv("data/metadata.csv")
## Warning: Missing column names filled in: 'X24' [24]
## Parsed with column specification:
## cols(
##   .default = col_double(),
##   sample_code = col_character(),
##   Pond = col_character(),
##   rawdata = col_character(),
##   size = col_character(),
##   layer = col_character(),
##   NH4 = col_integer(),
##   PO4 = col_integer(),
##   X24 = col_character(),
##   methanotrophs = col_character(),
##   sample = col_character()
## )
## See spec(...) for full column specifications.
metadata$sample_name <- gsub("_", "-", metadata$sample_code)
colnames(metadata)[(names(metadata) == "size")] <- "age"
rownames(metadata) <- metadata$sample_name
## Warning: Setting row names on a tibble is deprecated.

metadata.dataset_size <- read_csv("data/sample_name_reads")
## Parsed with column specification:
## cols(
##   sample_name = col_character(),
##   total_reads = col_integer()
## )

# add dataset sizes
metadata <- metadata %>%
  left_join(metadata.dataset_size)
## Joining, by = "sample_name"

# For each sample, read all result files

all_files = list.files("./metadomain_1M", full.names = TRUE)
big.table <- data.frame(rbindlist(lapply(all_files, function(x) getBigTable(paste(x)))))
big.table$n_reads <- as.integer(big.table$n_reads)
big.table$PFAM_length <- as.integer(big.table$PFAM_length)
big.table$PFAM_cov <- as.numeric(big.table$PFAM_cov)

# Group counts by sample and PFAM_ID, **save it**.
all.final.table <- big.table %>% 
  group_by(sample, PFAM_ID, PFAM_length) %>%
  summarise(sum(n_reads)) %>%
  ungroup()

colnames(all.final.table)[4] <- "read_count"

# Convert count table in matrix
all.final.matrix <- all.final.table %>%
  subset(select = -PFAM_length) %>%
  spread(PFAM_ID, read_count)

rownames.all.final.matrix <- all.final.matrix$sample
all.final.matrix <- data.matrix(subset(all.final.matrix, select = -sample)) 
rownames(all.final.matrix) <- rownames.all.final.matrix
```

### Overview of read mappings on PFAM HMMs of interest

Reported in file `results/read_mapping_carbon_HMM.csv`.


```r
data.frame(sample_name=row.names(all.final.matrix), total_mapped_reads=rowSums(all.final.matrix)) %>%
  left_join(metadata.dataset_size) %>%
  mutate("total_mapped_reads (%)"=total_mapped_reads*100/total_reads) %>%
  write.table(file="results/read_mapping_carbon_HMM.csv", quote=FALSE, sep="\t", row.names = FALSE)
```

```
## Joining, by = "sample_name"
```

```
## Warning: Column `sample_name` joining factor and character vector, coercing
## into character vector
```

```r
  #datatable(rownames = FALSE, filter = "top", options = list(pageLength = 30))
```

### NMDS on normalized data

Perform NMDS on normalized data. Data are normalized by rarefying to the sample with the lowest number of reads mapped on PFAM HMMs of interest and plotted with ggplot.


```r
# source: http://geoffreyzahn.com/nmds-example/
min_depth = min(rowSums(all.final.matrix))

# normalize data
all.final.matrix.rarefied <- as.matrix(round(rrarefy(all.final.matrix, min_depth)))

# calculate distance matrix
sample_dist.rarefied <- as.matrix((vegdist(all.final.matrix.rarefied, "bray")))

# perform actual NMDS
NMDS.rarefied <- metaMDS(sample_dist.rarefied)
## Run 0 stress 0.09887482 
## Run 1 stress 0.1882314 
## Run 2 stress 0.09887483 
## ... Procrustes: rmse 4.149699e-05  max resid 0.0001221684 
## ... Similar to previous best
## Run 3 stress 0.09887482 
## ... Procrustes: rmse 6.154622e-06  max resid 1.749881e-05 
## ... Similar to previous best
## Run 4 stress 0.09887482 
## ... New best solution
## ... Procrustes: rmse 8.73508e-06  max resid 3.214115e-05 
## ... Similar to previous best
## Run 5 stress 0.09887528 
## ... Procrustes: rmse 7.048768e-05  max resid 0.0002191241 
## ... Similar to previous best
## Run 6 stress 0.09887482 
## ... Procrustes: rmse 1.697827e-05  max resid 5.084294e-05 
## ... Similar to previous best
## Run 7 stress 0.1586289 
## Run 8 stress 0.09887482 
## ... Procrustes: rmse 1.681319e-06  max resid 4.969649e-06 
## ... Similar to previous best
## Run 9 stress 0.09887482 
## ... Procrustes: rmse 8.487312e-06  max resid 2.902382e-05 
## ... Similar to previous best
## Run 10 stress 0.09887482 
## ... Procrustes: rmse 3.76365e-06  max resid 1.688253e-05 
## ... Similar to previous best
## Run 11 stress 0.09887482 
## ... New best solution
## ... Procrustes: rmse 1.635403e-06  max resid 5.186226e-06 
## ... Similar to previous best
## Run 12 stress 0.09887483 
## ... Procrustes: rmse 3.100322e-05  max resid 0.000106676 
## ... Similar to previous best
## Run 13 stress 0.09887482 
## ... Procrustes: rmse 2.993185e-06  max resid 7.932806e-06 
## ... Similar to previous best
## Run 14 stress 0.09887488 
## ... Procrustes: rmse 8.190714e-05  max resid 0.0002684636 
## ... Similar to previous best
## Run 15 stress 0.09887482 
## ... Procrustes: rmse 4.091593e-06  max resid 1.234773e-05 
## ... Similar to previous best
## Run 16 stress 0.09887482 
## ... Procrustes: rmse 2.291966e-05  max resid 8.326452e-05 
## ... Similar to previous best
## Run 17 stress 0.09887482 
## ... Procrustes: rmse 2.461911e-06  max resid 5.782753e-06 
## ... Similar to previous best
## Run 18 stress 0.09887494 
## ... Procrustes: rmse 0.0001085749  max resid 0.0003572379 
## ... Similar to previous best
## Run 19 stress 0.09887486 
## ... Procrustes: rmse 5.529639e-05  max resid 0.0001721149 
## ... Similar to previous best
## Run 20 stress 0.09887486 
## ... Procrustes: rmse 6.578666e-05  max resid 0.0002044114 
## ... Similar to previous best
## *** Solution reached

# build a data frame with NMDS coordinates and metadata
#MDS1 = NMDS$points[,1]
#MDS2 = NMDS$points[,2]
NMDS.rarefied.df = data.frame(MDS1 = NMDS.rarefied$points[,1],
                              MDS2 = NMDS.rarefied$points[,2])
NMDS.rarefied.df$sample_name <- rownames(NMDS.rarefied.df)
NMDS.rarefied.df <- NMDS.rarefied.df %>%
  left_join(metadata)
## Joining, by = "sample_name"

# # Plot with ggplot
# ggplot(NMDS.rarefied.df, aes(x=MDS1, y=MDS2, col=age, shape=layer)) +
#  geom_point() +
#  geom_text(aes(label=sample_name),hjust=-0.15,vjust=0, size=3) +
#  #stat_ellipse() +
#  theme_bw() +
#  labs(title = "NMDS of read mapping on PFAM entries\nrelated to carbon degradation (normalized)")

# Plot with ggplot (new version)
NMDS.rarefied.df %>%
  mutate(age = ordered(as.factor(ifelse(NMDS.rarefied.df$age=="old", "Mature", ifelse(NMDS.rarefied.df$age=="medium", "Developing", "Emerging"))), levels = c("Emerging", "Developing", "Mature"))) %>%
  mutate(layer = ordered(as.factor(ifelse(NMDS.rarefied.df$layer=="epilimnion", "Epilimnion", ifelse(NMDS.rarefied.df$layer=="meta", "Metalimnion", "Hypolimnion"))), levels = c("Hypolimnion", "Metalimnion", "Epilimnion"))) %>%
  ggplot(aes(x=MDS1, y=MDS2, col=age, shape=layer)) +
 geom_point(size = 3) +
 scale_color_hue(l = 65, c = 100) +
 #geom_text(aes(label=sample_name),hjust=-0.15,vjust=0, size=3) +
 #stat_ellipse() +
 theme_bw() +
 labs(title = "NMDS of read mapping on PFAM entries\nrelated to carbon degradation (normalized)")
## Warning: Using shapes for an ordinal variable is not advised
```

![](Rmd_figs/Rmd-unnamed-chunk-4-1.png)<!-- -->

### Differential abundance analysis between ages in epilimnion samples

Generate input data for DESeq2, then create a CountDataSet object.


```r
# count table
de.epi.age.matrix <- metadata %>%
  subset(layer == "epilimnion") %>%
  inner_join(all.final.table, by = c("sample_name" = "sample")) %>%
  subset(select = c(sample_name, PFAM_ID, read_count)) %>%
  spread(PFAM_ID, read_count)

colnames.de.epi.age.matrix <- de.epi.age.matrix$sample_name
de.epi.age.matrix <- subset(de.epi.age.matrix, select = -sample_name) %>% t()
colnames(de.epi.age.matrix) <- colnames.de.epi.age.matrix

# metadata
de.epi.age.design <- metadata %>%
  subset(layer == "epilimnion") %>%
  inner_join(all.final.table, by = c("sample_name" = "sample")) %>%
  subset(select = c(sample_name, age)) %>%
  distinct()

de.epi.age.design <- de.epi.age.design[base::order(de.epi.age.design$sample_name),]

rownames.de.epi.age.design <- de.epi.age.design$sample_name
de.epi.age.design <- subset(de.epi.age.design, select = -sample_name)
rownames(de.epi.age.design) <- rownames.de.epi.age.design
## Warning: Setting row names on a tibble is deprecated.

# data need to be factor
de.epi.age.design$age <- factor(de.epi.age.design$age)

# Then create a CountDataSet object
ddsFullCountTable <- DESeqDataSetFromMatrix(
countData = de.epi.age.matrix,
colData = de.epi.age.design,
design = ~ age)
```

Run the DESeq2 pipeline and extract results for pairwise comparisons between ages


```r
dds <- DESeq(ddsFullCountTable)
## estimating size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing
res <- results(dds)
res.old_emerge    <- results(dds, contrast = c("age", "old", "emerge"))
res.old_medium    <- results(dds, contrast = c("age", "old", "medium"))
res.emerge_medium <- results(dds, contrast = c("age", "emerge", "medium"))
```

#### Old vs emerging


```r
plotMA.significant(res.old_emerge, cex = 0.6, main="DE analysis - MA plot: epilimnion, old vs emerging")
```

![](Rmd_figs/Rmd-unnamed-chunk-7-1.png)<!-- -->

```r
#DESeq2datatable(res.old_emerge)
DESeq2csv(res.old_emerge, file="results/DESeq2.epi.old_vs_emerging.csv")
```

#### Old vs medium


```r
plotMA.significant(res.old_medium, cex = 0.6, main="DE analysis - MA plot: epilimnion, old vs medium")
```

![](Rmd_figs/Rmd-unnamed-chunk-8-1.png)<!-- -->

```r
#DESeq2datatable(res.old_medium)
DESeq2csv(res.old_medium, file="results/DESeq2.epi.old_vs_medium.csv")
```

#### Emerging vs medium


```r
plotMA.significant(res.emerge_medium, cex = 0.6, main="DE analysis - MA plot: epilimnion, emerging vs medium")
```

![](Rmd_figs/Rmd-unnamed-chunk-9-1.png)<!-- -->

```r
#DESeq2datatable(res.emerge_medium)
DESeq2csv(res.emerge_medium, file="results/DESeq2.epi.emerging_vs_medium.csv")
```

### Differential abundance analysis between layers in old ponds

Generate input data for DESeq2.


```r
#### Generate input data for DESeq2:
### count table
de.old.matrix <- metadata %>%
  subset(age == "old") %>%
  inner_join(all.final.table, by = c("sample_name" = "sample")) %>%
  subset(select = c(sample_name, PFAM_ID, read_count)) %>%
  spread(PFAM_ID, read_count)

colnames.de.old.matrix <- de.old.matrix$sample_name
de.old.matrix <- subset(de.old.matrix, select = -sample_name) %>% t()
colnames(de.old.matrix) <- colnames.de.old.matrix

### metadata
de.old.design <- metadata %>%
  subset(age == "old") %>%
  inner_join(all.final.table, by = c("sample_name" = "sample")) %>%
  subset(select = c(sample_name, layer)) %>%
  distinct()

de.old.design <- de.old.design[base::order(de.old.design$sample_name),]

rownames.de.old.design <- de.old.design$sample_name
de.old.design <- subset(de.old.design, select = -sample_name)
rownames(de.old.design) <- rownames.de.old.design
## Warning: Setting row names on a tibble is deprecated.

# data need to be factor
de.old.design$layer <- factor(de.old.design$layer)

### Then create a CountDataSet object
de.old.ddsFullCountTable <- DESeqDataSetFromMatrix(
countData = de.old.matrix,
colData = de.old.design,
design = ~ layer)
```

Run the DESeq2 pipeline and extract results for pairwise comparisons between ages.


```r
de.old.dds <- DESeq(de.old.ddsFullCountTable)
## estimating size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing

de.old.res.epi_meta    <- results(de.old.dds, contrast = c("layer", "epilimnion", "meta"))
de.old.res.epi_hypo    <- results(de.old.dds, contrast = c("layer", "epilimnion", "hypo"))
de.old.res.hypo_meta   <- results(de.old.dds, contrast = c("layer", "hypo", "meta"))
```

#### Epilimnion vs metalimnion


```r
plotMA.significant(de.old.res.epi_meta, cex = 0.6, main="DA analysis - MA plot: old ponds, epilimnion vs metalimnion")
```

![](Rmd_figs/Rmd-unnamed-chunk-12-1.png)<!-- -->

```r
#DESeq2datatable(de.old.res.epi_meta)
DESeq2csv(de.old.res.epi_meta, file="results/DESeq2.old.epi_vs_meta.csv")
```

#### Epilimnion vs hypolimnion


```r
plotMA.significant(de.old.res.epi_hypo, cex = 0.6, main="DA analysis - MA plot: old ponds, epilimnion vs hypolimnion")
```

![](Rmd_figs/Rmd-unnamed-chunk-13-1.png)<!-- -->

```r
#DESeq2datatable(de.old.res.epi_hypo)
DESeq2csv(de.old.res.epi_hypo, file="results/DESeq2.old.epi_vs_hypo.csv")
```

#### Hypolimnion vs metalimnion


```r
plotMA.significant(de.old.res.hypo_meta, cex = 0.6, main="DA analysis - MA plot: old ponds, hypolimnion vs metalimnion")
```

![](Rmd_figs/Rmd-unnamed-chunk-14-1.png)<!-- -->

```r
#DESeq2datatable(de.old.res.hypo_meta)
DESeq2csv(de.old.res.hypo_meta, file="results/DESeq2.old.hypo_vs_meta.csv")
```

### Differential abundance analysis between oxygen content in old ponds

Generate input data for DESeq2.


```r
#### Generate input data for DESeq2:
### count table
# de.oxy.old.matrix <- metadata %>%
#   #mutate(oxygen_content = ifelse(metadata$layer == "hypo", "anoxic", "oxic")) %>%
#   subset(age == "old") %>%
#   inner_join(all.final.table, by = c("sample_name" = "sample")) %>%
#   subset(select = c(sample_name, PFAM_ID, read_count)) %>%
#   spread(PFAM_ID, read_count)
# 
# colnames.de.oxy.old.matrix <- de.oxy.old.matrix$sample_name
# de.oxy.old.matrix <- subset(de.oxy.old.matrix, select = -sample_name) %>% t()
# colnames(de.oxy.old.matrix) <- colnames.de.oxy.old.matrix

### metadata
de.oxy.old.design <- metadata %>%
  mutate(oxygen_content = ifelse(metadata$layer == "hypo", "anoxic", "oxic")) %>%
  subset(age == "old") %>%
  inner_join(all.final.table, by = c("sample_name" = "sample")) %>%
  subset(select = c(sample_name, oxygen_content)) %>%
  distinct()

de.oxy.old.design <- de.oxy.old.design[base::order(de.oxy.old.design$sample_name),]

rownames.de.oxy.old.design <- de.oxy.old.design$sample_name
de.oxy.old.design <- subset(de.oxy.old.design, select = -sample_name)
rownames(de.oxy.old.design) <- rownames.de.oxy.old.design
## Warning: Setting row names on a tibble is deprecated.

# data need to be factor
de.oxy.old.design$oxygen_content <- factor(de.oxy.old.design$oxygen_content)

### Then create a CountDataSet object
de.oxy.old.ddsFullCountTable <- DESeqDataSetFromMatrix(
countData = de.old.matrix,
colData = de.oxy.old.design,
design = ~ oxygen_content)
```

Run the DESeq2 pipeline and extract results for pairwise comparisons between ages.


```r
de.oxy.old.dds <- DESeq(de.oxy.old.ddsFullCountTable)
## estimating size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing

de.oxy.old.res.oxic_anoxic    <- results(de.oxy.old.dds, contrast = c("oxygen_content", "oxic", "anoxic"))
# de.old.res.epi_hypo    <- results(de.old.dds, contrast = c("layer", "epilimnion", "hypo"))
# de.old.res.hypo_meta   <- results(de.old.dds, contrast = c("layer", "hypo", "meta"))
```

#### Oxic vs anoxic


```r
plotMA.significant(de.oxy.old.res.oxic_anoxic, cex = 0.6, main="DA analysis - MA plot: old ponds, oxic vs anoxic")
```

![](Rmd_figs/Rmd-unnamed-chunk-17-1.png)<!-- -->

```r
#DESeq2datatable(de.oxy.old.res.oxic_anoxic)
DESeq2csv(de.oxy.old.res.oxic_anoxic, file="results/DESeq2.old.oxic_vs_anoxic.csv")
```

## Mapping of reads on contigs from coassembly of all samples

Map reads on contigs from coassembly of all samples. Then we check mappings on genes which have at least one of the HMMs of interest. The idea is to consider as match a whole gene even if only a part of it has a match with a HMM profile. If a gene has more than one match in PFAM, counts will be averaged for all HMMs.

### Data preparation


```r
# get PFAM entry lengths so we will calculate tpm when parsing the table
pfam.lengths <- read.table("data/Pfam-A.hmm.dat.lengths", header = TRUE, stringsAsFactors = FALSE) %>% 
                  mutate(PFAM_length_nt=PFAM_length*3)

all_pfam_files = list.files('counts', pattern="*all*", full.names = TRUE)
big.table.pfam <- data.frame(rbindlist(lapply(all_pfam_files, function(x) getBigTable.tpm.pfam(paste(x), pfam.lengths))))

# Read table with PFAM of interest
carbohydrate_degradation_pfams_tveit <- read_csv("data/carbohydrate_degradation_pfams_tveit.csv", 
                                                 col_names = FALSE, comment = "#", trim_ws = TRUE)
## Parsed with column specification:
## cols(
##   X1 = col_character(),
##   X2 = col_character(),
##   X3 = col_character(),
##   X4 = col_character(),
##   X5 = col_character()
## )

# Extract relevant data and save to file
big.table.pfam %>%
  select(-c(n_reads,PFAM_length_nt)) %>%
  filter(PFAM_ID %in% carbohydrate_degradation_pfams_tveit$X3) %>%
  spread(PFAM_ID, tpm) %>%
  t() %>%
  as.data.frame %>%
  write.csv(file="results/read_contigs_mapping_carbon_HMM.csv", quote=FALSE, col.names = FALSE)
```

### Differential abundance analysis


```r
all_pfam_files.unique = list.files('counts', pattern="*unique*", full.names = TRUE)
big.matrix.pfam.unique <- data.frame(rbindlist(lapply(all_pfam_files.unique, function(x) getBigTable.general(paste(x), ColNames = c("PFAM_ID", "read_count"))))) %>%
                            filter(!base::grepl("^_", PFAM_ID)) %>%
                            mutate_at("PFAM_ID", str_trunc, width=7, ellipsis="") %>%
                            filter(PFAM_ID %in% carbohydrate_degradation_pfams_tveit$X3) %>%
                            spread(PFAM_ID, read_count)

rownames.big.matrix.pfam.unique <- big.matrix.pfam.unique$sample
big.matrix.pfam.unique <- data.matrix(subset(big.matrix.pfam.unique, select = -sample)) 
rownames(big.matrix.pfam.unique) <- rownames.big.matrix.pfam.unique
```

### NMDS on normalized data

Perform NMDS on normalized data. Data are normalized by rarefying to the sample with the lowest number of reads mapped on PFAM HMMs of interest and plotted with ggplot.


```r
# source: http://geoffreyzahn.com/nmds-example/
#min_depth = min(rowSums(big.matrix.pfam.unique))

# normalize data
big.matrix.pfam.unique.rarefied <- as.matrix(round(rrarefy(big.matrix.pfam.unique, min(rowSums(big.matrix.pfam.unique)))))

# calculate distance matrix
big.matrix.pfam.unique.rarefied.sample_dist.rarefied <- as.matrix((vegdist(big.matrix.pfam.unique.rarefied, "bray")))

# perform actual NMDS
big.matrix.pfam.unique.rarefied.NMDS.rarefied <- metaMDS(big.matrix.pfam.unique.rarefied.sample_dist.rarefied)
## Run 0 stress 0.1046963 
## Run 1 stress 0.1422265 
## Run 2 stress 0.1046963 
## ... New best solution
## ... Procrustes: rmse 8.767143e-06  max resid 3.4709e-05 
## ... Similar to previous best
## Run 3 stress 0.1849999 
## Run 4 stress 0.1374446 
## Run 5 stress 0.170844 
## Run 6 stress 0.1046963 
## ... Procrustes: rmse 5.873271e-06  max resid 1.331745e-05 
## ... Similar to previous best
## Run 7 stress 0.1046963 
## ... Procrustes: rmse 3.050662e-06  max resid 1.122303e-05 
## ... Similar to previous best
## Run 8 stress 0.1046709 
## ... New best solution
## ... Procrustes: rmse 0.004692634  max resid 0.01849141 
## Run 9 stress 0.1046708 
## ... New best solution
## ... Procrustes: rmse 8.267654e-06  max resid 3.359155e-05 
## ... Similar to previous best
## Run 10 stress 0.1046709 
## ... Procrustes: rmse 3.62784e-05  max resid 8.472414e-05 
## ... Similar to previous best
## Run 11 stress 0.1046963 
## ... Procrustes: rmse 0.004694195  max resid 0.01849693 
## Run 12 stress 0.1374446 
## Run 13 stress 0.1678184 
## Run 14 stress 0.1587945 
## Run 15 stress 0.1046709 
## ... Procrustes: rmse 7.291912e-05  max resid 0.0003215842 
## ... Similar to previous best
## Run 16 stress 0.1374446 
## Run 17 stress 0.1374446 
## Run 18 stress 0.1046964 
## ... Procrustes: rmse 0.004696003  max resid 0.01849489 
## Run 19 stress 0.1421272 
## Run 20 stress 0.1046963 
## ... Procrustes: rmse 0.004692808  max resid 0.01849935 
## *** Solution reached

# build a data frame with NMDS coordinates and metadata
big.matrix.pfam.unique.rarefied.NMDS.rarefied.df = data.frame(MDS1 = big.matrix.pfam.unique.rarefied.NMDS.rarefied$points[,1],
                                                              MDS2 = big.matrix.pfam.unique.rarefied.NMDS.rarefied$points[,2])
big.matrix.pfam.unique.rarefied.NMDS.rarefied.df$sample_name <- rownames(big.matrix.pfam.unique.rarefied.NMDS.rarefied.df)
big.matrix.pfam.unique.rarefied.NMDS.rarefied.df <- big.matrix.pfam.unique.rarefied.NMDS.rarefied.df %>%
                                                      left_join(metadata)
## Joining, by = "sample_name"

# Plot with ggplot
# ggplot(big.matrix.pfam.unique.rarefied.NMDS.rarefied.df, aes(x=MDS1, y=MDS2, col=age, shape=layer)) +
#  geom_point() +
#  geom_text(aes(label=sample_name),hjust=-0.15,vjust=0, size=3) +
#  #stat_ellipse() +
#  theme_bw() +
#  labs(title = "NMDS of read mapping on PFAM entries annotated on contigs\nrelated to carbon degradation (normalized)")

big.matrix.pfam.unique.rarefied.NMDS.rarefied.df %>%
    mutate(age = ordered(as.factor(ifelse(NMDS.rarefied.df$age=="old", "Mature", ifelse(NMDS.rarefied.df$age=="medium", "Developing", "Emerging"))), levels = c("Emerging", "Developing", "Mature"))) %>%
    mutate(layer = ordered(as.factor(ifelse(NMDS.rarefied.df$layer=="epilimnion", "Epilimnion", ifelse(NMDS.rarefied.df$layer=="meta", "Metalimnion", "Hypolimnion"))), levels = c("Hypolimnion", "Metalimnion", "Epilimnion"))) %>%
    ggplot(aes(x=MDS1, y=MDS2, col=age, shape=layer)) +
 geom_point(size = 3) +
 scale_color_hue(l = 65, c = 100) + 
 #geom_text(aes(label=sample_name),hjust=-0.15,vjust=0, size=3) +
 #stat_ellipse() +
 theme_bw() +
 labs(title = "NMDS of read mapping on PFAM entries annotated on contigs\nrelated to carbon degradation (normalized)")
## Warning: Using shapes for an ordinal variable is not advised
```

![](Rmd_figs/Rmd-unnamed-chunk-20-1.png)<!-- -->

### Differential abundance analysis between ages in epilimnion samples

Generate input data for DESeq2, then create a CountDataSet object.


```r
# count table
contigs.de.epi.age.matrix <- data.frame(rbindlist(lapply(all_pfam_files.unique, function(x) getBigTable.general(paste(x), ColNames = c("PFAM_ID", "read_count"))))) %>%
                                    filter(!base::grepl("^_", PFAM_ID)) %>%
                                    mutate_at("PFAM_ID", str_trunc, width=7, ellipsis="") %>%
                                    #filter(PFAM_ID %in% carbohydrate_degradation_pfams_tveit$X3) %>%
                                    inner_join(subset(metadata, layer == "epilimnion"), by = c("sample" = "sample_name")) %>%
                                    subset(select = c(sample, PFAM_ID, read_count)) %>%
                                    spread(PFAM_ID, read_count)

# de.epi.age.matrix <- metadata %>%
#   subset(layer == "epilimnion") %>%
#   inner_join(all.final.table, by = c("sample_name" = "sample")) %>%
#   subset(select = c(sample_name, PFAM_ID, read_count)) %>%
#   spread(PFAM_ID, read_count)

colnames.contigs.de.epi.age.matrix <- contigs.de.epi.age.matrix$sample
contigs.de.epi.age.matrix <- subset(contigs.de.epi.age.matrix, select = -sample) %>% t()
colnames(contigs.de.epi.age.matrix) <- colnames.contigs.de.epi.age.matrix

# metadata
contigs.de.epi.age.design <- metadata %>%
  subset(layer == "epilimnion") %>%
  inner_join(data.frame(rbindlist(lapply(all_pfam_files.unique, function(x) getBigTable.general(paste(x), ColNames = c("PFAM_ID", "read_count"))))),
             by = c("sample_name" = "sample")) %>%
  subset(select = c(sample_name, age)) %>%
  distinct()

contigs.de.epi.age.design <- contigs.de.epi.age.design[base::order(contigs.de.epi.age.design$sample_name),]

rownames.contigs.de.epi.age.design <- contigs.de.epi.age.design$sample_name
contigs.de.epi.age.design <- subset(contigs.de.epi.age.design, select = -sample_name)
rownames(contigs.de.epi.age.design) <- rownames.contigs.de.epi.age.design
## Warning: Setting row names on a tibble is deprecated.

# data need to be factor
contigs.de.epi.age.design$age <- factor(contigs.de.epi.age.design$age)

# Then create a CountDataSet object
contigs.ddsFullCountTable <- DESeqDataSetFromMatrix(
countData = contigs.de.epi.age.matrix,
colData = contigs.de.epi.age.design,
design = ~ age)
```

Run the DESeq2 pipeline and extract results for pairwise comparisons between ages


```r
contigs.dds <- DESeq(contigs.ddsFullCountTable)
## estimating size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing
contigs.res <- results(contigs.dds)
contigs.res.old_emerge    <- results(contigs.dds, contrast = c("age", "old", "emerge"))
contigs.res.old_medium    <- results(contigs.dds, contrast = c("age", "old", "medium"))
contigs.res.emerge_medium <- results(contigs.dds, contrast = c("age", "emerge", "medium"))
```

#### Old vs emerging


```r
# plot DE analysis results
plotMA.significant(contigs.res.old_emerge, cex = 0.6, main="DE analysis - MA plot: epilimnion, old vs emerging")
```

![](Rmd_figs/Rmd-unnamed-chunk-23-1.png)<!-- -->

```r
# show results in dynamic table
#contigs.res.old_emerge[rownames(contigs.res.old_emerge) %in% carbohydrate_degradation_pfams_tveit$X3,] %>% DESeq2datatable()
# save results in csv file
dir.create("results", showWarnings = FALSE)
contigs.res.old_emerge[rownames(contigs.res.old_emerge) %in% carbohydrate_degradation_pfams_tveit$X3,] %>% write.csv(file="results/DESeq2.contigs.epi.old_vs_emerging.csv")
#DESeq2csv(contigs.res.old_emerge[rownames(contigs.res.old_emerge) %in% carbohydrate_degradation_pfams_tveit$X3,], file="DESeq2.contigs.epi.old_vs_emerging.csv")
```

#### Old vs medium


```r
# plot DE analysis results
plotMA.significant(contigs.res.old_medium, cex = 0.6, main="DE analysis - MA plot: epilimnion, old vs medium")
```

![](Rmd_figs/Rmd-unnamed-chunk-24-1.png)<!-- -->

```r
# show results in dynamic table
#contigs.res.old_medium[rownames(contigs.res.old_medium) %in% carbohydrate_degradation_pfams_tveit$X3,] %>% DESeq2datatable()
# save results in csv file
dir.create("results", showWarnings = FALSE)
contigs.res.old_medium[rownames(contigs.res.old_medium) %in% carbohydrate_degradation_pfams_tveit$X3,] %>% write.csv(file="results/DESeq2.contigs.epi.old_vs_medium.csv")
#DESeq2csv(contigs.res.old_medium[rownames(contigs.res.old_medium) %in% carbohydrate_degradation_pfams_tveit$X3,], file="DESeq2.contigs.epi.old_vs_emerging.csv")

#plotMA.significant(res.old_medium, cex = 0.6, main="DE analysis - MA plot: epilimnion, old vs medium")
#DESeq2datatable(res.old_medium)
#DESeq2csv(res.old_medium, file="DESeq2.epi.old_vs_medium.csv")
```

#### Emerging vs medium


```r
# plot DE analysis results
plotMA.significant(contigs.res.emerge_medium, cex = 0.6, main="DE analysis - MA plot: epilimnion, emerging vs medium")
```

![](Rmd_figs/Rmd-unnamed-chunk-25-1.png)<!-- -->

```r
# show results in dynamic table
#contigs.res.emerge_medium[rownames(contigs.res.emerge_medium) %in% carbohydrate_degradation_pfams_tveit$X3,] %>% DESeq2datatable()
# save results in csv file
contigs.res.emerge_medium[rownames(contigs.res.emerge_medium) %in% carbohydrate_degradation_pfams_tveit$X3,] %>% write.csv(file="results/DESeq2.contigs.epi.emerge_vs_medium.csv")

# plotMA.significant(res.emerge_medium, cex = 0.6, main="DE analysis - MA plot: epilimnion, emerging vs medium")
# DESeq2datatable(res.emerge_medium)
# DESeq2csv(res.emerge_medium, file="DESeq2.epi.emerging_vs_medium.csv")
```

### Differential abundance analysis between layers in old ponds

Generate input data for DESeq2.


```r
# #### Generate input data for DESeq2:
# ### count table
# de.old.matrix <- metadata %>%
#   subset(age == "old") %>%
#   inner_join(all.final.table, by = c("sample_name" = "sample")) %>%
#   subset(select = c(sample_name, PFAM_ID, read_count)) %>%
#   spread(PFAM_ID, read_count)
# 
# colnames.de.old.matrix <- de.old.matrix$sample_name
# de.old.matrix <- subset(de.old.matrix, select = -sample_name) %>% t()
# colnames(de.old.matrix) <- colnames.de.old.matrix
# 
# ### metadata
# de.old.design <- metadata %>%
#   subset(age == "old") %>%
#   inner_join(all.final.table, by = c("sample_name" = "sample")) %>%
#   subset(select = c(sample_name, layer)) %>%
#   distinct()
# 
# de.old.design <- de.old.design[base::order(de.old.design$sample_name),]
# 
# rownames.de.old.design <- de.old.design$sample_name
# de.old.design <- subset(de.old.design, select = -sample_name)
# rownames(de.old.design) <- rownames.de.old.design
# 
# # data need to be factor
# de.old.design$layer <- factor(de.old.design$layer)
# 
# ### Then create a CountDataSet object
# de.old.ddsFullCountTable <- DESeqDataSetFromMatrix(
# countData = de.old.matrix,
# colData = de.old.design,
# design = ~ layer)

############################ STARTS HERE
# count table
contigs.de.old.matrix <- data.frame(rbindlist(lapply(all_pfam_files.unique, function(x) getBigTable.general(paste(x), ColNames = c("PFAM_ID", "read_count"))))) %>%
                                    filter(!base::grepl("^_", PFAM_ID)) %>%
                                    mutate_at("PFAM_ID", str_trunc, width=7, ellipsis="") %>%
                                    inner_join(subset(metadata, age == "old"), by = c("sample" = "sample_name")) %>%
                                    subset(select = c(sample, PFAM_ID, read_count)) %>%
                                    spread(PFAM_ID, read_count)

# de.epi.age.matrix <- metadata %>%
#   subset(layer == "epilimnion") %>%
#   inner_join(all.final.table, by = c("sample_name" = "sample")) %>%
#   subset(select = c(sample_name, PFAM_ID, read_count)) %>%
#   spread(PFAM_ID, read_count)

colnames.contigs.de.old.matrix <- contigs.de.old.matrix$sample
contigs.de.old.matrix <- subset(contigs.de.old.matrix, select = -sample) %>% t()
colnames(contigs.de.old.matrix) <- colnames.contigs.de.old.matrix

# metadata
contigs.de.old.design <- metadata %>%
  subset(age == "old") %>%
  inner_join(data.frame(rbindlist(lapply(all_pfam_files.unique, function(x) getBigTable.general(paste(x), ColNames = c("PFAM_ID", "read_count"))))),
             by = c("sample_name" = "sample")) %>%
  subset(select = c(sample_name, layer)) %>%
  distinct()

contigs.de.old.design <- contigs.de.old.design[base::order(contigs.de.old.design$sample_name),]

rownames.contigs.de.old.design <- contigs.de.old.design$sample_name
contigs.de.old.design <- subset(contigs.de.old.design, select = -sample_name)
rownames(contigs.de.old.design) <- rownames.contigs.de.old.design
## Warning: Setting row names on a tibble is deprecated.

# data need to be factor
contigs.de.old.design$layer <- factor(contigs.de.old.design$layer)

# Then create a CountDataSet object
contigs.de.old.ddsFullCountTable <- DESeqDataSetFromMatrix(
countData = contigs.de.old.matrix,
colData = contigs.de.old.design,
design = ~ layer)
```

Run the DESeq2 pipeline and extract results for pairwise comparisons between ages.


```r
contigs.de.old.dds <- DESeq(contigs.de.old.ddsFullCountTable)
## estimating size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing
## -- replacing outliers and refitting for 89 genes
## -- DESeq argument 'minReplicatesForReplace' = 7 
## -- original counts are preserved in counts(dds)
## estimating dispersions
## fitting model and testing

contigs.de.old.res.epi_meta    <- results(contigs.de.old.dds, contrast = c("layer", "epilimnion", "meta"))
contigs.de.old.res.epi_hypo    <- results(contigs.de.old.dds, contrast = c("layer", "epilimnion", "hypo"))
contigs.de.old.res.hypo_meta   <- results(contigs.de.old.dds, contrast = c("layer", "hypo", "meta"))
```

#### Epilimnion vs metalimnion


```r
# plot DE analysis results
plotMA.significant(contigs.de.old.res.epi_meta, cex = 0.6, main="DA analysis - MA plot: old ponds, epilimnion vs metalimnion")
```

![](Rmd_figs/Rmd-unnamed-chunk-28-1.png)<!-- -->

```r
# show results in dynamic table
#contigs.de.old.res.epi_meta[rownames(contigs.de.old.res.epi_meta) %in% carbohydrate_degradation_pfams_tveit$X3,] %>% DESeq2datatable()
# save results in csv file
#DESeq2csv(contigs.de.old.res.epi_meta, file="DESeq2.old.epi_vs_meta.csv")
contigs.de.old.res.epi_meta[rownames(contigs.de.old.res.epi_meta) %in% carbohydrate_degradation_pfams_tveit$X3,] %>% write.csv(file="results/DESeq2.contigs.old.epi_vs_meta.csv")
```

#### Epilimnion vs hypolimnion


```r
# plot DE analysis results
plotMA.significant(contigs.de.old.res.epi_hypo, cex = 0.6, main="DA analysis - MA plot: old ponds, epilimnion vs hypolimnion")
```

![](Rmd_figs/Rmd-unnamed-chunk-29-1.png)<!-- -->

```r
# show results in dynamic table
#contigs.de.old.res.epi_hypo[rownames(contigs.de.old.res.epi_hypo) %in% carbohydrate_degradation_pfams_tveit$X3,] %>% DESeq2datatable()
# save results in csv file
#DESeq2csv(contigs.de.old.res.epi_meta, file="DESeq2.old.epi_vs_meta.csv")
contigs.de.old.res.epi_hypo[rownames(contigs.de.old.res.epi_hypo) %in% carbohydrate_degradation_pfams_tveit$X3,] %>% write.csv(file="results/DESeq2.contigs.old.epi_vs_hypo.csv")

# plotMA.significant(contigs.de.old.res.epi_hypo, cex = 0.6, main="DA analysis - MA plot: old ponds, epilimnion vs hypolimnion")
# DESeq2datatable(contigs.de.old.res.epi_hypo)
# DESeq2csv(contigs.de.old.res.epi_hypo, file="DESeq2.old.epi_vs_hypo.csv")
```

#### Hypolimnion vs metalimnion


```r
# plot DE analysis results
plotMA.significant(contigs.de.old.res.hypo_meta, cex = 0.6, main="DA analysis - MA plot: old ponds, hypolimnion vs metalimnion")
```

![](Rmd_figs/Rmd-unnamed-chunk-30-1.png)<!-- -->

```r
# show results in dynamic table
#contigs.de.old.res.hypo_meta[rownames(contigs.de.old.res.hypo_meta) %in% carbohydrate_degradation_pfams_tveit$X3,] %>% DESeq2datatable()
# save results in csv file
#DESeq2csv(contigs.de.old.res.epi_meta, file="DESeq2.old.epi_vs_meta.csv")
contigs.de.old.res.hypo_meta[rownames(contigs.de.old.res.hypo_meta) %in% carbohydrate_degradation_pfams_tveit$X3,] %>% write.csv(file="results/DESeq2.contigs.old.hypo_vs_meta.csv")

# plotMA.significant(de.old.res.hypo_meta, cex = 0.6, main="DA analysis - MA plot: old ponds, hypolimnion vs metalimnion")
# DESeq2datatable(de.old.res.hypo_meta)
# DESeq2csv(de.old.res.hypo_meta, file="DESeq2.old.hypo_vs_meta.csv")
```

### Differential abundance analysis between oxic and anoxic layers in old ponds

**Epilimnion** and **metalimnion** can be considered as **oxic** layers, while **hypolimnion** can be considered **anoxic**, according to environmental parameters.

Generate input data for DESeq2.


```r
############################ STARTS HERE
### The matrix is already built in the previous analysis
# count table
# contigs.de.old.matrix <- data.frame(rbindlist(lapply(all_pfam_files.unique, function(x) getBigTable.general(paste(x), ColNames = c("PFAM_ID", "read_count"))))) %>%
#                                     filter(!base::grepl("^_", PFAM_ID)) %>%
#                                     mutate_at("PFAM_ID", str_trunc, width=7, ellipsis="") %>%
#                                     inner_join(subset(metadata, age == "old"), by = c("sample" = "sample_name")) %>%
#                                     subset(select = c(sample, PFAM_ID, read_count)) %>%
#                                     spread(PFAM_ID, read_count)
# 
# colnames.contigs.de.old.matrix <- contigs.de.old.matrix$sample
# contigs.de.old.matrix <- subset(contigs.de.old.matrix, select = -sample) %>% t()
# colnames(contigs.de.old.matrix) <- colnames.contigs.de.old.matrix

### METADATA
# Generate additional column with "oxic" or "anoxic" status according to layer 
# metadata
contigs.oxy.de.old.design <- metadata %>%
  mutate(oxygen_content = ifelse(metadata$layer == "hypo", "anoxic", "oxic")) %>%
  subset(age == "old") %>%
  inner_join(data.frame(rbindlist(lapply(all_pfam_files.unique, function(x) getBigTable.general(paste(x), ColNames = c("PFAM_ID", "read_count"))))),
             by = c("sample_name" = "sample")) %>%
  subset(select = c(sample_name, oxygen_content)) %>%
  distinct()

contigs.oxy.de.old.design <- contigs.oxy.de.old.design[base::order(contigs.oxy.de.old.design$sample_name),]

rownames.contigs.oxy.de.old.design <- contigs.oxy.de.old.design$sample_name
contigs.oxy.de.old.design <- subset(contigs.oxy.de.old.design, select = -sample_name)
rownames(contigs.oxy.de.old.design) <- rownames.contigs.oxy.de.old.design
## Warning: Setting row names on a tibble is deprecated.

# data need to be factor
contigs.oxy.de.old.design$oxygen_content <- factor(contigs.oxy.de.old.design$oxygen_content)

# Then create a CountDataSet object
contigs.oxy.de.old.ddsFullCountTable <- DESeqDataSetFromMatrix(
countData = contigs.de.old.matrix,
colData = contigs.oxy.de.old.design,
design = ~ oxygen_content)
```

Run the DESeq2 pipeline and extract results for pairwise comparisons between oxygen content.


```r
contigs.oxy.de.old.dds <- DESeq(contigs.oxy.de.old.ddsFullCountTable)
## estimating size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing
## -- replacing outliers and refitting for 93 genes
## -- DESeq argument 'minReplicatesForReplace' = 7 
## -- original counts are preserved in counts(dds)
## estimating dispersions
## fitting model and testing

contigs.oxy.de.old.res <- results(contigs.oxy.de.old.dds, contrast = c("oxygen_content", "oxic", "anoxic"))
```

#### Oxic vs anoxic


```r
# plot DE analysis results
plotMA.significant(contigs.oxy.de.old.res, cex = 0.6, main="DA analysis - MA plot: old ponds, oxic vs anoxic")
```

![](Rmd_figs/Rmd-unnamed-chunk-33-1.png)<!-- -->

```r
# show results in dynamic table
#contigs.oxy.de.old.res[rownames(contigs.oxy.de.old.res) %in% carbohydrate_degradation_pfams_tveit$X3,] %>% DESeq2datatable()
# save results in csv file
#DESeq2csv(contigs.de.old.res.epi_meta, file="DESeq2.old.epi_vs_meta.csv")
contigs.oxy.de.old.res[rownames(contigs.oxy.de.old.res) %in% carbohydrate_degradation_pfams_tveit$X3,] %>% write.csv(file="results/DESeq2.contigs.oxy.de.old.oxic_vs_anoxic.csv")
```

## Statistical tests on abundance of carbon degradation related PFAMs 

**Testing if there is a difference in total abundance abundance of degradation related PFAMs in the oxic layer across the ponds representing different ages**

Perform a one-way ANOVA.


```r
aov_carbon_PFAM_oxic_old <- data.frame(sample_name=row.names(all.final.matrix), total_mapped_reads=rowSums(all.final.matrix)) %>%
  left_join(metadata.dataset_size) %>%
  mutate("total_mapped_reads_percent"=total_mapped_reads*100/total_reads) %>%
  left_join(metadata, by = ("sample_name")) %>%
  subset(age == "old") %>%
  subset(layer != "hypo") %>%
  dplyr::select(age, total_mapped_reads_percent)

aov_carbon_PFAM_oxic_medium <- data.frame(sample_name=row.names(all.final.matrix), total_mapped_reads=rowSums(all.final.matrix)) %>%
  left_join(metadata.dataset_size) %>%
  mutate("total_mapped_reads_percent"=total_mapped_reads*100/total_reads) %>%
  left_join(metadata, by = ("sample_name")) %>%
  subset(age == "medium") %>%
  subset(layer != "hypo") %>%
  dplyr::select(age, total_mapped_reads_percent)

aov_carbon_PFAM_oxic_emerge <- data.frame(sample_name=row.names(all.final.matrix), total_mapped_reads=rowSums(all.final.matrix)) %>%
  left_join(metadata.dataset_size) %>%
  mutate("total_mapped_reads_percent"=total_mapped_reads*100/total_reads) %>%
  left_join(metadata, by = ("sample_name")) %>%
  subset(age == "emerge") %>%
  subset(layer != "hypo") %>%
  dplyr::select(age, total_mapped_reads_percent)

aov.dataframe <- rbind(aov_carbon_PFAM_oxic_old, aov_carbon_PFAM_oxic_medium, aov_carbon_PFAM_oxic_emerge)
aov.dataframe$age <- base::as.factor(aov.dataframe$age)

#ggboxplot(aov.dataframe, x = "age", y = "total_mapped_reads_percent", color = "age")

# ANOVA test with assumption of equal variances
summary(aov(total_mapped_reads_percent ~ age, data = aov.dataframe))
```

```
##             Df   Sum Sq   Mean Sq F value Pr(>F)  
## age          2 0.001207 0.0006037   3.289 0.0654 .
## Residuals   15 0.002753 0.0001835                 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

One-way ANOVA didn't point out any significant difference, although the p-value is slighlty above the significance value (p-value = 0.0654).

**Testing if there is a difference in total abundance abundance of degradation related PFAMs in the oxic vs anoxic layer in the old ponds**


```r
carbon_PFAM_oxic_old <- data.frame(sample_name=row.names(all.final.matrix), total_mapped_reads=rowSums(all.final.matrix)) %>%
  left_join(metadata.dataset_size) %>%
  mutate("total_mapped_reads_percent"=total_mapped_reads*100/total_reads) %>%
  left_join(metadata, by = ("sample_name")) %>%
  subset(age == "old") %>%
  subset(layer != "hypo") %>%
  dplyr::select(total_mapped_reads_percent)

carbon_PFAM_anoxic_old <- data.frame(sample_name=row.names(all.final.matrix), total_mapped_reads=rowSums(all.final.matrix)) %>%
  left_join(metadata.dataset_size) %>%
  mutate("total_mapped_reads_percent"=total_mapped_reads*100/total_reads) %>%
  left_join(metadata, by = ("sample_name")) %>%
  subset(age == "old") %>%
  subset(layer == "hypo") %>%
  dplyr::select(total_mapped_reads_percent)

data.frame(oxygen_content="oxic", percent_read_mapped=carbon_PFAM_oxic_old$total_mapped_reads_percent) %>%
  rbind(data.frame(oxygen_content="anoxic", percent_read_mapped=carbon_PFAM_anoxic_old$total_mapped_reads_percent)) %>%
  ggboxplot(x="oxygen_content", y="percent_read_mapped", color="oxygen_content")
```

![](Rmd_figs/Rmd-unnamed-chunk-35-1.png)<!-- -->

```r
t_test_old_oxic_vs_anoxic <- t.test(carbon_PFAM_oxic_old, carbon_PFAM_anoxic_old)

t_test_old_oxic_vs_anoxic
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  carbon_PFAM_oxic_old and carbon_PFAM_anoxic_old
## t = -12.824, df = 12.376, p-value = 1.634e-08
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.03344556 -0.02375891
## sample estimates:
##  mean of x  mean of y 
## 0.06608594 0.09468817
```

Welch Two Sample t-test suggested a significant difference (p-value = 1.634-08).

