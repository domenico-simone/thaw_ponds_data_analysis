library(data.table)
library(tidyverse)
library(optparse)

option_list = list(
    make_option(c("-a", "--annotation_file"), type="character", default=NULL,
        help="cds annotation file [default= %default] (MANDATORY)", metavar="character"),
    make_option(c("-b", "--hmmscan_res"), type="character", default=NULL,
        help="hmmscan result table, domtblout format [default= %default] (MANDATORY)"),
    make_option(c("-o", "--output_file"), type="character", default=NULL,
        help="gff file with annotation of hmmscan hits  [default= %default] (MANDATORY)"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

annotation.file <- opt$annotation_file
hmmscan.res <- opt$hmmscan_res
outfile <- opt$output_file

### for interactive testing
# setwd("/home/adm2/Dropbox/Research/slu/thaw_ponds/thaw_ponds/thawpondsdata")
# annotation.file <- "prodigal/thawponds_assembly.cds.out"
# hmmscan.res <- "hmmer/split500000/thawponds_assembly.cds.split500000.all.hmmer_pfam.domtblout"
# outfile <- "prodigal/thawponds_assembly.cds.hmmer_pfam.gff.toy"

# Parse CDS annotation

#a <- fread("grep -v '^#' thawponds_assembly.cds.out",
#a <- fread(cmd = paste("grep -v '^#'", annotation.file),
a <- fread(paste("grep -v '^#'", annotation.file),
            sep = "\t",
            stringsAsFactors = FALSE,
            select = c(1,4,5,7),
            col.names = c("contigID",
                        "start",
                        "end",
                        "strand")) %>% group_by(contigID) %>%
                                    mutate(count = sequence(n())) %>%
                                    mutate(CDS_ID = paste(contigID, count, sep = "_")) %>%
                                    subset(select = -count) %>% ungroup()

# Read hmmscan result table

#b <- fread("grep -v '^#' thawponds_assembly.cds.split500000.all.hmmer_pfam.domtblout",
b <- fread(paste("grep -v '^#'", hmmscan.res),
            sep = " ",
            stringsAsFactors = FALSE,
            select = c(1,5,18,19),
            col.names = c("CDS_ID",
                        "PFAM_ID",
                        "hit_start",
                        "hit_end"))

# Conditionally join to hmmscan hits to obtain gff file

c <- b %>%
  left_join(a) %>%
  mutate(c_start=if_else(strand == "+", start+(hit_start-1)*3, end-hit_end*3)) %>%
  mutate(c_end=if_else(strand == "+", start+hit_end*3, end-(hit_start-1)*3)) %>%
  mutate(software="hmmer_3.2.1") %>% mutate(feature="Domain") %>% mutate(unk1=".") %>% mutate(last=paste0("PFAM_ID=", PFAM_ID, ";CDS_ID=", CDS_ID)) %>%
  base::subset(select = c(contigID, software, feature, c_start, c_end, unk1, strand, unk1, last))

# Write out gff file

write.table(c, file = outfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)