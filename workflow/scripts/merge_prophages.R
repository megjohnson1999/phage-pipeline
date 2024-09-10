library(readr)
library(tidyverse)
library(Biostrings)

# Function to check if 2 regions overlap
check_overlap <- function(region1, region2) {
  if (region1[2] >= region2[1] && region1[1] <= region2[2]) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

genomad_path <- file.path(snakemake@input[["genomad"]], 
                                          "final_filtered_contigs_find_proviruses", 
                                          "final_filtered_contigs_provirus.tsv")
genomad <- read_tsv(genomad_path) %>%
  separate_wider_delim(source_seq, '_', names=c('Node', 'contig', 'x', 'length_c', 'y', 'cov')) %>%
  as.data.frame() %>%
  select(contig, start, end) %>%
  mutate(tool = "genomad")


phispy_path <- file.path(snakemake@input[["phispy"]],
                         "prophage.tsv") 
phispy <- read_tsv(phispy_path) %>%
  separate_wider_delim('Prophage number', '_', names=c('pp', 'pp_number')) %>%
  separate_wider_delim(Contig, '_', names=c('Node', 'contig')) %>%
  as.data.frame() %>%
  rename(start = Start,
         end = Stop) %>%
  mutate(tool = "phispy")

# Set up empty data frame to hold unique-to-phispy prophages
phispy_unique <- data.frame(matrix(ncol=4, nrow=0))
cols <- c("contig", "start", "end", "tool")
colnames(phispy_unique) <- cols
shared_contigs <- intersect(genomad$contig, phispy$contig)

# Get all phispy prophages on unique contigs
for (row in 1:nrow(phispy)) {
  if (!(phispy$contig[row] %in% genomad$contig)) {
    phispy_unique <- rbind(phispy_unique, phispy[row, ])
  }
}

# Get phispy prophages on shared prophages but with no overlap with any genomad prophage
for (c in shared_contigs){
  genomad_shared <- genomad %>%
    filter(contig == c)
  phispy_shared <- phispy %>%
    filter(contig == c)
  for (i in 1:nrow(phispy_shared)){
    overlap_found <- FALSE
    for (j in 1:nrow(genomad_shared)){
      r1 <- c(phispy_shared$start[i], phispy_shared$end[i])
      r2 <- c(genomad_shared$start[j], genomad_shared$end[j])
      if (check_overlap(r1, r2) == TRUE){
        overlap_found <- TRUE
        break
      }
    }
    if (overlap_found == FALSE){
      phispy_unique <- bind_rows(phispy_unique, phispy_shared[i, ])
    }
  }
}


# Get pp_number for unique-to-phispy contigs
pp_num <- as.numeric(phispy_unique$pp_number)
fasta <- readDNAStringSet(file.path(snakemake@input[["phispy"]], 
                                          "phage.fasta"))
fasta_unique <- fasta[pp_num]
writeXStringSet(fasta_unique, filepath = snakemake@output[["fasta"]])


phispy_unique <- phispy_unique %>%
  select(contig, start, end, tool)

# Get taxonomy output in the right format
CAT_path <- file.path(snakemake@input[["CAT"]],
                         "contig.taxonomy")
CAT <- read_tsv(CAT_path) %>%
  separate_wider_delim('# contig', '_', names=c('Node', 'contig', 'x', 'length_c', 'y', 'cov')) %>%  
  as.data.frame() %>%
  select(contig, superkingdom:species)

# Merge the tables
final_prophage_table <- rbind(genomad, phispy_unique) %>%
  mutate(contig = as.numeric(contig)) %>%
  arrange(contig) %>%
  merge(CAT, by='contig')
  
write.table(final_prophage_table, snakemake@output[["table"]], row.names=FALSE, sep="\t", quote=FALSE)
