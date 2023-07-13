suppressMessages(library(tidyverse))
suppressMessages(library(argparser))

p <- arg_parser("Use Trinotate report to create pcf and pcfl transcript sets")
p <- add_argument(p, "--report", help = "Path to Trinotate report")
p <- add_argument(p, "--fasta", help = "Path to raw trinity assembly")

argv <- parse_args(p)
print(argv)

nota <- read.delim(argv$report)

nota[,c("X.gene_id", "transcript_id", "sprot_Top_BLASTX_hit", "prot_coords", "Pfam")]  %>% 
  mutate(prot_len = gsub("\\[.\\]", "", prot_coords),
         prot_len = (as.numeric(sapply(strsplit(prot_len, split = "-"), "[", 2)) - as.numeric(sapply(strsplit(prot_len, split = "-"), "[", 1)) + 1)/3) -> g 
         
g %>% 
  unique() %>% 
  filter(Pfam != "." | sprot_Top_BLASTX_hit != "." | prot_len > 100) -> dd

library(seqinr)

fa <- read.fasta(file = argv$fasta)

for (i in seq_along(fa)) {
  if (attr(fa[[i]], "name") %in% dd$transcript_id) {
    write.fasta(sequences = fa[[i]], names = attr(fa[[i]], "Annot"), file.out = "trinity_out_dir.Trinity.pcf.fasta", 
                open = "a")
  }
}

arrange(dd, X.gene_id, desc(prot_len)) %>% group_by(X.gene_id) %>% filter(row_number() == 1) ->d3

for (i in seq_along(fa)) {
  if (attr(fa[[i]], "name") %in% d3$transcript_id) {
    write.fasta(sequences = fa[[i]], names = attr(fa[[i]], "Annot"), file.out = "trinity_out_dir.Trinity.pcfl.fasta", 
                open = "a")
  }
}



