suppressMessages(library(argparser))
suppressMessages(library(pafr))
suppressMessages(library(dplyr))
suppressMessages(library(Rsamtools))

p <- arg_parser("Purge assembly by checking minimap sel-mappings")
p <- add_argument(p, "--dir", help = "Path to the folder with the ID split assembly")
p <- add_argument(p, "--outdir", help = "Path to the output folder")
p <- add_argument(p, "--paf", help = "Path to minimap2 self alignment")
p <- add_argument(p, "--re", help = "Path to repeatmasker output")

argv <- parse_args(p)

paf <- read_paf(argv$paf)
re <- read.table(argv$re)

filter(paf, cm >= 18, s2 < 0.5*s1, tname != qname, nmatch / alen > 0.4) -> u
gr <- GRanges(seqnames = u$qname, IRanges(start = u$qstart, end = u$qend))
red <- GenomicRanges::reduce(gr, with.revmap = TRUE)
u$revmap_length <- NA
for (i in as.data.frame(red)$revmap) {
  for (j in i) {
    u$revmap_length[j] <- length(i)
  }
}

filter(u, revmap_length == 1, mapq == 60, alen > 500) -> u2

gru <- GRanges(seqnames = u2$qname, IRanges(start = u2$qstart, end = u2$qend))
grre <- GRanges(seqnames = re$V5, IRanges(start = re$V6, end = re$V7))
setd <- IRanges::intersect(gru, grre)
olaps <- findOverlaps(query = setd, subject = gru)
gru$rep_overlap <- 0
for (i in unique(subjectHits(olaps))) {
  tmp <- olaps[which(subjectHits(olaps) == i)]
  gru$rep_overlap[i] <- max(width(setd[queryHits(tmp)])) / width(gru)[i]
}
u2$longest_repeat_overlap <- gru$rep_overlap

filter(u2, longest_repeat_overlap < 0.4) -> u3

# 1. for each u3 line, select the shorter contig, 
# and extract the region to remove
u3$ID <- 1:nrow(u3)
to_remove <- vector(mode = "list", length = nrow(u3))
names(to_remove) <- u3$ID
for (i in u3$ID) {
  if (u3$qlen[i] >= u3$tlen[i]) {
    to_remove[[i]] <- as.data.frame(u3[i, c("tname", "tstart", "tend")])
      } else {
    to_remove[[i]] <- as.data.frame(u3[i, c("qname", "qstart", "qend")])
      }
  names(to_remove[[i]]) <- c("seqnames", "start", "end")
}
to_remove_df <- do.call(rbind.data.frame, to_remove)

# 2. consolidate the data.frame, by 
# a) removing duplicate hits (originating from dual alignmens),
# b) collapsing alignments on the same contig into 1 list element and 
# c) collapse overlapping alignment regions

# a)
to_remove_df$id1 = apply(to_remove_df, 1, function(x) paste0(x[1], "_", as.numeric(x[2]), collapse = ""))
to_remove_df$n1 = table(to_remove_df$id1)[to_remove_df$id1]
to_remove_df$n0 = table(to_remove_df$seqnames)[to_remove_df$seqnames]
to_remove[which(duplicated(to_remove_df$id1))] <- NA
to_remove_df$ID <- 1:nrow(to_remove_df)
# b)
filter(to_remove_df,  n0 > 1) %>% arrange(seqnames, start) %>% .$seqnames %>% unique -> multaln_contigs
for (i in multaln_contigs) {
  tmp <- filter(to_remove_df, seqnames == i) %>% arrange(start) %>% distinct(start, .keep_all = TRUE)
  ok <- tmp$ID[1]
  to_remove[[ok]] = tmp[, c("seqnames", "start", "end")]
  to_remove[tmp$ID[-1]] <- NA
  # c)
  GRanges(to_remove[[ok]]) %>%
    GenomicRanges::reduce() %>%
    as.data.frame() %>%
    .[c("seqnames", "start", "end")] -> to_remove[[ok]]
}

wh <- which(sapply(to_remove, function(x) !all(is.na(x))))
to_remove <- to_remove[wh]

# 3. for each interesting contig, read in the fasta,
# and create a multiple fasta file, leaving out the regions to be removed

for (g in to_remove) {
  n = g$seqnames[1]
  fas = FaFile(paste0(argv$dir, "scaffolds.id_", n, ".fasta"))
  seq = getSeq(fas, GRanges(g$seqnames[1], IRanges(start = 1, end = seqlengths(fas))))
  subseq(rep(seq, nrow(g)+1), start = c(1, g$end), end = c(g$start, width(seq))) -> out
  names(out) <- paste0(names(out), ":", c(1, g$end), "-", c(g$start, width(seq)))
  writeXStringSet(out, filepath = paste0(argv$outdir, n, "_deleted", nrow(g), ".fa"))
}
