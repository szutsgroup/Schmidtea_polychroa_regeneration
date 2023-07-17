suppressMessages(library(dplyr))
suppressMessages(library(argparser))

p <- arg_parser("Filter getaf output by number of zeroes, coverage limits and supporting reads")
p <- add_argument(p, "--input", help = "Path to getaf report")
p <- add_argument(p, "--zeromin", help = "Minimum number of zeroes in line")
p <- add_argument(p, "--zeromax", help = "Maximum number of zeroes in line")
p <- add_argument(p, "--covmin", help = "Minimum value of mean coverage in line")
p <- add_argument(p, "--covmax", help = "Maximum value of mean coverage in line")
p <- add_argument(p, "--maxsuppmin", help = "Minimum value of max support in line")
p <- add_argument(p, "--minnonzerosupp", help = "Minimum value of the smallest non zero support in line")

argv <- parse_args(p)

g <- read.delim(argv$input, sep = " ", stringsAsFactors = FALSE)
# first four columns as chrom, pos, ref, alt, and each sample has two more cols
samplenum = (dim(g)[2] - 4) / 2
firstsuppcol = 5
firstcovcol = 4 + samplenum + 1

apply(g[,firstsuppcol:(firstsuppcol + samplenum - 1)], 1, function(x) sum(x == 0)) -> g$Zerosum
apply(g[,firstsuppcol:(firstsuppcol + samplenum - 1)], 1, max) -> g$Maxsupp
apply(g[,firstsuppcol:(firstsuppcol + samplenum - 1)], 1, function(x) ifelse(max(x) == 0, NA, min(x[which(x > 0)]))) -> g$MinnonzeroSupp
apply(g[,firstcovcol:(firstcovcol + samplenum - 1)], 1, mean) -> g$Covmean

filter(g, Maxsupp > as.numeric(argv$maxsuppmin),
          Covmean > as.numeric(argv$covmin),
          Covmean < as.numeric(argv$covmax),
          Zerosum > as.numeric(argv$zeromin),
          Zerosum < as.numeric(argv$zeromax),
          MinnonzeroSupp > as.numeric(argv$minnonzerosupp)) -> gg
if (nrow(gg) > 0) {
suppressMessages(library(Rsamtools))

checkCigar <- function(f, pos) {
    cigar = gsub("([0-9]+[A-Z])", " \\1 ", f$cigar) %>% strsplit(" ") %>% lapply(function(x) x[which(x != "")])
    rr = c()
    for (i in seq_along(cigar)) {
        ret = -2
        tmp = cigar[[i]]
        co = 0
        for (j in tmp) {
            a = as.numeric(substr(j, 1, nchar(j)-1))
            b = substr(j, nchar(j), nchar(j))
            if (b == "M") co = co+a
            if (b == "D") ret = co-1
            if (b == "I") ret = co
            }
        if (ret == (pos - f$pos[i])) {
            rr = c(rr, TRUE)
        } else {
            rr = c(rr, FALSE)
        }
    }
    return(rr)
}

getMQ = function(file, chrom, pos, ref, alt) {
    c("mapq", "pos", "seq", "cigar") -> what
    GRanges(chrom, IRanges(start = pos, end = pos)) -> whi
    sb = scanBam(file, param = ScanBamParam(which = whi, what = what))
    f = sb[[1]]
    pos = as.numeric(pos)
    ref = as.character(ref)
    alt = as.character(alt)

    if (nchar(ref) == nchar(alt)) {
        ok <- which(substr(f$seq, pos-f$pos+1, pos-f$pos+1) == alt)
    } else {
        ok <- which(checkCigar(f, pos))
    }

    return(list(mutated_mapq =  mean(f$mapq[ok]), other_mapq = mean(f$mapq[-ok]) ))
}

gg$Maxsupp_sample_mutMQ = gg$Maxsupp_sample_refMQ = "NA"
for (i in 1:nrow(gg)) {
    Sample = gsub("_Supp", "", names(gg)[which(gg[i,firstsuppcol:(firstsuppcol + samplenum - 1)] == gg[i,"Maxsupp"])[1]+4])
    tmp = getMQ(paste0("../../aligned/", Sample, "_RMdup.bam"), gg$Chrom[i], gg$Pos[i], gg$Ref[i], gg$Alt[i])
    gg$Maxsupp_sample_mutMQ[i] <- round(tmp$mutated_mapq, digits = 2)
    gg$Maxsupp_sample_refMQ[i] <- round(tmp$other_mapq, digits = 2)
}
gg <- filter(gg, as.numeric(Maxsupp_sample_refMQ) > 40)

write.table(gg, gsub(".out", "_filter.out", argv$input), append = TRUE, quote = FALSE, row.names = FALSE, col.names = !file.exists(gsub(".out", "_filter.out", argv$input)))
}
