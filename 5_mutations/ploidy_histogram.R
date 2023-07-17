suppressMessages(library(dplyr))
suppressMessages(library(argparser))

p <- arg_parser("Draw histograms of S. polychroa SNP allele frequencies")
p <- add_argument(p, "--input", help = "Input table")
argv = parse_args(p)
g = read.delim(argv$input, header = FALSE)

if (grepl("control", argv$input)) samples <- c("Pc1", "Pc2", "Fc1", "Fc2", "Pc3", "Pc4", "Fc4A", "Fc4B")
if (grepl("regenerated", argv$input)) samples <- c("Pr1", "Pr2", "Pr3", "Pr4", "R4", "R3", "R1", "R2", "Fr1", "Fr4A", "Fr4B", "Fr3")
for (i in 3:length(g)) {
  tmp <- g[,i]
  sapply(tmp, function(x) as.numeric(strsplit(x, ",")[[1]])) %>% 
    apply(2, function(x) x[2]/(x[1]+x[2])) -> g[,i]
}

