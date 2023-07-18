suppressMessages(library(dplyr))
suppressMessages(library(argparser))

p <- arg_parser("Draw histograms of S. polychroa SNP allele frequencies")
p <- add_argument(p, "--input", help = "Input table")
argv = parse_args(p)
g = read.delim(argv$input, header = FALSE)

if (grepl("control", argv$input)) samples <- c("Pc1", "Pc2", "Fc1", "Fc2", "Pc3", "Pc4", "Fc4A", "Fc4B")
if (grepl("regenerated", argv$input)) samples <- c("Pr1", "Pr2", "Pr3", "Pr4", "R4", "R3", "R1", "R2", "Fr1", "Fr4A", "Fr4B", "Fr3")
names(g) <- c("Scaffold", "Pos", "Ref", "Alt", samples)

for (i in 5:length(g)) {
  tmp <- g[,i]
  sapply(tmp, function(x) as.numeric(strsplit(x, ",")[[1]])) %>% 
    apply(2, function(x) x[2]/(x[1]+x[2])) -> g[,i]
}

g$Zerosum <- apply(g[,5:length(g)], 1, function(x) sum(x == 0))
g$afsd <- apply(g[,3:10], 1, function(x) mean(abs((x-mean(x)))))
g$ZscoreMax <- apply(g[,3:10], 1, function(x) max((x-mean(x))/sd(x)))

pdf(paste0("/Projects/Planaria/Mutations/", gsub("", "", argv$input), "_ploidy.pdf"), width = 20, height = 13)
par(mfrow = c(3, 4))
for (i in 5:(length(g)-3)) {
  
  filter(g, Zerosum == 0) %>% 
    with(hist(get(paste0("V", i)), breaks = seq(0, 1, .01), plot = FALSE)) -> hh
  print(hh$counts)
  plot(hh$counts, type = "h", x = seq(.01, 1, .01), 
       xlab = "Allele frequency", las = 1, 
       ylab = "Frequency", main = samples[(i-2)], xaxt = "n")
  axis(side = 1, at = seq(0, 1, .1))
  points(hh$counts, type = "p", pch = 16, x = seq(.01, 1, .01))
  lowess(hh$counts, f = .1) -> lo
  lines(lo$y, col = "red", lwd = 4, x = seq(.01, 1, .01))
  abline(v = c(0.33, 0.67), lwd = 4, lty = 2, col = "purple")

}
dev.off()
