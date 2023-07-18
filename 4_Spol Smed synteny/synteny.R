suppressMessages(library(tidyverse))

# parse GMAP output tables
# S. mediterranea

cmd = "grep -c '>SM' Projects/Planaria/Synteny/gmap_output/SMEST_on_SMESG.gmap"
inp = system(cmd, intern = TRUE) %>% as.integer()
ll <- vector(mode = "list", length = inp)

con = base::file("Projects/Planaria/Synteny/gmap_output/SMEST_on_SMESG.gmap", "r")
co = 0
while ( TRUE ) {
  line = readLines(con, n = 1)
  if ( length(line) == 0 ) {
    break
  }
  if ( grepl(">SM", line) ) {
    chimarea = FALSE
    co = co + 1
    name = gsub(">", "", line)
    names(ll)[co] = gsub(" .*", "", name)
    line = readLines(con, n = 1)
    pathno = substr(line, 8, 23)
    if (grepl("Possible", pathno) ) {
      pathno = 2
      chimarea = TRUE
    } else {
      pathno = min(as.integer(gsub("\\):.*", "", pathno)), 5)
    }
    if (pathno > 0) {
      paths = matrix(NA, nrow = pathno, ncol = 5) %>% as.data.frame()
      names(paths) <- c("Scaffold", "Start", "End", "Cov", "Ori")
      line2 = readLines(con, n = pathno*10+4)
      if (sum(grepl("Non-intron", line2)) > 0) line2 <- line2[-which(grepl("Non-intron", line2))]
      if (sum(grepl("Translation", line2)) > 0) line2 <- line2[-which(grepl("Translation", line2))]
      pos_lines = line2[seq(4, pathno*10, 10)]
      cov_lines = line2[seq(6, pathno*10, 10)]
      line2[grepl("strand", line2)] -> ori
      gsub(".*\\(", "", ori) %>% gsub(" strand\\)", "", .) -> paths$Ori
      paths$Cov <- strsplit(cov_lines, split = " ") %>% sapply("[", 6) %>% as.numeric()
      paths[,1:3] <- pos_lines %>% strsplit(split = " ") %>% sapply("[", 6) %>% 
        gsub(":", " ", .) %>% gsub("\\.\\.", " ", .) %>% 
        gsub(",", "", .) %>% strsplit(" ") %>% do.call(rbind, .)
      } 
     
      
      
    } else {
      paths = NA
    }
    ll[[co]]$Path_count = pathno
    ll[[co]]$Paths <-  paths
    ll[[co]]$Chimaera <- chimarea
    if (co %% 1000 == 0) print(co)
  }
}
close(con)
smed <- ll

# S. polychroa

cmd = "grep -c '>dd' /Projects/Planaria/Synteny/gmap_output/ddSpolv4_on_bmSpolg1.gmap"
inp = system(cmd, intern = TRUE) %>% as.integer()
ll <- vector(mode = "list", length = inp)

con = base::file("/Projects/Planaria/Synteny/gmap_output/ddSpolv4_on_bmSpolg1.gmap", "r")
co = 0
while ( TRUE ) {
  line = readLines(con, n = 1)
  if ( length(line) == 0 ) {
    break
  }
  if ( grepl(">dd", line) ) {
    chimarea = FALSE
    co = co + 1
    name = gsub(">", "", line)
    names(ll)[co] = gsub(" .*", "", name)
    line = readLines(con, n = 1)
    pathno = substr(line, 8, 23)
    if (grepl("Possible", pathno) ) {
      pathno = 2
      chimarea = TRUE
    } else {
      pathno = min(as.integer(gsub("\\):.*", "", pathno)), 5)
    }
    if (pathno > 0) {
      paths = matrix(NA, nrow = pathno, ncol = 5) %>% as.data.frame()
      names(paths) <- c("Scaffold", "Start", "End", "Cov", "Ori")
      line2 = readLines(con, n = pathno*10+4)
      if (sum(grepl("Non-intron", line2)) > 0) line2 <- line2[-which(grepl("Non-intron", line2))]
      if (sum(grepl("Translation", line2)) > 0) line2 <- line2[-which(grepl("Translation", line2))]
      pos_lines = line2[seq(4, pathno*10, 10)]
      cov_lines = line2[seq(6, pathno*10, 10)]
      line2[grepl("strand", line2)] -> ori
      gsub(".*\\(", "", ori) %>% gsub(" strand\\)", "", .) -> paths$Ori
      paths$Cov <- strsplit(cov_lines, split = " ") %>% sapply("[", 6) %>% as.numeric()
      k <- pos_lines %>% strsplit(split = " ") %>% sapply("[", 6) %>% 
        gsub(":", " ", .) %>% gsub("\\.\\.", " ", .) %>% 
        gsub(",", "", .) %>% strsplit(" ")
      for (ik in seq_along(k)) {
        if (length(k[[ik]]) == 4) {
          k[[ik]] <- c(paste0(k[[ik]][1], ":", k[[ik]][2]), k[[ik]][3], k[[ik]][4])
        }
      }

      paths[,1:3] <- do.call(rbind.data.frame, k)
    } else {
      paths = NA
    }

    ll[[co]]$Path_count = pathno
    ll[[co]]$Paths <-  paths
    ll[[co]]$Chimaera <- chimarea
    if (co %% 1000 == 0) print(co)
  }
}
close(con)
spol <- ll

# read in orthology information parsed previously from PlanMine
ortho_tab <- read.delim("/Projects/Planaria/Synteny/spol_smed_ortho_tab.txt")

species_list <- c("spol", "smed")

# basic function for getting all mapped transcripts on a scaffold together with the homology links 
get_genes_for_scaffold <- function(species, scaffold) {
  get(species)[sapply(get(species), function(x) length(x$Paths) == 5)] -> ret
  ret[sapply(ret, function(x) any(x$Paths$Scaffold == scaffold))] -> ret
  do.call(rbind.data.frame, sapply(ret, "[", 2)) -> ret
  gsub("\\.Paths", "", rownames(ret)) -> ret$transcript_id
  ret <- ret[which(ret$Scaffold == scaffold),]
  rownames(ret) <- NULL
  ret$End <- as.numeric(ret$End)
  ret$Start <- as.numeric(ret$Start)
  ret <- ret[order(ret$Start),]
  if (species == "spol") {
    ret$Spol <- gsub("([0-9]*)(\\.[0-9]*)$", "\\1", ret$transcript_id)
    ret$tmp <- gsub("([0-9]*)(\\.[0-9]*)$", "\\1", ret$transcript_id)
    merge(ret, ortho_tab, all.x = TRUE) -> ret
    ret$tmp2 <- ret$SMEST
    smed[ret$tmp2[!is.na(ret$tmp2)]] %>% sapply("[", 2) %>% do.call(rbind.data.frame, .) -> k
  }  else {
    ret$SMEST <- gsub("(\\.[0-9]*)(\\.[0-9]*)$", "\\1", ret$transcript_id)
    ret$tmp <- gsub("(\\.[0-9]*)(\\.[0-9]*)$", "\\1", ret$transcript_id)
    merge(ret, ortho_tab, all.x = TRUE) -> ret
    ret$tmp2 <- ret$Spol
    spol[ret$tmp2[!is.na(ret$tmp2)]] %>% sapply("[", 2) %>% do.call(rbind.data.frame, .) -> k
  } 
  if (nrow(k) == 0) return("zero_homologs_on_scaffold")
  k$tmp2 <- gsub("\\.Paths.*", "", rownames(k))
  rownames(k) <- NULL
  merge(ret, k, by = "tmp2", all = TRUE) -> ret
  ret$Ori.x = ifelse(is.na(ret$tmp2), ifelse(ret$num_Spol_homologs == 0, "no_Spol_ortholog_in_PlanMine", "Multiple_Spol_orthologs"), ret$Ori.x)
  ret %>% 
    mutate(Scaffold.y = ifelse(is.na(Scaffold.y) & !is.na(Spol) & !is.na(SMEST), "not_mapped", Scaffold.y),
           Ori.y = ifelse(is.na(Ori.y) & Scaffold.y == "not_mapped", "not_mapped", Ori.y)) %>% 
    arrange(Start.x) %>%
    ungroup()%>%
    dplyr::select(-tmp, -tmp2) %>% filter(!is.na(Ori.x)) %>%
    mutate(Len = as.numeric(End.x) -as.numeric(Start.x), 
           inside = (lag(Len) < 5e4 & as.numeric(lag(End.x)) > as.numeric(Start.x)),
           inside = ifelse(is.na(inside), FALSE, inside)) %>%
    as.data.frame() -> ret
  return(ret)
}

# summarizing table for each Spol scaffold: linked Smed site count & gene count 
g <- spol[sapply(spol, function(x) all(x$Paths != "NA"))]
sapply(g, "[[", 2) %>% do.call(rbind.data.frame, .) -> g
unique(g$Scaffold) ->g
spol_scaf_tab <- data.frame(Scaffold = g, Smed_sites = NA, Homolog_count = NA, Filtered = NA)
for (i in g) {
  tmp <- get_genes_for_scaffold("spol", i) 
  if (class(tmp) == "character") {
    spol_scaf_tab[which(spol_scaf_tab$Scaffold == i), "Smed_sites"] <- 0
  } else {
    tmp2 <- filter(tmp, n_ratio < .1, !inside, Len > 300)
    spol_scaf_tab[which(spol_scaf_tab$Scaffold == i), "Smed_sites"] <- length(table(unique(tmp2$Scaffold.y)))
    spol_scaf_tab[which(spol_scaf_tab$Scaffold == i), "Homolog_count"] <- nrow(tmp2)
    spol_scaf_tab[which(spol_scaf_tab$Scaffold == i), "Filtered"] <- nrow(tmp) - nrow(tmp2)
  }
  
  print(spol_scaf_tab[which(spol_scaf_tab$Scaffold == i),])
}

# get scaffolds with 2 Smed links
filter(spol_scaf_tab, Smed_sites == 2)$Scaffold -> h
scaf_2hit <- data.frame(Scaffold = h, A_count = NA, B_count = NA, Filtered = NA)
for (i in h) {
  tmp <- get_genes_for_scaffold("spol", i) 
  if (class(tmp) != "character") {
    tmp2 <- filter(tmp, n_ratio < .1, !inside, Len > 300)
    scaf_2hit[which(scaf_2hit$Scaffold == i), "A_count"] <- sort(table(tmp2$Scaffold.y))[2]
    scaf_2hit[which(scaf_2hit$Scaffold == i), "B_count"] <- sort(table(tmp2$Scaffold.y))[1]
    scaf_2hit[which(scaf_2hit$Scaffold == i), "Filtered"] <- nrow(tmp) - nrow(tmp2)
  }
  print(scaf_2hit[which(scaf_2hit$Scaffold == i),])
}

# get scaffolds with 3 Smed links
filter(spol_scaf_tab, Smed_sites == 3)$Scaffold -> h
scaf_3hit <- data.frame(Scaffold = h, A_count = NA, B_count = NA, C_count = NA, Filtered = NA)
for (i in h) {
  tmp <- get_genes_for_scaffold("spol", i) 
  if (class(tmp) != "character") {
    tmp2 <- filter(tmp, n_ratio < .1, !inside, Len > 300)
    scaf_3hit[which(scaf_3hit$Scaffold == i), "A_count"] <- sort(table(tmp2$Scaffold.y))[3]
    scaf_3hit[which(scaf_3hit$Scaffold == i), "B_count"] <- sort(table(tmp2$Scaffold.y))[2]
    scaf_3hit[which(scaf_3hit$Scaffold == i), "C_count"] <- sort(table(tmp2$Scaffold.y))[1]
    scaf_3hit[which(scaf_3hit$Scaffold == i), "Filtered"] <- nrow(tmp) - nrow(tmp2)
  }
  print(scaf_3hit[which(scaf_3hit$Scaffold == i),])
}

# get scaffolds with 4 Smed links
filter(spol_scaf_tab, Smed_sites == 4)$Scaffold -> h
scaf_4hit <- data.frame(Scaffold = h, A_count = NA, B_count = NA, C_count = NA, D_count = NA, Filtered = NA)
for (i in h) {
  tmp <- get_genes_for_scaffold("spol", i) 
  if (class(tmp) != "character") {
    tmp2 <- filter(tmp, n_ratio < .1, !inside, Len > 300)
    scaf_4hit[which(scaf_4hit$Scaffold == i), "A_count"] <- sort(table(tmp2$Scaffold.y))[4]
    scaf_4hit[which(scaf_4hit$Scaffold == i), "B_count"] <- sort(table(tmp2$Scaffold.y))[3]
    scaf_4hit[which(scaf_4hit$Scaffold == i), "C_count"] <- sort(table(tmp2$Scaffold.y))[2]
    scaf_4hit[which(scaf_4hit$Scaffold == i), "D_count"] <- sort(table(tmp2$Scaffold.y))[1]
    scaf_4hit[which(scaf_4hit$Scaffold == i), "Filtered"] <- nrow(tmp) - nrow(tmp2)
  }
  print(scaf_4hit[which(scaf_4hit$Scaffold == i),])
}


