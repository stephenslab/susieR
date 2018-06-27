#!/usr/bin/env Rscript
library(dplyr)
library(readr)
library(magrittr)

#' CAVIAR I/O
write_caviar_sumstats <- function(z, prefix) {
  cfg = list(z=paste0(prefix,".z"),
             set=paste0(prefix,"_set"),
             post=paste0(prefix,"_post"),
             log=paste0(prefix,".log"))
  write.table(z,cfg$z,quote=F,col.names=F)
  return(cfg)
}

#' Run CAVIAR
#' https://github.com/fhormoz/caviar

run_caviar <- function(z, LD_file, args = "", prefix="data")
{
  cfg = write_caviar_sumstats(z, prefix)
  cmd = paste("CAVIAR", "-z", cfg$z, "-l", LD_file, "-o", prefix, args)
  dscrutils::run_cmd(cmd)
  if(!all(file.exists(cfg$post, cfg$set, cfg$log))) {
      stop("Cannot find one of the post, set, and log files")
  }
  
  log <- readLines(cfg$log)

  # read output tables
  snp <- read.delim(cfg$post)  
  stopifnot(ncol(snp) == 3)
  names(snp) <- c("snp", "snp_prob_set", "snp_prob")
  snp$snp <- as.character(snp$snp)
  snp <- rank_snp(snp)

  # `set` of snps
  set <- readLines(cfg$set)
  set_ordered <- left_join(data_frame(snp = set), snp, by = "snp") %>% 
    arrange(rank) %$% snp
  return(list(snp=snp, set=set_ordered))
}

rank_snp <- function(snp) {
  snp <- arrange(snp, -snp_prob) %>%
    mutate(
        rank = seq(1, n()),
        snp_prob_cumsum = cumsum(snp_prob) / sum(snp_prob)) %>%
    select(rank, snp, snp_prob, snp_prob_cumsum, snp_prob_set)
  return(snp)    
}

finemap_mcaviar <- function(zscore, LD_file, args, prefix) {
  if (is.null(dim(zscore))) {
      zscore = matrix(ncol=1,zscore)
  }
  return(parallel::mclapply(1:ncol(zscore), function(r)
          run_caviar(zscore[,r], LD_file, args, 
                     paste0(prefix, '_condition_', r)), 
                            mc.cores = min(8, ncol(zscore))))
}

eval(parse(text=commandArgs(T)))
dat = readRDS(input)
sumstats = dat$sumstats
ld = tempfile(fileext = ".ld")
write.table(cor(dat$data$X),ld,quote=F,col.names=F,row.names=F)
posterior = finemap_mcaviar(sumstats[1,,] / sumstats[2,,],
                            ld,
                            args, prefix=tempfile(fileext = ".caviar"))
saveRDS(posterior, paste0(output, '.rds'))
