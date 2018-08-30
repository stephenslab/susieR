#!/usr/bin/env Rscript
library(dplyr)
library(readr)
library(magrittr)

#' FINEMAP I/O
write_finemap_sumstats <- function(z, LD_file, n, k, prefix) {
  cfg = list(z=paste0(prefix,".z"),
             ld=LD_file,
             snp=paste0(prefix,".snp"),
             config=paste0(prefix,".config"),
             k=paste0(prefix,".k"),
             log=paste0(prefix,".log"),
             meta=paste0(prefix,".master"))
  write.table(z,cfg$z,quote=F,col.names=F)
  if (!is.null(k)) {
      write.table(t(k),cfg$k,quote=F,col.names=F,row.names=F)
      write("z;ld;snp;config;k;log;n-ind",file=cfg$meta)
      write(paste(cfg$z, cfg$ld, cfg$snp, cfg$config, cfg$k, cfg$log, n, sep=";"),
        file=cfg$meta,append=TRUE)
  } else {
      write("z;ld;snp;config;log;n-ind",file=cfg$meta)
      write(paste(cfg$z, cfg$ld, cfg$snp, cfg$config, cfg$log, n, sep=";"),
            file=cfg$meta,append=TRUE)
  }
  return(cfg)
}

#' Run FINEMAP version 1.1
#' http://www.christianbenner.com
## FIXME: read the finemapr implementation for data sanity check.
## Can be useful as a general data sanity checker (in previous modules)

run_finemap <- function(z, LD_file, n, k, args = "", prefix="data")
{
  cfg = write_finemap_sumstats(z, LD_file, n, k, prefix)
  cmd = paste("finemap --sss --log", "--in-files", cfg$meta, args)
  dscrutils::run_cmd(cmd)

  # read output tables
  snp = read.table(cfg$snp,header=TRUE,sep=" ")
  snp$snp = as.character(snp$snp)

  snp = rank_snp(snp)
  config = read.table(cfg$config,header=TRUE,sep=" ")

  # Only keep configurations with cumulative 95% probability
  # config = within(config, config_prob_cumsum <- cumsum(config_prob))
  # config = config[config$config_prob_cumsum <= 0.95,]

  # extract number of causal
  ncausal = finemap_extract_ncausal(cfg$log)
  return(list(snp=snp, set=config, ncausal=ncausal))
}

rank_snp <- function(snp) {
  snp <- arrange(snp, -snp_prob) %>%
    mutate(
        rank = seq(1, n()),
        snp_prob_cumsum = cumsum(snp_prob) / sum(snp_prob)) %>%
    select(rank, snp, snp_prob, snp_prob_cumsum, snp_log10bf)
  return(snp)
}

finemap_extract_ncausal <- function(logfile)
{
  lines <- grep("->", readLines(logfile), value = TRUE)
  lines <- gsub("\\(|\\)|>", "", lines)
  splits <- strsplit(lines, "\\s+")
  tab <- data.frame(
    ncausal_num = sapply(splits, function(x) as.integer(x[2])),
    ncausal_prob = sapply(splits, function(x) as.double(x[4])))
  tab <- mutate(tab, type = ifelse(duplicated(ncausal_num), "post", "prior"))
  return(tab)
}

finemap_mvar <- function(zscore, LD_file, n, k, args, prefix, parallel = FALSE) {
  if (is.null(dim(zscore))) {
      zscore = matrix(ncol=1,zscore)
  }

  single_core = function(r) 
      run_finemap(zscore[,r], LD_file, n, k, args, 
                  paste0(prefix, '_condition_', r))
  if (parallel)
      return(parallel::mclapply(1:ncol(zscore), function(r) single_core(r),
                                mc.cores = min(8, ncol(zscore))))
  else
      return(lapply(1:ncol(zscore), function(r) single_core(r)))
}

eval(parse(text=commandArgs(T)))
dat = readRDS(input)
sumstats = dat$sumstats
N = nrow(dat$data$X)
ld = tempfile(fileext = ".ld")
write.table(cor(dat$data$X),ld,quote=F,col.names=F,row.names=F)
posterior = finemap_mvar(sumstats[1,,] / sumstats[2,,],
                                        ld, N, k=NULL,
                                        args, prefix=tempfile(fileext = ".finemap"))
saveRDS(posterior, paste0(output, '.rds'))
