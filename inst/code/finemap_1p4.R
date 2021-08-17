#!/usr/bin/env Rscript
library(dplyr)
library(readr)
library(magrittr)

#' FINEMAP I/O
write_finemap_sumstats <- function(beta, se, LD_file, n, k, prefix) {
  cfg = list(z=paste0(prefix,".z"),
             ld=LD_file,
             snp=paste0(prefix,".snp"),
             config=paste0(prefix,".config"),
             cred=paste0(prefix, ".cred"),
             k=paste0(prefix,".k"),
             log=paste0(prefix,".log"),
             meta=paste0(prefix,".master"))
  se = replace(se, se == 0, 'nan')
  z = data.frame(chromosome="chr", position=seq(1, length(beta)), allele1='nan', allele2='nan', maf='nan', beta, se)
  z = cbind(rsid=z$position, z)
  write.table(z,cfg$z,quote=F,col.names=T,row.names=F)
  if (!is.null(k)) {
      write.table(t(k),cfg$k,quote=F,col.names=F,row.names=F)
      write("z;ld;snp;config;cred;n_samples;k;log",file=cfg$meta)
      write(paste(cfg$z, cfg$ld, cfg$snp, cfg$config, cfg$cred, n, cfg$k, cfg$log, sep=";"),
        file=cfg$meta,append=TRUE)
  } else {
      write("z;ld;snp;config;cred;n_samples;log",file=cfg$meta)
      write(paste(cfg$z, cfg$ld, cfg$snp, cfg$config, cfg$cred, n, cfg$log, sep=";"),
            file=cfg$meta,append=TRUE)
  }
  return(cfg)
}

#' Run FINEMAP version 1.4
#' http://www.christianbenner.com
## FIXME: read the finemapr implementation for data sanity check.
## Can be useful as a general data sanity checker (in previous modules)

run_finemap <- function(beta, se, LD_file, n, k, args = "", prefix="data")
{
  cfg = write_finemap_sumstats(beta, se, LD_file, n, k, prefix)
  cmd = paste("finemap --sss --log", "--in-files", cfg$meta, args)
  dscrutils::run_cmd(cmd)
  cfg$log = paste0(cfg$log, "_sss")

  # read output tables
  snp = read.table(cfg$snp,header=TRUE,sep=" ")
  snp$snp = as.character(snp$rsid)

  snp = rank_snp(snp)
  # we add snp-prob for backwards-compatability with code that used this script with FINEMAP v1.1
  snp$prob = snp$snp_prob
  config = read.table(cfg$config,header=TRUE,sep=" ")

  # Only keep configurations with cumulative 95% probability
  # config = within(config, config_prob_cumsum <- cumsum(config_prob))
  # config = config[config$config_prob_cumsum <= 0.95,]

  # extract number of causal
  ncausal = finemap_extract_ncausal(cfg$log)
  return(list(snp=snp, set=config, ncausal=ncausal))
}

rank_snp <- function(snp) {
  snp <- arrange(snp, -prob) %>%
    mutate(
        rank = seq(1, n()),
        prob_cumsum = cumsum(prob) / sum(prob)) %>%
    select(rank, snp, prob, prob_cumsum, log10bf)
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

finemap_mvar <- function(beta, se, LD_file, n, k, args, prefix, parallel = FALSE) {
  if (is.null(dim(beta))) {
      beta = matrix(ncol=1,beta)
  }
  if (is.null(dim(se))) {
      se = matrix(ncol=1,se)
  }

  single_core = function(r) 
      run_finemap(beta[,r], se[,r], LD_file, n, k, args, 
                  prefix=paste0(prefix, '_condition_', r))
  if (parallel)
      return(parallel::mclapply(1:ncol(beta), function(r) single_core(r),
                                mc.cores = min(8, ncol(beta))))
  else
      return(lapply(1:ncol(beta), function(r) single_core(r)))
}

eval(parse(text=commandArgs(T)))
dat = readRDS(input)
sumstats = dat$sumstats
N = nrow(dat$data$X)
ld = tempfile(fileext = ".ld")
ld_mat = cor(dat$data$X)
ld_mat[is.na(ld_mat)] = 'nan'
write.table(ld_mat,ld,quote=F,col.names=F,row.names=F)
posterior = finemap_mvar(sumstats[1,,], sumstats[2,,],
                                        ld, N, k=NULL,
                                        args, prefix=tempfile(fileext = ".finemap"))
saveRDS(posterior, paste0(output, '.rds'))
