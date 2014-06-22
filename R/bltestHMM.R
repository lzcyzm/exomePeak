

# Yuchen Zhang
# Email: jia.meng@hotmail.com
# differential analysis of RNA methylation MeRIP-Seq data
# based on bltest and HMM

bltestHMM <- function(untreated_ip, untreated_input, treated_ip, treated_input, 
                    untreated_ip_total, untreated_input_total, treated_ip_total, treated_input_total,
                    minimal_count_fdr =10,
                    Pi=matrix(c(0.9,0.1,0.1,0.9),byrow=TRUE, nrow=2),
                    delta=c(0.5,0.5),
                    pm=list(prob=c(0.1, 0.9)),
                    threshold_fdr=0.05, 
                    gene_label=rep(1,length(untreated_ip))
                    ) {
  # required library
  require(exomePeak)
  require(HiddenMarkov)
  
  # bltest
  result= bltest(untreated_ip, untreated_input, treated_ip, treated_input, untreated_ip_total,
                 untreated_input_total, treated_ip_total, treated_input_total,
                 minimal_count_fdr = minimal_count_fdr)
  labelbltest<- ( result$log.fdr<log(threshold_fdr) )
  
  # HMM section
  unique_gl = unique(gene_label)
  pos <- untreated_ip*0 + 1; # initialization the posterior probability
  for (i in 1:length(unique_gl)) {
    id = which(gene_label == unique_gl[i])
    labelbltest_i <- labelbltest[id]
    pos[id] <- .HMM_single(labelbltest_i,Pi,delta,pm)
  }
  
  # save HMM result
  log.fdr=log(p.adjust(pos, method = "fdr"))
  # adjust only sig reads count testing result
  m = untreated_ip + untreated_input + treated_ip + treated_input
  ID= which(m > minimal_count_fdr)
  log.fdr_sig=log(p.adjust(pos[ID], method = "fdr"))
  log.fdr[ID] = log.fdr_sig
  
  # remove infinity
  log.fdr=pmax(log.fdr,-1000)
  log.p=pmax(log(pos),-1000)
  
  # save result
  DIFF=list(log.fdr=log.fdr,log.p=log.p,log.fc=result$log.fc, bltest_result = result)
  return(DIFF)
}

# subfunction
.HMM_single <- function(labelbltest,Pi,delta,pm) {
  aa<-length(labelbltest)
  pn <- list(size=rep(1,aa))
  x <- dthmm(labelbltest, Pi, delta, "binom", pm, pn,discrete=TRUE)
  log <- capture.output({
    y <- BaumWelch(x);
  })
  pos <- y$u[,1]
  return(pos)
}

