\name{bltestHMM}
\alias{bltestHMM}
\title{
bltestHMM
}
\description{
combine bltest and HMM. Details coming soon.
}
\usage{
bltestHMM(untreated_ip, untreated_input, treated_ip, treated_input, 
                    untreated_ip_total, untreated_input_total, treated_ip_total, treated_input_total,
                    minimal_count_fdr =10,
                    Pi=matrix(c(0.9,0.1,0.1,0.9),byrow=TRUE, nrow=2),
                    delta=c(0.5,0.5),
                    pm=list(prob=c(0.1, 0.9)),
                    threshold_fdr=0.05, 
                    gene_label=rep(1,length(untreated_ip)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{untreated_ip}{ a vector of integers of n, which is the number of binding sites tested. Each element represents the number of reads fall into a binding site for the IP sample under untreated condition
%%     ~~Describe \code{untreated_ip} here~~
}
  \item{untreated_input}{ a vector of integers of n, which is the number of binding sites tested. Each element represents the number of reads fall into a binding site for the Input control sample under untreated condition
%%     ~~Describe \code{untreated_input} here~~
}
  \item{treated_ip}{ a vector of integers of n, which is the number of binding sites tested. Each element represents the number of reads fall into a binding site for the IP sample under treated condition
%%     ~~Describe \code{treated_ip} here~~
}
  \item{treated_input}{ a vector of integers of n, which is the number of binding sites tested. Each element represents the number of reads fall into a binding site for the Input control sample under treated condition
%%     ~~Describe \code{treated_input} here~~
}
  \item{untreated_ip_total}{ an integer, total number of reads for the IP sample under untreated condition
%%     ~~Describe \code{untreated_ip_total} here~~
}
  \item{untreated_input_total}{ an integer, total number of reads for the Input control sample under untreated condition
%%     ~~Describe \code{untreated_input_total} here~~
}
  \item{treated_ip_total}{ an integer, total number of reads for the IP sample under treated condition
%%     ~~Describe \code{treated_ip_total} here~~
}
  \item{treated_input_total}{  an integer, total number of reads for the Input control sample under treated condition
%%     ~~Describe \code{treated_input_total} here~~
}
  \item{minimal_count_fdr}{ an integer threshold, only the loci with reads more than this number are subjected for fdr calculation. default: 10
}
  \item{Pi}{ initial state transition probability matrix. default: matrix(c(0.9,0.1,0.1,0.9),byrow=TRUE, nrow=2)
}
  \item{delta}{ initial state probability. default: c(0.5,0.5)
}
  \item{pm}{ initial emission probability. default: c(0.1,0.9), i.e., the probability of observing a differential methylation state on an undifferential site is 0.1; and the probability for observing a differential mehtylation site on an differential methylation site is 0.9.
}

\item{threshold_fdr}{ the threshold used to convert bltest states to observations of HMM model. default: 0.05, i.e., if a site is differential methylated with significance level 0.05, than it is considered a differential methylation loci in HMM model.
}

\item{gene_label}{ if the counts are from multiple genes, then HMM is built on each gene. default: rep(1,length(untreated_ip)), i.e., the counts are all from the same gene.
}

}
\value{
The function returns a list of length 3, which contains the log(p-value), log(fdr) and log(fold change), respectively, from the test.
}

\details{
combine bltest and HMM.
The number of reads at the same locus for the aligned reads are
 counted by other packages, such as Rsamtools or HTseq-count.

}

\references{
Reference coming soon!
}
\author{ Yuchen Zhang  <490317018@qq.com>
%%  ~~who you are~~
}

\examples{
# example of using bltestHMM
untreated_ip = c(rpois(20,20),rpois(20,55),rpois(20,20), rpois(20,5)) 
# this is to simulate MeRIP-Seq data with 4 adjacent sections, 
# including 2 undifferential sites (US), and 2 differential methylation sites (DMS).
# the order of them is: US-DMS-US-DMS
untreated_input = rep(20,80)
treated_ip = rep(20,80)
treated_input = rep(20,80)
# sequencing depths
untreated_ip_total = 10^7
untreated_input_total = 10^7
treated_ip_total = 10^7
treated_input_total = 10^7
result <- bltestHMM(untreated_ip, untreated_input, 
                    treated_ip, treated_input, 
                    untreated_ip_total, untreated_input_total, 
                    treated_ip_total, treated_input_total)
# prediction
prediction <- (rbind( hmm=result$log.fdr,
  bltest=result$bltest_result$log.fdr) < log(0.05))
# comparison
comparison<- rbind(setting = c(rep(FALSE,20),rep(TRUE,20),rep(FALSE,20),rep(TRUE,20)), 
                   prediction);
hmm_perf <- 1-mean(abs(comparison[1,]-comparison[2,]))
bltest_perf <- 1-mean(abs(comparison[1,]-comparison[3,]))
result <- data.frame(blestHMM_performance = hmm_perf, blest_performance=bltest_perf)
print(result) # show the performance of two methods
# print(comparison) # show the deteailed prediction.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ tatistical Inference }
