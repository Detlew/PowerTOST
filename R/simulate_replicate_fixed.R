#--------------------------------------------------------------------
# Functions to simulate replicate cross-over data
# fixed effects model (EMA)
# 
# Author: dlabes
#--------------------------------------------------------------------

# we are simulating without period effects and subject effects
# the last is astonishing for me but gives the same results as
# simulating subject effects with sD=0
# args:
# seqs = vector of sequences (f.i. c("TRTR", "RTRT")) 
# nseq = vector of number of subjects in sequences
# muR = mean of logvals for treatment "R" (arbitrary)
# ldiff = difference µT-µR (= log(GMR))
# s2wT, s2wR = within-subject variances for "T" or "R"
#
# prepare a dataframe with subject, sequence, period, logval
prep_data2 <- function(seqs, nseq, muR=log(10), ldiff, s2wT, s2wR)
{
  n_seqs <- length(seqs)
  ns <- integer(length=n_seqs)
  # if only one n is given repeat it for the sequences
  if (length(nseq)==1) ns <- rep(nseq, times=n_seqs)  else ns <- nseq
  data <- data.frame()
  for (i in seq_along(seqs))
  { 
    pers <- nchar(seqs[i])
    pervals <- rep(0,times=pers)
    tmt <- strsplit(seqs[i], split="")
    pervals <- matrix(0,nseq[i],pers)
    pdata <- data.frame(pervals)
    pdata$seqno    <- i
    pdata$sequence <- seqs[i]
    for(j in 1:ncol(pervals))
    {
      nam <- paste("P",j,sep="")
      names(pdata)[j] <- nam
    } 
    data <- rbind(data,pdata)
  } 
  # make the long form
  ldata <- reshape(data, v.names="logval", varying=list(1:(ncol(data)-2)),
                   timevar="period", idvar="subject", direction="long")
  # order by subject, period
  ldata     <- ldata[order(ldata$subject,ldata$period),]
  ldata$tmt <- substr(ldata$sequence,ldata$period,ldata$period)
  ldata     <- ldata[,c("subject","seqno","sequence","period","tmt","logval")]
  nvals     <- nrow(ldata)
  ldata$logval <- sim_data2_y(data_tmt=ldata$tmt, muR=muR, ldiff=ldiff, 
                              s2wT=s2wT, s2wR=s2wR)
               
  return(ldata)
}

# get new values for logval
sim_data2 <- function(data, muR=log(10), ldiff, s2wT, s2wR)
{
  nvals     <- nrow(data)
  data$logval <- ifelse(data$tmt=="T",
                        ldiff + rnorm(n=nvals, mean=0, sd=sqrt(s2wT)),
                        rnorm(n=nvals, mean=0, sd=sqrt(s2wR))) + muR
  data
}

# Bens modification to give only the vector of simulated log values back
# args:
# data_tmt= vector of treatments by subject and period coded as "T" or "R"
# muR = mean of logvals for treatment "R"
# ldiff = difference µT-µR
# s2wT, s2wR = within-subject variances for "T" or "R"
sim_data2_y <- function(data_tmt, muR=log(10), ldiff, s2wT, s2wR)
{
  nvals <- length(data_tmt)
  ifelse(data_tmt == "T", ldiff + rnorm(n=nvals, mean=0, sd=sqrt(s2wT)), 
         rnorm(n=nvals, mean=0, sd=sqrt(s2wR))) + muR
}
