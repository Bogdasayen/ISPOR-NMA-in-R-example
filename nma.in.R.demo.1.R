# Example of NMA in R
# Howard Thom 11-October-2020

# Part of ISPOR EU 2020 Workshop titled:
# WHY R WE BOTHERING? THE WHY AND HOW OF USING R FOR INTEGRATED ANALYSIS AND MODELLING

# Data based on "Intra-Cavity Lavage and Wound Irrigation for Prevention of Surgical Site Infection: Systematic Review and Network Meta-Analysis"
# Thom H, Norman G, Welton NJ, Crosbie, EJ, Blazeby J, Dumville JC.
# Surg Infect (Larchmt) 2020 Apr 29. doi: 10.1089/sur.2019.318

# This script can integrated with a CEA analysis as MCMC sampes are exported by bugs() function

# Load data (processed in another R script)
load("irrigation.data.v2.rda")

library(R2OpenBUGS)


# NICE DSU TSD 2
source("model.binomial.logistic.1.R")
# Chaimani 2013: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0076654
source("mtm.networkplot.fun.R")
# Load script to generate odds ratio cross table
source("cross.or.1.R")

# Draw the network plot
x<-bugs.data$t
for(i in 1:dim(x)[1])
{
  x[i,]<-i
}
t1<-c(x[!is.na(bugs.data$t)])
t2<-c(bugs.data$t[!is.na(bugs.data$t)])
jpeg(width=550,file="network.plot.jpg",sep="")
mtm.networkplot.fun(t1=t1,t2=t2,percomparison=FALSE,
                    nameoftreatments=paste(1:bugs.data$nt,t.names))
dev.off()

# Number of MCMC chains and samples
n.chains<-2		 # 2
num.sims=30000*n.chains 
burn.in=30000*n.chains	

# Create initial values for MCMC simulation
inits1<-list(d=c(NA,rep(1,bugs.data$nt-1)),mu=rep(0.5,bugs.data$ns))
inits2<-list(d=c(NA,rep(0.5,bugs.data$nt-1)),mu=rep(0.25,bugs.data$ns))
bugs.inits<-list(inits1,inits2)

# Call OpenBUGS
bugs.object.fe<-bugs(data=bugs.data,inits=bugs.inits,
                     parameters.to.save=c("d","or","rk","totresdev"),
                     model=model.binomial.logistic.fe,clearWD=TRUE,summary.only=FALSE,
                     n.iter=(num.sims+burn.in),n.burnin=burn.in,
                     n.chains=n.chains,bugs.seed=1,debug=TRUE)

# Format the odds ratios
cross.or.fe<-cross.or(bugs.object=bugs.object.fe,t.names=t.names,med=TRUE)
write.csv(x=cross.or.fe,file="cross.or.fe.csv")






