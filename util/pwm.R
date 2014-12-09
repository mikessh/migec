#
# Plots position-weight matrices for UMI sequences based on data precomputed by Histogram util
#
# usage:
# $cd histogram_output_dir/
# $RScript pwm.R
#

require(seqLogo) # available @ bioconductor

# read in
logo <- function(prefix) {
   df<-read.table(paste(prefix,".txt",sep=""))
   rownames(df)<-df[,1]
   df<-df[,2:ncol(df)]

   # build seqlogo
   df.p <- makePWM(df)
   pdf(paste(prefix,".pdf",sep=""))
   seqLogo(df.p, ic.scale = F)
   dev.off()

   con <- file(paste(prefix,"-stats.txt",sep=""))
   sink(con, append=TRUE)
   sink(con, append=TRUE, type="output")

   # stats
   print(paste("IC =",sum(df.p@ic)))
   #entropy
   h<-sum(2-df.p@ic)
   print(paste("H =",h))
   #correlation with pos
   print(paste("R(Hi,i) =", cor(-df.p@ic,1:length(df.p@ic))))
   #number of variants
   print(paste("N_obs =",2^h))
   #theoretical number of variants
   print(paste("N_exp =",2^(2*ncol(df))))

   sink()
}

logo("pwm-summary")
logo("pwm-summary-units")