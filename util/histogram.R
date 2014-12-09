#
# Plots a fancy oversequencing histogram from tables pre-computed by Histogram util
#
# usage:
# $cd histogram_output_dir/
# $RScript histogram.R
#

require(ggplot2); require(reshape)

build_df <- function(stat, units, sweep_df = NULL) {
      df<- read.table(paste(stat, units, ".txt", sep= ""), header=FALSE, comment="", sep = "\t")
      x<-df[1,3:ncol(df)]
      id<-df[2:nrow(df), 1] 
      
      if (is.null(sweep_df)) {
         sweep_df <- df
      }
      
      df[2:nrow(df),3:ncol(df)] <- sweep(df[2:nrow(df),3:ncol(df)], 1, rowSums(sweep_df[2:nrow(df),3:ncol(df)]), FUN = "/")
      
      df.m <- melt(df[2:nrow(df), c(1,3:ncol(df))])
      df.m$x <- as.vector(t(x[df.m$variable]))
      df.m$s <- rep(stat, nrow(df.m))
      list(sweep_df, df.m)
   }
   
for (units in c("", "-units")) {  
   dfo <- build_df("overseq", units)
   dfc <- build_df("collision1", units, dfo[[1]])     
   df <- rbind(dfo[[2]], dfc[[2]])
   
   pdf(paste("histogram", units, ".pdf", sep= ""))
   
   print(
      ggplot(df, aes(x=x,y=value))+geom_smooth()+
         scale_x_log10(name = "MIG size, reads", expand=c(0,0), limits=c(1, 10000), breaks = c(1:10,100,1000,10000), labels = c("1", rep("", 8), 10, 100, 1000, 10000), oob=scales::rescale_none)+
         scale_y_continuous(name = "",expand=c(0,0), limits = c(0, max(df$value)), oob=scales::rescale_none)+theme_bw()+theme(panel.grid.minor = element_blank()) + facet_grid(s~.)
   )
   
   dev.off()  
}