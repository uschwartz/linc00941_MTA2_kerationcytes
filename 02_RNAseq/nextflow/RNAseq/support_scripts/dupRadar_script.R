
args   <- commandArgs(TRUE)



library("dupRadar")

bamDuprm<-args[1]
gtf <- args[2]
stranded <- as.integer(args[3])
paired   <-  ifelse(args[4],TRUE, FALSE)
threads  <- 4 

dm <- analyzeDuprates(bamDuprm,gtf,stranded,paired,threads)

pdf("density_plot_coding.pdf")
  duprateExpDensPlot(DupMat=dm)
dev.off()

write.table(dm, file="duplication_table_coding.txt", sep="\t", quote=F,
            row.names=F)

pdf("boxplot_dups_coding.pdf")
  duprateExpBoxplot(DupMat=dm) 
dev.off()

pdf("smoothed_density_plot_coding.pdf")
  duprateExpPlot(DupMat=dm)
dev.off()

pdf("expression_bins_distribution_coding.pdf")
  readcountExpBoxplot(DupMat=dm)
dev.off()
