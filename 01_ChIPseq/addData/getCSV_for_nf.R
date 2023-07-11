pathSET<-"/Volumes/PromisePegasus/_Service_/S026_siChIP_E2F6_MTA2/data/FC1_221116/"

setwd(pathSET)
library(stringr)

fwd<-grep("R1_001.fastq.gz",list.files(), value = T)
rev<-grep("R2_001.fastq.gz",list.files(), value = T)

##extract information

name.split<-str_replace(fwd,"Input_", "Input-") %>%
    str_split_fixed( "_",8)

names<-paste(name.split[,3], name.split[,4], sep="_")
replicate<-name.split[,4]
antibody<-str_split_fixed(name.split[,4],"-",2)[,1]

input.match<-str_split_fixed(name.split[,4],"-",2)[,2]
input<-(name.split[match(antibody, input.match),4])
input[is.na(input)]<-""



df<-data.frame(name=names,
               replicate,
               path_fwdReads=paste0(pathSET,fwd),
               path_revReads=paste0(pathSET,rev),
               antibody,
               input)

write.csv(df, file="../../script/addData/fastqIN.csv",
          quote=F, row.names = F)
