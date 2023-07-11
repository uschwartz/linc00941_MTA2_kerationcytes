#!/usr/bin/env Rscript

library(stringr)
library(DESeq2)

csv<-list.files(pattern = "csv")
count.table<-read.delim(file=csv, skip=1)

# get features of name
id<-str_remove(colnames(count.table)[7:ncol(count.table)],".bam")
condition<-str_split_fixed(id,"_",6)[,1]
replica<-str_split_fixed(id,"_",6)[,2]
antibody<- str_split_fixed(replica,"\\.",2)[,1]


colnames(count.table)[7:ncol(count.table)]<-id
counts<-count.table[,7:ncol(count.table)]
rownames(counts)<-as.character(count.table$Geneid)

## prepare meta Data
#
meta.DATA<-data.frame(name=colnames(counts),
                      condition=condition,
                      replica,
                      antibody,
                      row.names = colnames(counts))

meta.DATA

## create DESeq obj
dds.obj <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = meta.DATA,
                                  design= ~condition)


## differential gene expression likelihhod ratio test
dds.diff <- DESeq(dds.obj)

res<-results(dds.diff,independentFiltering=FALSE )


samp<-str_remove(csv, ".csv")
sink(file = paste(samp,"summary.txt", sep="_"))
summary(res)
sink()


sig<-table(subset(res,padj < 0.1)$log2FoldChange>0)
save(sig, file = "sig.rda")
