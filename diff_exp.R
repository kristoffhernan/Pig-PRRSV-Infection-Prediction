library(tidyverse)
library(DESeq2)



# setwd("C:/GenoTwin/project_9")


# read in metadata
meta <- read.table("./metadata.csv", sep="\t", header=TRUE)

# read in rna data
df<-read.table("./allfile_counts.csv", sep="\t", header=TRUE)

# remove extra column
df <- subset(df, select = -c(no_feature,ambiguous,too_low_aQual,not_aligned,alignment_not_unique, dpi) )
df <- subset(df, select = -c(X45360240,X45360241,X45360242,X45360243,X45360244,X45360245,X45360246,X45360247) )

# drop rows with NA
df <- drop_na(df)

rownames(df) <- df$X
df <- subset(df, select = -c(X) )

df <- t(df)


# only keep metadata for samples we have data for
colnames <- colnames(df)
meta <- subset(meta, X %in% colnames)


# DESEQ

dds <- DESeqDataSetFromMatrix(countData = df,
                              colData = meta,
                              design = ~ postInf)


DE <- DESeq(dds)

# deseq normalization
normalized_counts <- counts(DE, normalized=TRUE)
write.table(normalized_counts, file="./normalized_counts.txt", sep="\t", quote=F, col.names=NA)

# results
res <- results(DE)
write.table(res, file="./deFull.txt", sep="\t", quote=F, col.names=NA)


# only keep with adjusted p value less than 0.1
resSig <- res[ which(res$padj < 0.1 ), ]

# most over and under expressed
head( resSig[ order( resSig$log2FoldChange ), ] )
tail( resSig[ order( resSig$log2FoldChange ), ] )

write.table(resSig, file="./deGenes.txt", sep="\t", quote=F, col.names=NA)









