#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
path <-  args[1]

library(dada2)
marker <- basename(path)


fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_R1"), `[`, 1)

# jpeg(file=file.path(path, paste(marker, ".jpg", sep = "")) )
# plotQualityProfile(fnFs[1:2])
# dev.off()


filtFs <- file.path(path,  paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path,  paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, minLen=100,
                     compress=TRUE, multithread=TRUE, matchIDs=TRUE, verbose=TRUE)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap=20, maxMismatch=1, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)
uniquesToFasta(getUniques(seqtab), fout=file.path(path, paste(marker, "uniqueSeqs.fasta", sep="")), ids=paste0("Seq", seq(length(getUniques(seqtab)))))
write.csv(t(seqtab), file.path(path, paste(marker, ".seqtab.csv", sep="")))
save.image(file=file.path(path,paste(marker, "uniqueSeqs.Dada2.RData", sep="")))



