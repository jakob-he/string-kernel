library(kebabs)

#preprocessing
sequences <- read.table("sequences.txt")
dnaseqs <- DNAStringSet(sequences[,1])
names(dnaseqs) <- paste("S", 1:length(dnaseqs), sep="")
mot <- motifKernel(c("TCAGCA","TGCTGA","ATGCA.A","T.TGCAT"), normalized=FALSE)

## applying the Kernel
start_time <- Sys.time()
km <- mot(dnaseqs)
end_time <- Sys.time()

end_time-start_time

