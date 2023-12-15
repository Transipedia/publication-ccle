#!/usr/bin/env Rscript

library("Biostrings")
require(dplyr)

# transform a fasta file into dataframe
fasta2dataframe <- function(fasta){
	s <- readDNAStringSet(fasta)
	id <- names(s)
	seq <- c()
	#print(toString(s[[1]]))
	for (i in 1:length(s)){
		seq[i] <- toString(s[[i]])
	}
	id_seq=data.frame(id,seq)
	return(id_seq)
}

# extract substrings of length k from a sequence
extractSubstrings <- function(long_string, k) {
	n <- nchar(long_string)
	indexes <- 1:(n - k + 1)
	substrings <- unlist(lapply(indexes, function(i) substring(long_string, i, i + k - 1)))
	return(substrings)
}

# get one sequence complexity based on the 3-mer diversity
# the computed score is between 0 and 1
# fa = one sequence (kmers or contigs)
complexity <- function(fa, k=3) {
	sub <- extractSubstrings(fa, k)
	div=length(unique(sub))
	nmers=length(sub)
	#print(head(paste(fa, div/nmers))) #div/nmers = complexity
	return(c(fa, div/nmers))
}

# complexity from Daniel :
#for (i in 1:length(fa)) {
#  si=fa[i]
#  sub <- extract_substrings(si, k)
#  div=length(unique(sub))
#  nmers=length(sub)
#  print(paste(si, div/nmers) #div/nmers est la complexité.
#}

# write a dataframe to fasta format
writeFasta <-function(data, filename){
	fastaLines = c()
	for (rowNum in 1:nrow(data)){
		fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"id"], sep = "")))
		fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
	}
	fileConn<-file(filename)
	writeLines(fastaLines, fileConn)
	close(fileConn)
}

######### -------------------- Main -------------------- #########

# check the input arguments
## args[1] = fasta file; in our case = kmerator output fasta file
## args[2] = file name for output full table (with kmers name, sequence and complexity score);
## args[3] = file name for output fasta file (containing only the kmers that pass the complexity score threshold)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else {
	if (length(args)==1) {
  		# default output files
  		args[2] = "out.tsv" # output table
  		args[3] = "out.fa" # output fasta file
	}
	else if (length(args)==2) {
		args[3] = "out.fa"
	}
}

# get the kmer sequences into tab format from input fasta file

df_fasta <- fasta2dataframe(args[1])

# get the complexity score for each kmer and write the output full table

df <- data.frame(id=df_fasta$id,do.call(rbind,lapply(df_fasta$seq, complexity)))
colnames(df)[2:3] <- c("seq","score")
write.table(df, file = args[2], quote=FALSE, sep = "\t", row.names=FALSE)

# generate a new fasta file with only the kmers with complexity above threshold

threshold <- 0.55
df_filter <- df[which(df$score>threshold),]
#print(head(df_filter))

writeFasta(df_filter, args[3])