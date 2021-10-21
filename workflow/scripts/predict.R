
library(seqinr)
library(pracma)
library(monomvn)
source("workflow/scripts/seqBlasso.R")

mafft_file <- as.character(snakemake@input["mafft"])
blasso_file <- as.character(snakemake@input["blasso"])
pos_file   <- as.character(snakemake@input["positions"])
amino_feature_file <- as.character(snakemake@input["amino_feature_file"])
output_file <- as.character(snakemake@output)

position <- as.matrix(read.table(pos_file))
# blasso object
load(blasso_file)
train_range <- 1:length(blasso$sample_names)

###############################################
# Inspired by workflow/scripts/seqBlasso.R    #
###############################################

features <- read.table(amino_feature_file, row.names = 1, header = T, sep = ",")
aa_feat <- features / outer(rep(1, nrow(features)), apply(abs(features), 2, max))

fasta <- load_fasta_data(mafft_file, position)
sample_names <- names(fasta$fasta)

X <- fasta$seq.bin
n_all <- nrow(X)
aseq_sub <- fasta$fasta.seq
seq_len <- fasta$seq.len

X <- do.call("cbind", lapply(1:ncol(aa_feat), function(i) {
	matrix(aa_feat[aseq_sub, i], n_all, seq_len)
}))
for (i in 1:nrow(X)) {
	X[i,] <- c(t(matrix(X[i,], ncol(aseq_sub), ncol(aa_feat))))
}
wavelength <- cbind(1,X) %*% blasso$beta_mean

###############################################

pred_wavelen_sample <- data.frame(id = sample_names[-train_range], wavelength = wavelength[-train_range])
write.table(pred_wavelen_sample, file = output_file, quote = F, row.names = F, sep = "\t")
