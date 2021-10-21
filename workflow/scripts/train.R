
library(seqinr)
library(pracma)
library(monomvn)

source("workflow/scripts/seqBlasso.R")

wavelength <- as.character(snakemake@input["wavelength"])
positions  <- as.character(snakemake@input["positions"])
sequences  <- as.character(snakemake@input["sequences"])
zero       <- as.character(snakemake@input["zero"])
amino_acid_feature <- as.character(snakemake@input["amino_acid_feature"])

n_sample <- as.numeric(snakemake@params["n_sample"])
dim_prescreening <- as.numeric(snakemake@params["dim_prescreening"])
if (dim_prescreening < 0) {
	dim_prescreening <- Inf
}

output <- as.character(snakemake@output)

blasso <- train_blasso(sequences,
	wavelength,
	pos_file = positions,
	amino_feature_file = amino_acid_feature,
	base_wl_file = zero,
	n_sample = n_sample,
	dim_prescreening = dim_prescreening
)

save(result, file = output)
