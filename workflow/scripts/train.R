
library(seqinr)
library(pracma)
library(monomvn)

source("workflow/scripts/seqBlasso.R")

input <- snakemake@input

wavelength <- input$wavelength
positions  <- input$positions
sequences  <- input$sequences
zero       <- input$zero
amino_acid_feature <- input$amino_acid_feature
output <- as.character(snakemake@output)

n_sample <- as.numeric(snakemake@params$n_sample)
dim_prescreening <- as.numeric(snakemake@params$dim_prescreening)
if (dim_prescreening < 0) {
	dim_prescreening <- Inf
}

blasso <- train_blasso(sequences,
	wavelength,
	pos_file = positions,
	amino_feature_file = amino_acid_feature,
	base_wl_file = zero,
	n_sample = n_sample,
	dim_prescreening = dim_prescreening
)

save(blasso, file = output)
