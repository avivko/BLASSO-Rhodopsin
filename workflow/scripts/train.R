
library(seqinr)
library(pracma)
library(monomvn)

snakemake@source("seqBlasso.R")

with(snakemake@input, {
    wavelength <<- wavelength
    positions  <<- positions
    sequences  <<- sequences
    zero       <<- zero
    amino_acid_feature <<- amino_acid_feature
})
with(snakemake@params, {
    n_sample    <<- as.numeric(n_sample)
    dim_prescreening <<- as.numeric(dim_prescreening)
})
output <- as.character(snakemake@output)

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
