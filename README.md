Snakemake workflow for BLASSO-for-Rhodopsins
============================================

The workflow provides a wrapper for the [BLASSO-for-Rhodopsins](http://www-als.ics.nitech.ac.jp/~karasuyama/BLASSO-for-Rhodopsins/) code for prediction of microbial rhodopsin absorption spectra.
All the credit goes to the authors of the original software, see [Exploration of natural red-shifted rhodopsins using a machine learning-based Bayesian experimental design](https://www.nature.com/articles/s42003-021-01878-9). BLASSO-for-Rhodopsins was downloaded from [http://www-als.ics.nitech.ac.jp/~karasuyama/BLASSO-for-Rhodopsins/BLASSO-Rhodopsin.zip](http://www-als.ics.nitech.ac.jp/~karasuyama/BLASSO-for-Rhodopsins/BLASSO-Rhodopsin.zip), archived on 2021-02-13. The original `Data` folder is located in `original/Data`, the `seqBlasso.R` script is located under `workflow/scripts/` and is used verbatim.

By default, the workflow does predictions for target sequences in `targets/*.fasta`. If the training file `training/blasso.RData` does not exist -- it will be re-recreated, which might take some time.

Note that if when using conda you encounter an error about missing `libRlapack.so`, you have to create a symlink: `ln -s liblapack.so .snakemake/conda/5f7bb0321af85efc62bd7af40c446e45/lib/libRlapack.so`.
