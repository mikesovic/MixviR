MixviR can be installed with...
devtools::install_github("mikesovic/MixviR/MixviR_vX.X")

Primary function in MixviR is call_mutations(). It requires the 'samp.dir' argument, which is a directory that stores one or more csv files with the columns...
POS, REF, ALT, DP, ALT_COUNT

Such csv files can be genereated with scripts in 'variant_calling' dir.

Note that as of v1.2, only genomes with a single chromosome are supported.

