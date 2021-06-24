# writes a basic table other scripts that scan several datasets load
# provides mapping between paper names and internal dir codes

library(readr)
library(tibble)

#################
### CONSTANTS ###
#################

# names of datasets (paired list)
datasets <- tibble(
    name_paper = c(
        'Large sample size sim.',
        'Small sample size sim.',
        'Family structure sim.',
        'Human Origins',
        'HGDP',
        '1000 Genomes',
        'Human Origins sim.',
        'HGDP sim.',
        '1000 Genomes sim.'
    ),
    name_dir = c(
        'sim-n1000-k10-f0.1-s0.5-g1',
        'sim-n100-k10-f0.1-s0.5-g1',
        'sim-n1000-k10-f0.1-s0.5-g20',
        'HoPacAll_ld_prune_1000kb_0.3_maf-0.01',
        'hgdp_wgs_autosomes_ld_prune_1000kb_0.3_maf-0.01',
        'all_phase3_filt-minimal_ld_prune_1000kb_0.3_maf-0.01',
        'HoPacAll_ld_prune_1000kb_0.3_maf-0.01_sim',
        'hgdp_wgs_autosomes_ld_prune_1000kb_0.3_maf-0.01_sim',
        'all_phase3_filt-minimal_ld_prune_1000kb_0.3_maf-0.01_sim'
    ),
    type = c(
        'Admix.',
        'Admix.',
        'Admix.+Pedig.',
        'Real',
        'Real',
        'Real',
        'Tree',
        'Tree',
        'Tree'
    ),
    # colors for several plots
    col = c(1, 2, 3, 4, 5, 6, 4, 5, 6),
    # line types, for MAF plots only
    lty = c(2, 2, 2, 1, 1, 1, 2, 2, 2)
)

# move to where the data is
setwd( '../data/' )

# write table to a file
write_tsv( datasets, file = 'datasets.txt' )
