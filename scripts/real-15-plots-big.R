# this script reads big table of AUC and RMSD estimates
# version that combines three datasets in a single panel

library(optparse)
library(scales) # for transparency
library(readr)
library(tibble)
library(ochoalabtools)
# shared code with another plotting function that spans several datasets
source('plot-auc-rmsd.R')

#################
### CONSTANTS ###
#################

# full page size
## genetics dimensions
## width max 20cm = 7.87402 in
## height max 25cm = 9.84252 in
height_max <- 9.84252
width_max <- 7.87402
# margin for panels with titles
mar1 <- c(1.5, 3, 2.5, 0) + 0.2
# margin for panels without titles
mar2 <- c(1.5, 3, 0, 0) + 0.2

# input file name (big table), per dataset
file_table <- 'sum.txt'
# hardcoded params
rep_max <- 50
r_max <- 90
pcs <- 0 : r_max # global needed
# pure plink version is only PCA version, add LMM too
method_to_label <- list(
    'pca-plink-pure' = 'Fixed effects (PCA)',
    gcta = 'Mixed effects (LMM)'
)
# extract methods from table itself
methods <- names( method_to_label ) # not from table, but from hardcoded map, always lists PCA first!
# hardcoded same order as method_to_label
method_cols <- c(
    'red',
    'blue'
)
# names of datasets (paired list)
# simulated datasets are default
datasets <- tibble(
    name_short = c(
        'Large sample size sim.',
        'Small sample size sim.',
        'Family structure sim.'
    ),
    name_long = c(
        'sim-n1000-k10-f0.1-s0.5-g1',
        'sim-n100-k10-f0.1-s0.5-g1',
        'sim-n1000-k10-f0.1-s0.5-g20'
    )
)
datasets_real <- tibble(
    name_short = c(
        'Human Origins',
        'HGDP',
        '1000 Genomes'
    ),
    name_long = c(
        'HoPacAll_ld_prune_1000kb_0.3',
        'hgdp_wgs_autosomes_ld_prune_1000kb_0.3',
        'all_phase3_filt-minimal_ld_prune_1000kb_0.3_thinned-0.1'
    )
)
# output names
name_out <- 'rmsd-auc-sim'
name_out_real <- 'rmsd-auc-real'

# default legend position
legend_pos <- 'topright'
## # hack for small sample size simulation only
## if (name == 'sim-n100-k10-f0.1-s0.5-g1')
##     legend_pos <- 'bottomleft'


############
### ARGV ###
############

# define options
option_list = list(
    make_option("--real", action = "store_true", default = FALSE, 
                help = "Process real datasets (default: simulated datasets)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# switch to process real datasets if requested
if ( opt$real ) {
    datasets <- datasets_real
    name_out <- name_out_real
}

################
### DATA/FIG ###
################

# move to where the data is
setwd( '../data/' )

# start PDF
fig_start(name_out, width = width_max, height = height_max)
# add lower margin, so inner margins can be even smaller
par( oma = c(1.5, 0, 0, 0) )
# 6 panels/rows
# since some have titles and others don't, looks better to have this uneven sizing
layout(
    cbind( 1:6 ),
    heights = rep( c(1.24, 1), 3 )
)
# change tick mark frequency on x axis
par( lab = c(10, 3, 7) )

# load each of our datasets
for ( index_dataset in 1 : nrow( datasets ) ) {
    name <- datasets$name_long[ index_dataset ]
    setwd( name )

    # read the big table!
    tib <- read_tsv(
        file_table,
        col_types = 'ciiddd'
    )

    # gather data into lists, best for boxplots
    # store in bigger structure
    data <- tibble_to_lists_rmsd_auc( tib )
    
    # top panel: RMSD
    # set margin for titles
    par( mar = mar1 )
    # main plot
    lineplots_rmsd_auc_one_panel( data$rmsd, lab_rmsd, r_max, main = datasets$name_short[ index_dataset ] )
    # add panel letter
    panel_letter( toupper( letters[ index_dataset ] ) )
    
    # add legend to first panel only
    if ( index_dataset == 1 ) {
        legend(
            legend_pos,
            unlist( method_to_label ),
            text.col = method_cols,
            bty = 'n'
        )
        # add second legend explaining quartiles, etc
        # NOTE: only small sample size sim has diff legend_pos == 'bottomleft'
        legend_pos_quants <- 'right' # if ( legend_pos == 'topright' ) 'bottomright' else 'topright'
        leg_mean_quarts( alpha_q, alpha_e, x = legend_pos_quants, cex = 1 )
    }

    # bottom panel: AUC
    # reduce top margin for next panel
    par( mar = mar2 )
    # main plot
    lineplots_rmsd_auc_one_panel( data$auc, lab_auc, r_max, guide_max = TRUE )

    # go back down, for next dataset
    setwd( '..' )
}

# add outer margin label
mtext(
    "Number of PCs (r)",
    side = 1,
    line = 0.5,
    adj = 0.55,
    outer = TRUE
)
fig_end()