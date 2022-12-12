# this script reads big table of AUC and RMSD estimates
# version that combines three datasets in a single panel

# TODO:
# - switch to `ochoalabtools` solution for fig dimensions

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
width_max <- fig_width()
height_max <- width_max # square?
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
    'pca-plink-pure' = 'PCA',
    gcta = 'LMM'
)
# hardcoded same order as method_to_label
method_cols <- c(
    'red',
    'blue'
)
# names of datasets
# simulated datasets are default
names_dir <- c(
    'sim-n1000-k10-f0.1-s0.5-g1',
    'sim-n100-k10-f0.1-s0.5-g1',
    'sim-n1000-k10-f0.1-s0.5-g20'
)
names_dir_real <- c(
    'HoPacAll_ld_prune_1000kb_0.3_maf-0.01',
    'hgdp_wgs_autosomes_ld_prune_1000kb_0.3_maf-0.01_geno-0.1',
    'tgp-nygc-autosomes_ld_prune_1000kb_0.3_maf-0.01'
)
names_dir_real_sim <- paste0( names_dir_real, '_sim' )
# output names
name_out <- 'rmsd-auc-sim'
name_out_real <- 'rmsd-auc-real'
name_out_real_sim <- 'rmsd-auc-real-sim'

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--real", action = "store_true", default = FALSE, 
                help = "Process real datasets (default: simulated datasets)"),
    make_option("--real_sim", action = "store_true", default = FALSE, 
                help = "Process real-sim datasets (default: simulated datasets)"),
    make_option("--fes", action = "store_true", default = FALSE, 
                help = "Use FES instead of RC trait model"),
    make_option("--herit", type = "double", default = 0.8, 
                help = "heritability", metavar = "double"),
    make_option("--m_causal_fac", type = "double", default = 10,
                help = "Proportion of individuals to causal loci", metavar = "double"),
    make_option("--env1", type = "double", default = NA,
                help = "Variance of 1st (coarsest) level of environment (non-genetic) effects (default NA is no env)", metavar = "double"),
    make_option("--env2", type = "double", default = NA,
                help = "Variance of 2nd (finest) level of environment (non-genetic) effects (default NA is no env)", metavar = "double"),
    make_option(c('-l', "--labs"), action = "store_true", default = FALSE, 
                help = "Include LMM with labels data")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# switch to process real datasets if requested
if ( opt$real ) {
    # can't have both of these be on at the same time
    if ( opt$real_sim )
        stop( 'Cannot use "--real" and "--real_sim" at the same time!' )
    # proceed now
    names_dir <- names_dir_real
    name_out <- name_out_real
} else if ( opt$real_sim ) {
    names_dir <- names_dir_real_sim
    name_out <- name_out_real_sim
}
# get values
fes <- opt$fes
m_causal_fac <- opt$m_causal_fac
herit <- opt$herit
env1 <- opt$env1
env2 <- opt$env2
labs <- opt$labs

if ( labs ) {
    # add a third method
    method_to_label <- c( method_to_label, list( 'gcta-labs' = 'LMM lab.' ) )
    method_cols <- c( method_cols, 'green' )
}

# extract methods from table itself
methods <- names( method_to_label ) # not from table, but from hardcoded map, always lists PCA first!


############
### DATA ###
############

# move to where the data is
setwd( '../data/' )

# get panel names from here
datasets <- read_tsv( 'datasets.txt', col_types = 'cccii' )
# subset and reorder as desired
indexes <- match( names_dir, datasets$name_dir )
datasets <- datasets[ indexes, ]

# to load real datasets and get back down easily
dir_orig <- getwd()

# store data here
data <- list()

# if fes is true, move to directory containing input and outputs
dir_out <- ''
if ( fes )
    dir_out <- paste0( dir_out, 'fes/' )
if ( m_causal_fac != 10 )
    dir_out <- paste0( dir_out, 'm_causal_fac-', m_causal_fac, '/' )
if ( herit != 0.8 )
    dir_out <- paste0( dir_out, 'h', herit, '/' )
if ( !is.na( env1 ) )
    dir_out <- paste0( dir_out, 'env', env1, '-', env2, '/' )

for ( name in datasets$name_dir ) {
    setwd( name )
    # move there in this case
    if ( dir_out != '' )
        setwd( dir_out )

    # read the big table!
    tib <- read_tsv(
        file_table,
        col_types = 'ciiddd'
    )

    # gather data into lists, best for boxplots
    # store in bigger structure
    data[[ name ]] <- tibble_to_lists_rmsd_auc( tib )
    
    # go back down, for next dataset
    setwd( dir_orig )
}

###########
### FIG ###
###########

# move there in this case
if ( dir_out != '' ) {
    # create directory if it didn't already exist
    if ( !dir.exists( dir_out ) )
        dir.create( dir_out, recursive = TRUE )
    # now move in
    setwd( dir_out )
}

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

# plot each of our datasets
for ( index_dataset in 1 : nrow( datasets ) ) {
    # get preloaded and preprocessed data
    data_i <- data[[ datasets$name_dir[ index_dataset ] ]]
    
    # top panel: RMSD
    # set margin for titles
    par( mar = mar1 )
    # main plot
    lineplots_rmsd_auc_one_panel( data_i$rmsd, lab_rmsd, r_max, main = datasets$name_paper[ index_dataset ] )
    # add panel letter
    panel_letter( toupper( letters[ index_dataset ] ), adj = -0.05 )

    # add both legends to first panel only
    if ( index_dataset == 1 ) {
        leg_methods( )
        leg_mean_quarts( alpha_q, alpha_e )
    }

    # bottom panel: AUC
    # reduce top margin for next panel
    par( mar = mar2 )
    # main plot
    lineplots_rmsd_auc_one_panel( data_i$auc, lab_auc, r_max, guide_max = TRUE )
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
