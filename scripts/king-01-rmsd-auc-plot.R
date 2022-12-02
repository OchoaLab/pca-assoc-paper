# this script reads big table of AUC and RMSD estimates
# version that combines three datasets in a single panel

library(optparse)
library(scales) # for transparency
library(readr)
library(tibble)
library(ochoalabtools)
# shared code with another plotting function that spans several datasets
# here merely for lab_rmsd and lab_auc
source('plot-auc-rmsd.R')

#################
### CONSTANTS ###
#################

# input file name (big table), per dataset
file_table <- 'sum.txt'
# hardcoded params
r_pca <- 20
r_lmm <- 0
fes <- TRUE
# P-value threshold for Wilcoxon 2-sample paired 1-tail comparisons
p_cut <- 0.01
# pure plink version is only PCA version, add LMM too
method_to_label <- list(
    'pca-plink-pure' = 'PCA',
    gcta = 'LMM'
)
# extract methods from table itself
methods <- names( method_to_label ) # not from table, but from hardcoded map, always lists PCA first!
n_methods <- length( methods )
# hardcoded same order as method_to_label
method_cols <- c(
    'red',
    'blue'
)
# names of datasets
# base is real datasets (for getting "paper" names), but we're really looking at king-cutoff versions
names_dir <- c(
    'HoPacAll_ld_prune_1000kb_0.3_maf-0.01',
    'hgdp_wgs_autosomes_ld_prune_1000kb_0.3_maf-0.01_geno-0.1',
    'tgp-nygc-autosomes_ld_prune_1000kb_0.3_maf-0.01'
)
# suffix to actual dirs
suffix <- '_king-cutoff-4'
# output name
name_out <- paste0( 'rmsd-auc', suffix )

############
### ARGV ###
############

# define options
option_list = list(
    make_option(c("-r", "--rep"), type = "integer", default = 50,
                help = "Max replicates", metavar = "int"),
    make_option("--herit", type = "double", default = 0.8, 
                help = "heritability", metavar = "double"),
    make_option("--m_causal_fac", type = "double", default = 10,
                help = "Proportion of individuals to causal loci", metavar = "double"),
    make_option("--env1", type = "double", default = NA,
                help = "Variance of 1st (coarsest) level of environment (non-genetic) effects (default NA is no env)", metavar = "double"),
    make_option("--env2", type = "double", default = NA,
                help = "Variance of 2nd (finest) level of environment (non-genetic) effects (default NA is no env)", metavar = "double")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
rep_max <- opt$rep
m_causal_fac <- opt$m_causal_fac
herit <- opt$herit
env1 <- opt$env1
env2 <- opt$env2

# specify location of files to process, as many levels as needed
dir_out <- ''
if ( fes )
    dir_out <- paste0( dir_out, 'fes/' )
if ( m_causal_fac != 10 )
    dir_out <- paste0( dir_out, 'm_causal_fac-', m_causal_fac, '/' )
if ( herit != 0.8 )
    dir_out <- paste0( dir_out, 'h', herit, '/' )
if ( !is.na( env1 ) )
    dir_out <- paste0( dir_out, 'env', env1, '-', env2, '/' )


# hack version, originally `tibble_to_lists_rmsd_auc` from plot-auc-rmsd.R but simplified to single PC per method, plus check it's the expected value
tibble_to_lists_rmsd_auc_king <- function( tib ) {
    # initialize lists
    data_rmsd <- list()
    data_auc <- list()
    for ( method in methods ) {
        # get label for figure
        method_label <- method_to_label[[ method ]]
        # subset big table
        tib_i <- tib[ tib$method == method, ]
        # check that all data is as desired
        # assumes first method is PCA, second is LMM (assumption made elsewhere too)
        if ( method == methods[1] )
            stopifnot( all( tib_i$pc == r_pca ) )
        if ( method == methods[2] )
            stopifnot( all( tib_i$pc == r_lmm ) )
        # copy sublists to big lists
        data_rmsd[[ method_label ]] <- tib_i$rmsd
        data_auc[[ method_label ]] <- tib_i$auc
    }
    # return data in yet another list
    return(
        list(
            rmsd = data_rmsd,
            auc = data_auc
        )
    )
}

# also adapted for the slightly different structure here, from `report_cross_method` in real-13-stats.R
report_cross_method_king <- function( method_to_vals, metric ) {
    # here we haven't taken |srmsd| before function, so do that here before everything else
    if ( metric == 'rmsd' )
        method_to_vals <- lapply( method_to_vals, abs )
    
    # decide which performed best, by mean
    method_to_mean <- sapply( method_to_vals, mean )
    if ( metric == 'auc' ) {
        method_best <- method_to_label[[ which.max( method_to_mean ) ]]
        alternative <- 'l' # alternative is lesser than max
    } else if ( metric == 'rmsd' ) {
        method_best <- method_to_label[[ which.min( method_to_mean ) ]]
        alternative <- 'g' # alternative is greater than min
    }
    # to get direction right, need to identify worst method too
    method_worst <- setdiff( unlist( method_to_label ), method_best )
    
    # and decide if they're statistically significantly different or not
    sig <- wilcox.test(
        method_to_vals[[ method_worst ]],
        method_to_vals[[ method_best ]],
        paired = TRUE,
        alternative = alternative
    )$p.value < p_cut

    if ( !sig )
        method_best <- 'tie'

    return( method_best )
}



#################
### LOAD DATA ###
#################

# move to where the data is
setwd( '../data/' )

# remember this location, using global path, for easily returning across multiple levels when done
dir_base <- getwd()

# get panel names from here
datasets <- read_tsv( 'datasets.txt', col_types = 'cccii' )
# subset and reorder as desired
indexes <- match( names_dir, datasets$name_dir )
datasets <- datasets[ indexes, ]

# load each of our datasets
# all data goes in `data`
data <- list()
for ( index_dataset in 1 : nrow( datasets ) ) {
    name <- datasets$name_dir[ index_dataset ]
    # add suffix to actually read king-cutoff filtered version's data
    name <- paste0( name, suffix )
    setwd( name )
    # move in one more level in this case
    if ( dir_out != '' )
        setwd( dir_out )

    # read the big table!
    tib <- read_tsv(
        file_table,
        col_types = 'ciiddd'
    )

    # gather data into lists, best for boxplots
    # store in bigger structure, by index here!
    data_i <- tibble_to_lists_rmsd_auc_king( tib )
    
    # do the stats!
    message( name )
    message( 'RMSD: ', report_cross_method_king( data_i$rmsd, 'rmsd' ) )
    message( 'AUC: ', report_cross_method_king( data_i$auc, 'auc' ) )

    data[[ index_dataset ]] <- data_i
    
    # return to base when done with this dataset
    setwd( dir_base )
}

############
### PLOT ###
############

# calculate shared ylim's
ylim_rmsd <- range( -srmsd_cut, srmsd_cut, sapply( data, function (x) range( unlist( x$rmsd ) ) ) ) # include tolerance band
ylim_auc <- range( 0, sapply( data, function (x) range( unlist( x$auc ) ) ) ) # include zero

# new level to this hierarchy
if ( dir_out != '' ) {
    # create if it didn't already exist
    if( !dir.exists( dir_out ) )
        dir.create( dir_out, recursive = TRUE )
    # now move in there
    setwd( dir_out )
}

# start PDF
width <- fig_width( )
fig_start(
    name_out,
    width = width,
    height = width * 2 / 3
)
# margin for panels with titles
mar1 <- c(1.5, 1.5, 2.5, 0) + 0.2
# margin for panels without titles
mar2 <- c(1.5, 1.5, 0, 0) + 0.2
# add lower margin, so inner margins can be even smaller
# same with left (only outer needs a big margin)
par( oma = c(1.5, 1.5, 0, 0) )
# draw by columns (to have AUCs on top, RMSDs on bottom, each column is a different dataset)
# since some have titles and others don't, looks better to have this uneven sizing
layout(
    matrix( 1:6, nrow = 2, ncol = 3 ),
    heights = c(1.15, 1)
)

# process each of our datasets
for ( index_dataset in 1 : nrow( datasets ) ) {
    data_i <- data[[ index_dataset ]]

    # top panel: RMSD
    # set margin for titles
    par( mar = mar1 )
    # start blank plot
    plot(
        NA,
        main = datasets$name_paper[ index_dataset ],
        xlab = '',
        ylab = '',
        xlim = c(0.5, n_methods + 0.5), # this is boxplot's default
        ylim = ylim_rmsd,
        axes = FALSE # boxplot(add=TRUE) overplots these, let it
    )
    # add y-label on leftmost panels only, make sure it can go on outer margin space
    if ( index_dataset == 1 )
        mtext( lab_rmsd, side = 2, line = 1.5, xpd = NA )
    # mark zero line
    abline( h = 0, lty = 2, col = 'gray' )
    # mark band where inflation is practically low enough
    polygon(
        c(0, n_methods + 1, n_methods + 1, 0),
        c( srmsd_cut, srmsd_cut, -srmsd_cut, -srmsd_cut ),
        col = alpha( 'gray', alpha_q ), # the darker of the two areas to see it more easily
        border = FALSE # no border lines
    )
    # main plot add after bg gray bands
    boxplot(
        data_i$rmsd,
        range = 0, # wiskers extend to outliers
        border = method_cols,
        col = alpha( method_cols, alpha_q ), # inside box, to match big plots in transparency
        add = TRUE
    )
    # add panel letter
    panel_letter( toupper( letters[ index_dataset ] ) )

    # bottom panel: AUC
    # reduce top margin for next panel
    par( mar = mar2 )
    # main plot
    boxplot(
        data_i$auc,
        xlab = '',
        ylab = '',
        ylim = ylim_auc, # c(0, max( unlist( data_i$auc ) ) ), # include zero
        range = 0, # wiskers extend to outliers
        border = method_cols,
        col = alpha( method_cols, alpha_q ) # inside box, to match big plots in transparency
    )
    # mark zero line
    abline( h = 0, lty = 2, col = 'gray' )
    # add y-label on leftmost panels only, make sure it can go on outer margin space
    if ( index_dataset == 1 )
        mtext( lab_auc, side = 2, line = 1.5, xpd = NA )
}

# add outer margin label
mtext(
    'Assoc. Model',
    side = 1,
    line = 0.5,
    adj = 0.55,
    outer = TRUE
)
fig_end()
