# this script reads big table of AUC and RMSD estimates

library(optparse)
library(scales) # for transparency
library(readr)
library(ochoalabtools)
# shared code with another plotting function that spans several datasets
source('plot-auc-rmsd.R')

# constants
# output file name (big table)
file_table <- 'sum.txt'
rep_max <- 50

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--bfile", type = "character", default = NA, 
                help = "Base name for input plink files (.BED/BIM/FAM)", metavar = "character"),
    make_option("--complete", action = "store_true", default = FALSE, 
                help = "Plot only complete replicates (useful for partial runs)"),
    make_option("--fes", action = "store_true", default = FALSE, 
                help = "Use FES instead of RC trait model"),
    make_option("--herit", type = "double", default = 0.8, 
                help = "heritability", metavar = "double"),
    make_option("--m_causal_fac", type = "double", default = 10,
                help = "Proportion of individuals to causal loci", metavar = "double")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name <- opt$bfile
fes <- opt$fes
m_causal_fac <- opt$m_causal_fac
herit <- opt$herit

# stop if name is missing
if ( is.na(name) )
    stop('`--bfile` terminal option is required!')

# final versions for papers
# pure plink version is only PCA version, add LMM too
method_to_label <- list(
    'pca-plink-pure' = 'PCA',
    gcta = 'LMM'
)
name_base <- 'sum-'
# hardcoded same order as method_to_label
method_cols <- c(
    'red',
    'blue'
)

# move to where the data is
setwd( '../data/' )
setwd( name )

# if fes is true, move to directory containing input and outputs
dir_out <- ''
if ( fes )
    dir_out <- paste0( dir_out, 'fes/' )
if ( m_causal_fac != 10 )
    dir_out <- paste0( dir_out, 'm_causal_fac-', m_causal_fac, '/' )
if ( herit != 0.8 )
    dir_out <- paste0( dir_out, 'h', herit, '/' )
# move there
setwd( dir_out )

# read the big table!
tib <- read_tsv(
    file_table,
    col_types = 'ciiddd'
)

# extract methods from table itself
methods <- names( method_to_label ) # not from table, but from hardcoded map, always lists PCA first!
# and PCs too, numerically
pcs <- sort( as.numeric( unique( tib$pc ) ) )

# a hack to only plot complete reps
if ( opt$complete ) {
    for ( rep in 1 : rep_max ) {
        # subset for a bit
        tib2 <- tib[ tib$rep == rep, ]
        rep_good <- TRUE # boolean that remembers if things were ok for this rep
        # navigate main methods
        for ( method in methods ) {
            # now just count cases, must have 91 to be complete!
            if ( sum( tib2$method == method ) != 91 )
                rep_good <- FALSE
        }
        # remove rep if it wasn't good!
        if ( !rep_good )
            tib <- tib[ tib$rep != rep, ]
    }
}

################
### NUM REPS ###
################

# this is for internal use only, to monitor progress
# (when complete, all method/pc cases will have all replicates and this'll be entirely boring)

# want a table that counts, for every method and PC, the number of replicates so far
# direct way (old C style)
# initialize matrix
counts <- matrix(
    NA,
    nrow = length(methods),
    ncol = length(pcs),
    dimnames = list(
        methods = methods,
        pcs = pcs
    )
)
# navigate and fill cases
for ( i in 1 : nrow( counts ) ) {
    method <- methods[ i ]
    # subset big table
    tib_i <- tib[ tib$method == method, ]
    for ( j in 1 : ncol( counts ) ) {
        pc <- pcs[ j ]
        # this is the count we want
        counts[ i, j ] <- sum( tib_i$pc == pc )
    }
}
# don't plot anything if this is trivial
if ( any( counts < rep_max ) ) {
    # simple plot of progress (shows partial run status, not all replicates complete)
    fig_start(
        paste0( name_base, 'num-reps' ),
        width = 12
    )
    barplot(
        counts,
        beside = TRUE,
        xlab = '# PCs',
        ylab = 'replicates',
        legend.text = unlist( method_to_label ),
        col = method_cols
    )
    fig_end()
}

################
### NUM FAIL ###
################

# counts cases that failed to yield any results (should be exclusively GCTA under small sample sizes, where model did not converge)

# want a table that counts, for every method and PC, the number of replicates so far
# direct way (old C style)
# initialize matrix
counts <- matrix(
    0,
    nrow = length(methods),
    ncol = length(pcs),
    dimnames = list(
        methods = methods,
        pcs = pcs
    )
)
# navigate and fill cases
for ( i in 1 : nrow( counts ) ) {
    method <- methods[ i ]
    # subset big table
    tib_i <- tib[ tib$method == method, ]
    for ( j in 1 : ncol( counts ) ) {
        pc <- pcs[ j ]
        # this is the count we want
        counts[ i, j ] <- sum( is.na( tib_i$rmsd[ tib_i$pc == pc ] ) )
    }
}
# don't plot anything if this is trivial
if ( any( counts > 0 ) ) {
    # simple plot of counts
    fig_start(
        paste0( name_base, 'num-fail' ),
        width = 12
    )
    barplot(
        counts,
        beside = TRUE,
        xlab = '# PCs',
        ylab = '# Failures',
        legend.text = unlist( method_to_label ),
        col = method_cols,
        args.legend = list(
            x = 'topleft',
            bty = 'n'
        )
    )
    fig_end()
}

#########################
### AUC/RMSD BOXPLOTS ###
#########################

# default legend position
legend_pos <- 'topright'
# hack for small sample size simulation only
if ( name == 'sim-n100-k10-f0.1-s0.5-g1' )
    legend_pos <- 'bottomleft'

# gather data into lists, best for boxplots
data <- tibble_to_lists_rmsd_auc( tib )

# actually make plot
# for main figure output (goes on papers, presentations)
lineplots_rmsd_auc(
    name_out = paste0( name_base, 'rmsd-auc' ),
    data_rmsd = data$rmsd,
    data_auc = data$auc,
    r_max = max( pcs ),
    legend_pos = legend_pos
)
