# an earlier figure plotted inflation factors vs RMSD in a single dataset, revealing a clean curve
# is the curve the same across datasets with different sample sizes, etc?

library(readr)
library(tibble)
library(ochoalabtools)

# constants
# output file name (big table)
file_table <- 'sum.txt.gz'
# for main figure output (goes on papers)
name_out <- 'sum-rmsd-auc'

# names of datasets (paired list)
datasets <- tibble(
    name_short = c(
        'large',
        ## 'small',
        ## 'family',
        'HO'
    ),
    name_long = c(
        'sim-n1000-k10-f0.1-s0.5-g1',
        ## 'sim-n100-k10-f0.1-s0.5-g1',
        ## 'sim-n1000-k10-f0.1-s0.5-g20',
        'HoPacAll_ld_prune_1000kb_0.3'
    ),
    col = 1:2
)

# data we want to collect
# separate both so axis ranges are easier to extract
rmsds <- list()
lambdas <- list()

# move to where the data is
setwd( '../data/' )

# load each dataset
for ( i in 1 : nrow( datasets ) ) {
    # enter dir
    setwd( datasets$name_long[ i ] )
    
    # read the big table!
    tib <- read_tsv(
        file_table,
        col_types = 'ciiddd'
    )
    
    # save only the two columns of interest
    rmsds[[ i ]] <- tib$rmsd
    lambdas[[ i ]] <- tib$lambda
    
    # go back down
    setwd( '..' )
}

######################
### LAMBDA vs RMSD ###
######################

# plotting labels
lab_rmsd <- expression( bold( SRMSD[p] ) )
lab_lambda <- expression( bold( paste("Inflation Factor (", lambda, ")") ) )

# find common data range
range_rmsd <- range( unlist( rmsds ) )
range_lambda <- range( unlist( lambdas ) )

# plot in base data dir
fig_start(
    'sum-rmsd-vs-lambda'
)
# start base plot
plot(
    NA,
    xlim = range_rmsd,
    ylim = range_lambda,
    xlab = lab_rmsd,
    ylab = lab_lambda,
    log = 'y'
)
# add per-dataset points
for ( i in 1 : nrow( datasets ) ) {
    points(
        rmsds[[ i ]],
        lambdas[[ i ]],
        pch = '.',
        col = datasets$col[ i ]
    )
}
# legend
legend(
    'topleft',
    datasets$name_short,
    text.col = datasets$col,
    pch = NA,
    bty = 'n'
)
fig_end()
