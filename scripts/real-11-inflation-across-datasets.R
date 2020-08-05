# an earlier figure plotted inflation factors vs RMSD in a single dataset, revealing a clean curve
# is the curve the same across datasets with different sample sizes, etc?

library(readr)
library(tibble)
library(ochoalabtools)

# FOR FAILED EXPERIMENT AT END OF SCRIPT
# load new functions from external scripts
dir_orig <- getwd()
setwd("../../scripts") # scripts of main GAS project
source('rmsd.R')
setwd( dir_orig ) # go back to where we were

# constants
# output file name (big table)
file_table <- 'sum.txt'
# for main figure output (goes on papers)
name_out <- 'sum-rmsd-auc'

# names of datasets (paired list)
datasets <- tibble(
    name_short = c(
        'large',
        'small',
        'family',
        'HO'
    ),
    name_long = c(
        'sim-n1000-k10-f0.1-s0.5-g1',
        'sim-n100-k10-f0.1-s0.5-g1',
        'sim-n1000-k10-f0.1-s0.5-g20',
        'HoPacAll_ld_prune_1000kb_0.3'
    ),
    col = 1:4
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
range_rmsd <- range( unlist( rmsds ), na.rm = TRUE )
range_lambda <- range( unlist( lambdas ), na.rm = TRUE )

# plot in base data dir
fig_start(
    'sum-rmsd-vs-lambda',
    mar_t = 1,
    mar_r = 0.3
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
# guide lines
abline( v = 0, lty = 2, col = 'gray' )
abline( h = 1, lty = 2, col = 'gray' )
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
    'bottomright',
    datasets$name_short,
    text.col = datasets$col,
    pch = NA,
    bty = 'n'
)
#fig_end()

#############

# test a hypothesis about generating this curve
# failed experiment, curves don't track well for high RMSD, lambda, so meh

# assume the data is coming from a chi-square with different degrees of freedom
m <- 1000 # samples (numbers of p-vals)
df <- 1 # in GWAS this is always 1 (genotype vector is added in H1, removed in H0, all else the same)
ps <- ( 1 : m ) / ( m + 1 ) # percentiles to retrieve, so we barely miss 0 and 1 edges
x_m <- qchisq( 0.5, df = df ) # median statistic

# data we want
n <- 100 # resolution of curve
# deltas in a log scale
deltas <- ( -n : n ) / n * 0.8
# in real scale
deltas <- 10^( deltas )
l <- length( deltas )
lambdas_sim <- vector( 'numeric', l )
rmsds_sim <- vector( 'numeric', l )
for ( i in 1 : l ) {
    message( 'i: ', i )
    delta <- deltas[i]
    # this generates the alternate chi-squared quantiles
    xs <- qchisq( ps, df = delta )
    # then runs them through the cumulative distribution (standard version) to get the (bad) pvalues
    pvals <- pchisq( xs, df = df )
    # lambdas (directly on raw chi-squared stats, not through p-values)
    lambdas_sim[i] <- median( xs ) / x_m
    # compute error metric (RMSD)
    rmsds_sim[i] <- rmsd( pvals, ps )
    # add sign depending on delta
    if ( delta < 1 )
        rmsds_sim[i] <- - rmsds_sim[i]
}

lines(
    rmsds_sim,
    lambdas_sim
)

fig_end()
