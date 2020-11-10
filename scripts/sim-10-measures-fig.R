# this script creates an illustration, using the simulated data, of the SRMSD and AUCPR measures starting from more traditional, easy-to-understand p-p and PR plots
# assumes p-values are "archived"

library(readr)
library(simtrait) # for pval_srmsd, pval_aucpr, pval_infl
library(ochoalabtools) # for plots

# the replicate to use makes figure vary a lot, need to choose carefully
# set via terminal to tweak more quickly
args <- args_cli()
rep <- args[1]

# constants
# use this simulation only for illustration
name <- 'sim-n1000-k10-f0.1-s0.5-g1'
# these have weird names, simplify mapping
method_pca <- 'pca-plink-pure'
method_lmm <- 'gcta'
# dir for archived data, for sims it's present on viiiaX6 only
dir_archive <- '~/dbs2/PCA/'

# define curves to look at
# good case (LMM with zero PCs)
# bad case inflated (PCA with zero PCs)
# bad case deflated (LMM with 90 PCs)
methods <- c(
    method_lmm,
    method_pca,
    method_lmm
)
method_names <- 1:length(methods)
rs <- c( 0, 0, 90 )
colors <- c(
    'green',
    'red',
    'blue'
)

# load causal_indexes (in ordinary data location)
setwd( '../data/' )
setwd( name )
# remember where we were
dir_orig <- getwd()
# actually go into rep
dir_rep <- paste0( 'rep-', rep )
setwd( dir_rep )
load( 'simtrait.RData' )


# now move to where the archived p-value data is
setwd( dir_archive )
setwd( name )
setwd( dir_rep )

read_pvals <- function( method, n_pcs, m_loci = NA ) {
    # file to read
    file_pvals <- paste0( 'pvals_', method, '_', n_pcs, '.txt.gz' )
    
    # if there's no input file, die!
    if ( !file.exists( file_pvals ) )
        stop( 'P-value file missing: ', file_pvals )

    # read the file
    pvals <- read_lines(
        file_pvals,
        na = 'NA' # use this to prevent warnings
    )

    # Null cases are not acceptable for this plot
    if ( length(pvals) == 1 && pvals == 'NULL' )
        stop( 'Got NULL p-value file!' )

    # make sure length is correct! (optional)
    if ( !is.na( m_loci ) && length( pvals ) != m_loci )
        stop( 'File has ', length( pvals ), ' p-values, expected ', m_loci )
    
    # convert strings to numeric
    pvals <- as.numeric( pvals )
    
    # return those p-values
    return( pvals )
}

# read p-values in loop, shared by all panels
pvals <- vector( 'list', length(methods) )
for ( i in 1 : length(methods) ) {
    pvals[[i]] <- read_pvals( methods[i], rs[i] )
}

# switch back to ordinary data location
setwd( dir_orig )

# start figure
fig_start(
    'measures-illustration',
    width = 4,
    height = 4.6,
    mar_t = 2 # for panel letters
)

# create panels
par( mfrow = c(2, 2) )

# for first two panels

# start plotting the first panel right away
plot(
    NA,
    type = 'n',
    xlim = c(0, 1),
    ylim = c(0, 1),
    xlab = 'Expected null p-value', # quantiles
    ylab = 'Observed null p-value'
)
# plot reference line in background
abline( 0, 1, col = 'gray', lty = 2 )
# add a nice panel letter
panel_letter('A')
# add a legend
legend(
    'bottomright',
    title = 'Methods',
    legend = method_names,
    lty = 1,
    col = colors,
    cex = 0.7
)

srmsd <- vector( 'numeric', length(methods) )
for ( i in 1 : length(methods) ) {
    # calculate SRMSD_p
    obj <- pval_srmsd(pvals[[i]], causal_indexes, detailed = TRUE)
    # record value for second panel
    srmsd[i] <- obj$srmsd
    # plot curve
    lines(
        obj$pvals_null,
        obj$pvals_unif,
        col = colors[i]
    )
    
    # calculate inflation factor, just for reference
    # just print out, will not be in figure
    lambda <- pval_infl( pvals[[i]] )
    message( 'Inflation factor: ', methods[i], ': ', lambda )
}

# second panel is just a bar plot
barplot(
    srmsd,
    names.arg = method_names,
    col = colors,
    xlab = 'Methods',
    ylab = expression( bold( SRMSD[p] ) )
)
# add a nice panel letter
panel_letter('B')

# third panel: AUC plot
# in this case we don't start it until inside the loop (for first element only)
aucpr <- vector( 'numeric', length(methods) )
for ( i in 1 : length(methods) ) {
    # calculate AUC_PR
    obj <- pval_aucpr(pvals[[i]], causal_indexes, curve = TRUE)
    # record value for next panel
    aucpr[i] <- obj$auc.integral
    # plot curve
    plot(
        obj,
        color = colors[i],
        auc.main = FALSE,
        main = '',
        add = (i > 1)
    )
    # add panel letter first time only
    if ( i == 1 )
        panel_letter('C')
}

# last panel is just a bar plot
barplot(
    aucpr,
    names.arg = method_names,
    col = colors,
    xlab = 'Methods',
    ylab = expression( bold( AUC[PR] ) )
)
# add a nice panel letter
panel_letter('D')

fig_end()
