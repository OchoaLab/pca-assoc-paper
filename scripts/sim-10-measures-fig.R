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
model_pca <- 'pca-plink-pure'
model_lmm <- 'gcta'
# dir for archived data, for sims it's present on viiiaX6 only
dir_archive <- '~/dbs2/PCA/'

# define curves to look at
# good case (LMM with zero PCs)
# bad case deflated (LMM with 90 PCs)
# bad case inflated (PCA with zero PCs)
models <- c(
    model_lmm,
    model_lmm,
    model_pca
)
model_names <- paste0( 'M', 1:length(models) )
rs <- c( 0, 90, 0 )
colors <- c(
    'green',
    'blue',
    'red'
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

read_pvals <- function( model, n_pcs, m_loci = NA ) {
    # file to read
    file_pvals <- paste0( 'pvals_', model, '_', n_pcs, '.txt.gz' )
    
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
pvals <- vector( 'list', length(models) )
for ( i in 1 : length(models) ) {
    pvals[[i]] <- read_pvals( models[i], rs[i] )
}

# switch back to ordinary data location
setwd( dir_orig )

panel_letter_adj <- -0.5

# start figure
width <- fig_width() / 2
fig_start(
    'measures-illustration',
    width = width,
    height = width * 1.18,
    mar_t = 2, # for panel letters
    mar_l = 3.5, # more space in this case, everything is squished
    mar_b = 3.5 
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
    xlab = 'Expected null p-val.', # quantiles
    ylab = 'Observed null p-val.'
)
# plot reference line in background
abline( 0, 1, col = 'gray', lty = 2 )
# add a nice panel letter
panel_letter( 'A', adj = panel_letter_adj )
# add a legend
legend(
    'topleft',
    title = 'Assoc. Models',
    legend = model_names,
    lty = 1,
    col = colors,
    cex = 0.7,
    bty = 'n'
)

srmsd <- vector( 'numeric', length(models) )
# plot backwards so green is on top
for ( i in rev( 1 : length(models) ) ) {
    # calculate SRMSD_p
    obj <- pval_srmsd(pvals[[i]], causal_indexes, detailed = TRUE)
    # record value for second panel
    srmsd[i] <- obj$srmsd
    # plot curve
    lines(
        obj$pvals_unif,
        obj$pvals_null,
        col = colors[i]
    )
    
    # calculate inflation factor, just for reference
    # just print out, will not be in figure
    lambda <- pval_infl( pvals[[i]] )
    message( 'Inflation factor: ', models[i], ': ', lambda )
}

# second panel is just a bar plot
barplot(
    srmsd,
    names.arg = model_names,
    col = colors,
    xlab = 'Assoc. Models',
    ylab = expression( bold( SRMSD[p] ) ),
    ylim = c( min( srmsd ), round( max( srmsd ), 1 ) ) # make axis nicer this way
)
# add a nice panel letter
panel_letter( 'B', adj = panel_letter_adj )

# third panel: AUC plot
# in this case we don't start it until inside the loop (for first element only)
aucpr <- vector( 'numeric', length(models) )
# plot backwards so green is on top
for ( i in rev( 1 : length(models) ) ) {
    # calculate AUC_PR
    obj <- pval_aucpr(pvals[[i]], causal_indexes, curve = TRUE)
    # record value for next panel
    aucpr[i] <- obj$auc.integral
    # add for all but last model (first in reverse order)
    add <- i < length( models )
    # plot curve
    plot(
        obj,
        color = colors[i],
        auc.main = FALSE,
        main = '',
        add = add
    )
}
panel_letter( 'C', adj = panel_letter_adj )

# last panel is just a bar plot
barplot(
    aucpr,
    names.arg = model_names,
    col = colors,
    xlab = 'Assoc. Models',
    ylab = expression( bold( AUC[PR] ) )
)
# add a nice panel letter
panel_letter( 'D', adj = panel_letter_adj )

fig_end()
