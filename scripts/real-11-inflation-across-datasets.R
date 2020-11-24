# an earlier figure plotted inflation factors vs RMSD in a single dataset, revealing a clean curve
# is the curve the same across datasets with different sample sizes, etc?

library(readr)
library(dplyr) # for bind_rows
library(tibble)
library(ochoalabtools)

## # FOR FAILED EXPERIMENT AT END OF SCRIPT
## library(simtrait) # for rmsd

# constants
# output file name (big table)
file_table <- 'sum.txt'
# for main figure output (goes on papers)
name_out <- 'sum-rmsd-auc'
# methods to keep in analysis
methods <- c('pca-plink-pure', 'gcta')
# report on the RMSD predicted for lambda = 1.05
lambda_cut <- 1.05

# names of datasets (paired list)
datasets <- tibble(
    name_short = c(
        'Large sample size sim.',
        'Small sample size sim.',
        'Family structure sim.',
        'Human Origins',
        'HGDP',
        '1000 Genomes'
    ),
    name_long = c(
        'sim-n1000-k10-f0.1-s0.5-g1',
        'sim-n100-k10-f0.1-s0.5-g1',
        'sim-n1000-k10-f0.1-s0.5-g20',
        'HoPacAll_ld_prune_1000kb_0.3',
        'hgdp_wgs_autosomes_ld_prune_1000kb_0.3',
        'all_phase3_filt-minimal_ld_prune_1000kb_0.3_thinned-0.1'
    ),
    col = 1:6
)
# color for fit curve
col_fit_sigmoid <- 'gray40'
lty_fit_sigmoid <- 5
col_fit_loglin <- 'gray'
lty_fit_loglin <- 2
# and guide lines
col_guides <- 'gray95'
lty_guides <- 1

# move to where the data is
setwd( '../data/' )

# big table of interest
# initialize this way, it'll grow correctly
tib_main <- NULL

# load each dataset
for ( i in 1 : nrow( datasets ) ) {
    # enter dir
    setwd( datasets$name_long[ i ] )
    
    # read the big table!
    tib <- read_tsv(
        file_table,
        col_types = 'ciiddd'
    )

    # subset to use only the two methods we talk about in the paper
    tib <- tib[ tib$method %in% methods, ]

    # recall the dataset of origin
    tib$dataset <- datasets$name_short[ i ]

    # boring idea, obviously GCTA was always negative and PCA was always positive SRMSD
##     # color by method
##     tib$col <- 'blue' # default
## #    tib$col[ tib$method == methods[1] ] <- 'blue'
##     tib$col[ tib$method == methods[2] ] <- 'red' # GCTA is red

    # color by dataset
    tib$col <- datasets$col[ i ]
    
    # concatenate into bigger table
    tib_main <- bind_rows( tib_main, tib )
    
    # go back down
    setwd( '..' )
}

######################
### LAMBDA vs RMSD ###
######################

# plotting labels
lab_rmsd <- expression( bold( SRMSD[p] ) )
lab_lambda <- expression( bold( paste("Inflation Factor (", lambda, ")") ) )

# sigmoid fit
# remove NAs for simplicity
# test on one of these, it should be the same for all other metrics
indexes_missing <- is.na( tib_main$rmsd )
stopifnot( all( is.na( tib_main$lambda[ indexes_missing ] ) ) )
stopifnot( all( is.na( tib_main$auc[ indexes_missing ] ) ) )
# subset now
tib_main <- tib_main[ !indexes_missing, ]

# NOTE: x and y are reversed from actual plot
y <- tib_main$rmsd
#x <- log(tib_main$lambda)
x <- tib_main$lambda
x2 <- sort(x)

# fit full data (lambda < 1 seems to mess things up though)
## obj <- nls( y ~ a * (x^b - 1) / (x^b + 1), start = list(a = 1/2, b = 1) )
## Nonlinear regression model
##   model: y ~ a * (x^b - 1)/(x^b + 1)
##    data: parent.frame()
##      a      b 
## 0.5392 0.6046 
##  residual sum-of-squares: 0.3812
## Number of iterations to convergence: 5 
## Achieved convergence tolerance: 2.558e-06

# fit on top data only
indexes <- x > 1
xp <- x[ indexes ]
yp <- y[ indexes ]
objp <- nls( yp ~ a * (xp^b - 1) / (xp^b + 1), start = list(a = 1/2, b = 1) )
#print( objp )
print( coef( objp ) )
## Nonlinear regression model
##   model: yp ~ a * (xp^b - 1)/(xp^b + 1)
##    data: parent.frame()
##      a      b 
## 0.5481 0.6382 
##  residual sum-of-squares: 0.01682
## Number of iterations to convergence: 4 
## Achieved convergence tolerance: 3.601e-07

# report on the RMSD predicted for lambda = 1.05
# lambda is in linear scale here
rmsd_cut_sigmoid <- predict(objp, list(xp = lambda_cut) )
message(
    'threshold map (sigmoidal): ',
    'lambda = ', lambda_cut, ', ',
    'RMSD = ', signif( rmsd_cut_sigmoid, 3 )
)

# invert relationship, failed :(
# y ~ a * (x^b - 1) / (x^b + 1)
# y * (x^b + 1) ~ a * (x^b - 1)
# y * x^b + y ~ a * x^b - a
# (y-a) * x^b ~ - (y + a)
# x^b ~ (a + y) / (a - y)
# x ~ ( (a + y) / (a - y) )^{1/b}
# x ~ ( (a + y) / (a - y) )^c, c = 1/b
## obj2 <- nls( x ~ ( (a + y) / (a - y) )^c, start = list(a = 1/2, c = 1.6) )
## Error in numericDeriv(form[[3L]], names(ind), env) : 
##   Missing value or an infinity produced when evaluating the model
## obj2 <- nls( x ~ ( (a + y) / (a - y) )^c, start = list(a = 0.6, c = 1) )

# test plot
## plot( x, y, pch = '.', log = 'x' )
## #plot( x, y, pch = '.' )
## lines( x2, predict(obj, list(x = x2) ), col = 2 )
## #lines( x2, fl(x2, 1, 1), col = 2 )
## lines( x2, predict(objp, list(xp = x2) ), col = 3 ) # best fit!

# find common data range
range_rmsd <- range( tib_main$rmsd, na.rm = TRUE )
range_lambda <- range( tib_main$lambda, na.rm = TRUE )

# symetrize lambda range (in log scale)
log_max_lambda <- max( abs( log( range_lambda ) ) )
range_lambda <- exp( c( -log_max_lambda, log_max_lambda ) )
# same with SRMSD, though no log needed here
max_srmsd <- max( abs( range_rmsd ) )
range_rmsd <- c( -max_srmsd, max_srmsd )

# uniformly spaced points for curve
xp <- exp( log_max_lambda * (-100 : 100) / 100 )

# plot in base data dir
fig_start(
    'sum-rmsd-vs-lambda',
    mar_t = 1,
    mar_r = 0.3
)
# try to change default y labeling
par( lab = c(3, 5, 7) )
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
abline( v = 0, lty = lty_guides, col = col_guides )
abline( h = 1, lty = lty_guides, col = col_guides )
# plot reversed to match actual plot (regression had x and y flipped)
lines( predict(objp, list(xp = xp) ), xp, lty = lty_fit_sigmoid, col = col_fit_sigmoid ) # using uniformly spaced points
#lines( predict(objp, list(xp = x2) ), x2, lty = lty_fit_sigmoid, col = col_fit_sigmoid ) # using actual data

# slope prediction at lambda = 1
# (linear approx)
lm_m <- 1 / ( coef( objp )[1] * coef( objp )[2] / 2)
## message( 'linear approx: lambda = RMSD * ', signif( lm_m, 3)  )
message( 'log-linear approx: log(lambda) = RMSD * ', signif( lm_m, 3)  )
## abline(
##     1,
##     lm_m,
##     lty = lty_fit_loglin,
##     col = col_fit_loglin,
##     untf = TRUE # untransform because it's a log plot
## )
# log-linear approx
lines(
    log( xp ) / lm_m,
    xp,
    lty = lty_fit_loglin,
    col = col_fit_loglin
)
# report on the RMSD predicted for lambda = 1.05
# lambda is in linear scale here too
rmsd_cut_loglin <- log( lambda_cut ) / lm_m
message(
    'threshold map (log-linear): ',
    'lambda = ', lambda_cut, ', ',
    'RMSD = ', signif( rmsd_cut_loglin, 3 )
)


# randomize rows so last dataset doesn't just overlap previous datasets
tib_main <- tib_main[ sample( nrow(tib_main) ), ]
# add data on top
points(
    tib_main$rmsd,
    tib_main$lambda,
    col = tib_main$col,
    pch = '.'
)
# legend
legend(
    'bottomright',
    datasets$name_short,
    title = 'Dataset',
    text.col = datasets$col,
    pch = NA,
    bty = 'n',
    cex = 0.5
)
# second legend just for fit
legend(
    'topleft',
    c(
        expression( bold( 'Data' ) ),
        expression( bold( 'Sigmoid fit' ) ),
        expression( bold( paste( 'Log-linear approx at ', lambda == 1 ) ) )
    ),
    col = c('black', col_fit_sigmoid, col_fit_loglin),
    lty = c(NA, lty_fit_sigmoid, lty_fit_loglin),
    pch = c('.', NA),
    bty = 'n',
    cex = 0.5
)
fig_end()

## #############

## # test a hypothesis about generating this curve
## # failed experiment, curves don't track well for high RMSD, lambda, so meh

## # assume the data is coming from a chi-square with different degrees of freedom
## m <- 1000 # samples (numbers of p-vals)
## df <- 1 # in GWAS this is always 1 (genotype vector is added in H1, removed in H0, all else the same)
## ps <- ( 1 : m ) / ( m + 1 ) # percentiles to retrieve, so we barely miss 0 and 1 edges
## x_m <- qchisq( 0.5, df = df ) # median statistic

## # data we want
## n <- 100 # resolution of curve
## # deltas in a log scale
## deltas <- ( -n : n ) / n * 0.8
## # in real scale
## deltas <- 10^( deltas )
## l <- length( deltas )
## lambdas_sim <- vector( 'numeric', l )
## rmsds_sim <- vector( 'numeric', l )
## for ( i in 1 : l ) {
##     message( 'i: ', i )
##     delta <- deltas[i]
##     # this generates the alternate chi-squared quantiles
##     xs <- qchisq( ps, df = delta )
##     # then runs them through the cumulative distribution (standard version) to get the (bad) pvalues
##     pvals <- pchisq( xs, df = df )
##     # lambdas (directly on raw chi-squared stats, not through p-values)
##     lambdas_sim[i] <- median( xs ) / x_m
##     # compute error metric (RMSD)
##     rmsds_sim[i] <- rmsd( pvals, ps )
##     # add sign depending on delta
##     if ( delta < 1 )
##         rmsds_sim[i] <- - rmsds_sim[i]
## }

## lines(
##     rmsds_sim,
##     lambdas_sim
## )

## fig_end()
