# an earlier figure plotted inflation factors vs RMSD in a single dataset, revealing a clean curve
# is the curve the same across datasets with different sample sizes, etc?

library(readr)
library(dplyr) # for bind_rows
library(ochoalabtools)

## # FOR FAILED EXPERIMENT AT END OF SCRIPT
## library(simtrait) # for rmsd

#################
### CONSTANTS ###
#################

# to toggle between fitting the whole data or just the bottom half
fit_top_half_only <- TRUE

# output file name (big table)
file_table <- 'sum.txt'
# directory name, needed in one mode (weird var name just stuck)
dir_phen <- 'fes'
dir_low <- 'm_causal_fac-27/h0.3'
dir_env <- paste0( dir_low, '/env0.3-0.2' )
# output name
name_out_lambda <- 'rmsd-vs-lambda'
name_out_tie <- 'rmsd-vs-tie'
name_out_power <- 'auc-vs-power'
# methods to keep in analysis
methods <- c('pca-plink-pure', 'gcta')
# report on the RMSD predicted for lambda = 1.05
lambda_cut <- 1.05
srmsd_cut <- 0.01

# only thing compared to other plots is we want this alternate order
name_dir_order <- c(
    'sim-n1000-k10-f0.1-s0.5-g1',
    'sim-n100-k10-f0.1-s0.5-g1',
    'sim-n1000-k10-f0.1-s0.5-g20',
    'HoPacAll_ld_prune_1000kb_0.3_maf-0.01',
    'HoPacAll_ld_prune_1000kb_0.3_maf-0.01_sim',
    'hgdp_wgs_autosomes_ld_prune_1000kb_0.3_maf-0.01_geno-0.1',
    'hgdp_wgs_autosomes_ld_prune_1000kb_0.3_maf-0.01_geno-0.1_sim',
    'tgp-nygc-autosomes_ld_prune_1000kb_0.3_maf-0.01',
    'tgp-nygc-autosomes_ld_prune_1000kb_0.3_maf-0.01_sim'
)

# color for fit curve
col_fit_sigmoid <- 'gray40'
lty_fit_sigmoid <- 2
# and guide lines
col_guides <- 'gray95'
lty_guides <- 1

############
### DATA ###
############

# move to where the data is
setwd( '../data/' )

# read datasets info (names for inputs and output, colors, line types)
datasets <- read_tsv( 'datasets.txt', col_types = 'cccii' )
# reorder as desired
indexes <- match( name_dir_order, datasets$name_dir )
datasets <- datasets[ indexes, ]

# big table of interest
# initialize this way, it'll grow correctly
data <- NULL

# load each dataset
for ( i in 1 : nrow( datasets ) ) {
    # enter dir
    setwd( datasets$name_dir[ i ] )

    for ( fes in c(FALSE, TRUE) ) {
        # move in one more level in this case
        if ( fes )
            setwd( dir_phen )

        for ( trait in c('hi', 'lo', 'env') ) {
            # hi means do nothing
            dir_trait <- ''
            if ( trait == 'lo' ) {
                dir_trait <- dir_low
            } else if ( trait == 'env' ) {
                dir_trait <- dir_env
            }
            if ( dir_trait != '' ) {
                # stop if it doesn't exist (won't for real-sim)
                if ( !dir.exists( dir_trait ) )
                    next
                # to get back easily
                dir_pre_trait <- getwd()
                # actually go there
                setwd( dir_trait )
            }
            
            # read the big table!
            tib <- read_tsv(
                file_table,
                col_types = 'ciiddd'
            )

            # subset to use only the two methods we talk about in the paper
            tib <- tib[ tib$method %in% methods, ]

            # recall the dataset of origin
            tib$dataset <- datasets$name_paper[ i ]
            # and the trait simulation type (in a shorter-hand notation)
            tib$trait <- if ( fes ) 'inv' else 'rand'
            
            ## # color by method
            ## # (boring idea, obviously GCTA was always negative and PCA was always positive SRMSD)
            ## tib$col <- 'blue' # default
            ## tib$col[ tib$method == methods[2] ] <- 'red' # GCTA is red

            # color by dataset
            tib$col <- datasets$col[ i ]
            
            # concatenate into bigger table
            data <- bind_rows( data, tib )

            # get back if needed
            if ( dir_trait != '' )
                setwd( dir_pre_trait )
        }
        
        # move back one more level in this case
        if ( fes )
            setwd( '..' )
    }
    # go back down
    setwd( '..' )
}

# randomize rows so last dataset doesn't just overlap previous datasets
data <- data[ sample( nrow(data) ), ]

###########################
### MISSINGNESS FILTERS ###
###########################

# remove all NAs, can't plot or fit

# power has slightly higher missingness than the rest of the stats
indexes_missing <- is.na( data$power_calib2 )
# subset now
data <- data[ !indexes_missing, ]
# there should be no more missingness after this
stopifnot( !anyNA( data ) )

# also remove infinities (observed once for lambda only)
indexes_finite <- is.finite( data$lambda )
data <- data[ indexes_finite, ]

######################
### LAMBDA vs RMSD ###
######################

# sigmoid fit
# NOTE: x and y are reversed from actual plot
y <- data$rmsd
#x <- log(data$lambda)
x <- data$lambda
x2 <- sort(x)

if ( fit_top_half_only ) {
    # fit on top data only
    indexes <- x > 1
    xp <- x[ indexes ]
    yp <- y[ indexes ]
} else {
    # use all data (copy to the same variables so downstream code works both ways)
    xp <- x
    yp <- y
}
objp <- nls( yp ~ a * (xp^b - 1) / (xp^b + 1), start = list(a = 1/2, b = 1) )
print( coef( objp ) )

# extract coefficients for predictions
a_fit <- coef( objp )[1]
b_fit <- coef( objp )[2]

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
# this is the inverse function
srmsd_to_lambda <- function ( y, a, b ) {
    ( (a + y) / (a - y) ) ^ ( 1 / b )
}
lambda_cut_by_srmsd <- srmsd_to_lambda( srmsd_cut, a_fit, b_fit )
message(
    'Inverse threshold map (sigmoidal): ',
    'RMSD = ', srmsd_cut, ', ',
    'lambda = ', signif( lambda_cut_by_srmsd, 3 )
)


# test plot
## plot( x, y, pch = '.', log = 'x' )
## #plot( x, y, pch = '.' )
## lines( x2, predict(obj, list(x = x2) ), col = 2 )
## #lines( x2, fl(x2, 1, 1), col = 2 )
## lines( x2, predict(objp, list(xp = x2) ), col = 3 ) # best fit!

###################
### PLOT LAMBDA ###
###################

# find common data range
range_rmsd <- range( data$rmsd, na.rm = TRUE )
range_lambda <- range( data$lambda, na.rm = TRUE )

# symetrize lambda range (in log scale)
log_max_lambda <- max( abs( log( range_lambda ) ) )
range_lambda <- exp( c( -log_max_lambda, log_max_lambda ) )
# same with SRMSD, though no log needed here
max_srmsd <- max( abs( range_rmsd ) )
range_rmsd <- c( -max_srmsd, max_srmsd )

# uniformly spaced points for curve
xp <- exp( log_max_lambda * (-100 : 100) / 100 )

# plotting labels, shared across plots
lab_rmsd <- expression( bold( SRMSD[p] ) )

# plot in base data dir
width <- fig_width() / 2
fig_start(
    name_out_lambda,
    mar_t = 1,
    mar_r = 0.3,
    width = width,
    height = width
)
# try to change default y labeling
par( lab = c(3, 5, 7) )
# start base plot
plot(
    NA,
    xlim = range_rmsd,
    ylim = range_lambda,
    xlab = lab_rmsd,
    ylab = 'Inflation Factor',
    log = 'y'
)
# guide lines
abline( v = 0, lty = lty_guides, col = col_guides )
abline( h = 1, lty = lty_guides, col = col_guides )
# plot reversed to match actual plot (regression had x and y flipped)
lines( predict(objp, list(xp = xp) ), xp, lty = lty_fit_sigmoid, col = col_fit_sigmoid ) # using uniformly spaced points
#lines( predict(objp, list(xp = x2) ), x2, lty = lty_fit_sigmoid, col = col_fit_sigmoid ) # using actual data


# add data on top
points(
    data$rmsd,
    data$lambda,
    col = data$col,
    pch = '.'
)
# legend
legend(
    'bottomright',
    datasets$name_paper,
    title = 'Dataset',
    text.col = datasets$col,
    pch = NA,
    bty = 'n',
    cex = 0.7
)
# second legend just for fit
legend(
    'topleft',
    c(
        expression( bold( 'Data' ) ),
        expression( bold( 'Sigmoid fit' ) )
    ),
    col = c('black', col_fit_sigmoid),
    lty = c(NA, lty_fit_sigmoid),
    pch = c('.', NA),
    bty = 'n',
    cex = 0.7
)
fig_end()

## #############

## # NOTE: should have used non-centrality parameters!  (instead of different degree of freedom!)

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

####################
### TYPE I ERROR ###
####################

# log: 2 (1e-2) looks perfect, the rest start to look crummy because of all the zeroes
# linear: all look hard to interpret, meh
# so stick with log version of 1e-2
log10_alpha <- 2

# same dimensions as lambda fig
fig_start(
    name_out_tie,
    mar_t = 1,
    mar_r = 0.3,
    width = width,
    height = width
)
plot(
    data$rmsd,
    data[[ paste0( 'type_1_err', log10_alpha ) ]],
    col = data$col,
    pch = '.',
    xlab = lab_rmsd,
    ylab = 'Type I error rate',
    log = 'y'
)
abline( h = 10^(-log10_alpha), lty = 2, col = 'gray' )
abline( v = 0, lty = 2, col = 'gray' )
# legend
legend(
    'bottomright',
    datasets$name_paper,
    title = 'Dataset',
    text.col = datasets$col,
    pch = NA,
    bty = 'n',
    cex = 0.7
)
fig_end()


#############
### POWER ###
#############

# force symmetric and leave space for legend
ylim <- range( data$power_calib4 )

# same dimensions as lambda fig
fig_start(
    name_out_power,
    mar_t = 1,
    mar_r = 0.3,
    width = width,
    height = width
)
plot(
    data$auc,
    data$power_calib4,
    col = data$col,
    pch = '.',
    xlab = expression( bold( AUC[PR] ) ),
    ylab = 'Calibrated power',
    xlim = ylim,
    ylim = ylim
)
abline( 0, 1, lty = 2, col = 'gray' )
# legend
legend(
    'bottomright',
    datasets$name_paper,
    title = 'Dataset',
    text.col = datasets$col,
    pch = NA,
    bty = 'n',
    cex = 0.7
)
fig_end()
