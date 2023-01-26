# this script reads big table of AUC and RMSD estimates, produces stats (whether anything is significantly different than the best, per metric).

library(readr)
library(ochoalabtools)
library(tibble)
library(dplyr)
library(optparse)

# constants
# output file name (big table)
file_table <- 'sum.txt'
# P-value threshold for Wilcoxon 2-sample paired 1-tail comparisons
# NOTE: becomes Bonferroni-corrected below
p_cut <- 0.01
# threshold at which SRMSDs are essentially good enough
srmsd_cut <- 0.01
# hardcode for all datasets
pcs <- 0L : 90L
# for a length validation
n_reps <- 50L
# number of trait models (FES, RC), hardcoded, for Bonferroni correction
n_traits <- 2L

# final versions for papers
# pure plink version is only PCA version, add LMM too
methods <- tibble(
    code = c( 'pca-plink-pure', 'gcta' ),
    name = c( 'PCA', 'LMM' )
)
metrics <- c('rmsd', 'auc')

############
### ARGV ###
############

# define options
option_list = list(
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

# get values
m_causal_fac <- opt$m_causal_fac
herit <- opt$herit
env1 <- opt$env1
env2 <- opt$env2
labs <- opt$labs

if ( labs ) {
    # add a third method
    methods <- bind_rows( methods, tibble( code = 'gcta-labs', name = 'LMM lab.' ) )
}

##################
### STAT TESTS ###
##################

get_best_pcs <- function(tib, metric) {
    # NOTE: input tibble should already be subset to be for a single method, but for several PCs and replicates per PC
    # extract this column we call very often
    tib_metric <- tib[[ metric ]]

    # get summary mean performance per PC
    obj_agg <- aggregate( tib_metric, list(tib$pc), mean)
    # make sure means are aligned to `pcs`
    stopifnot( obj_agg$pc == pcs )
    # get the only data we wanted out of this aggregation
    means <- obj_agg$x
    
    # this PC was closest to perfection
    if ( metric == 'rmsd' ) {
        # best RMSDs are near zero in absolute value
        r <- pcs[ which.min( means ) ]
    } else if ( metric == 'auc' ) {
        # best AUCs are largest
        r <- pcs[ which.max( means ) ]
    }
    # extract data corresponding to that PC only
    vals <- tib_metric[ tib$pc == r ]
    if( length( vals ) != n_reps )
        stop( 'Length of `vals` not ', n_reps, ' for ', metric, ', r = ', r )
    
    # return that plus distribution of best values, for additional comparisons
    return(
        list(
            r = r,
            vals = vals
        )
    )
}

# applies overall statistical test between two distributions
test_distrs <- function( vals, vals_best, metric ) {

    # defaults for trivial case
    pval <- 1
    reverse <- FALSE

    # identify a common trivial case
    # often the two distributions to compare are actually the same (particularly LMM r=0 case)
    # self comparison always results in p-value of 1
    # avoid wilcox.test since it gives NaN instead because there's zeroes
    if ( !all( vals == vals_best, na.rm = TRUE ) ) {

        # we got lazy in coding, always test worst vs best regardless of order provided
        if ( metric == 'rmsd' ) {
            alternative <- 'g' # alternative is greater than min
            # this says `vals` is actually better than `vals_best`
            if ( mean( vals_best, na.rm = TRUE ) > mean( vals, na.rm = TRUE ) )
                reverse <- TRUE
        } else if ( metric == 'auc' ) {
            alternative <- 'l' # alternative is lesser than max
            # this says `vals` is actually better than `vals_best`
            if ( mean( vals_best, na.rm = TRUE ) < mean( vals, na.rm = TRUE ) )
                reverse <- TRUE
        }

        # apply reversal if needed
        if ( reverse ) {
            vals_tmp <- vals
            vals <- vals_best
            vals_best <- vals_tmp
        }
        
        # test if they are significantly different than the best or not
        # use paired version
        pval <- wilcox.test( vals, vals_best, paired = TRUE, alternative = alternative )$p.value
    }
    
    # return p-value and whether data was reversed or not
    return(
        list(
            pval = pval,
            reverse = reverse
        )
    )
}

test_calib <- function( vals, metric ) {
    calib <- NA # keep it this way for AUC
    if ( metric == 'rmsd' ) 
        calib <- mean( vals, na.rm = TRUE ) < srmsd_cut
    return( calib )
}

# gathers more loops essentially and massages data to look as needed for output table
process_dataset <- function( name, fes ) {
    ######################
    ### LOAD / PROCESS ###
    ######################
    
    # read the big table!
    # load from wherever current location is
    tib <- read_tsv(
        file_table,
        col_types = 'ciiddd'
    )

    # translate method names to nicer ones for paper
    # replace method names for paper
    for ( i in 1 : nrow( methods ) )
        tib$method[ tib$method == methods$code[ i ] ] <- methods$name[ i ]
    
    # RMSD data is signed, but in the comparisons of this entire script we don't care about the sign/direction
    # for simplicity, take absolute value now for all data (all PCs and all methods)!
    tib$rmsd <- abs( tib$rmsd )

    #############
    ### TESTS ###
    #############
    
    # separate LMM and PCA data
    # (each all pcs, reps, metrics, so still a tibble)
    tib_lmm <- tib[ tib$method == 'LMM', ]
    tib_pca <- tib[ tib$method == 'PCA', ]
    
    # get data for LMM r=0
    # (all reps, metrics, so still a tibble)
    tib_lmm_0 <- tib_lmm[ tib_lmm$pc == 0, ]

    if ( labs ) {
        # one more set of tests to add!
        # (all reps, metrics, so still a tibble)
        tib_lmm_labs_0 <- tib[ tib$method == 'LMM lab.' & tib$pc == 0, ]
    }

    # output tibble to grow
    output <- NULL

    # begin processing each metric separately
    for ( metric in metrics ) {
        # get data for r=0, now a vector
        lmm_0_vals <- tib_lmm_0[[ metric ]]
        
        # identify best r for LMM, which is metric-dependent
        obj <- get_best_pcs( tib_lmm, metric )
        lmm_b_vals <- obj$vals
        lmm_b_r <- obj$r
        
        # find best PCA model too (here never consider r=0)
        obj <- get_best_pcs( tib_pca, metric )
        pca_b_vals <- obj$vals
        pca_b_r <- obj$r

        # perform some tests
        # LMM best vs LMM r=0
        obj <- test_distrs( lmm_0_vals, lmm_b_vals, metric )
        # by definition shouldn't have reversals here
        stopifnot( ! obj$reverse )
        lmm_b_lmm_0_pval <- obj$pval

        # PCA best vs LMM r=0
        # here we don't know which one is actually better
        obj <- test_distrs( pca_b_vals, lmm_0_vals, metric )
        pca_b_lmm_0_pval <- obj$pval
        pca_b_lmm_0_best <- if ( obj$reverse ) 'PCA' else 'LMM'

        # PCA best vs LMM best
        obj <- test_distrs( pca_b_vals, lmm_b_vals, metric )
        pca_b_lmm_b_pval <- obj$pval
        pca_b_lmm_b_best <- if ( obj$reverse ) 'PCA' else 'LMM'
        
        # store all data of interest
        output_row <- tibble(
            name_paper = name,
            trait = if ( fes ) 'FES' else 'RC',
            metric = metric,
            # is lmm r=0 calibrated?
            lmm_0_calib = test_calib( lmm_0_vals, metric ),
            # determine if lmm with best r is better than r=0 or not
            lmm_b_r = lmm_b_r,
            lmm_b_calib = test_calib( lmm_b_vals, metric ),
            lmm_b_lmm_0_pval = lmm_b_lmm_0_pval,
            lmm_b_lmm_0_sig = lmm_b_lmm_0_pval < p_cut,
            # compare lmm r=0 to pca with best r too
            pca_b_r = pca_b_r,
            pca_b_calib = test_calib( pca_b_vals, metric ),
            pca_b_lmm_0_pval = pca_b_lmm_0_pval,
            pca_b_lmm_0_best = pca_b_lmm_0_best,
            pca_b_lmm_0_sig = pca_b_lmm_0_pval < p_cut,
            # lastly, overkill comparison of best lmm to best pca (only non-redundant under env, given previous results)
            pca_b_lmm_b_pval = pca_b_lmm_b_pval,
            pca_b_lmm_b_best = pca_b_lmm_b_best,
            pca_b_lmm_b_sig = pca_b_lmm_b_pval < p_cut
        )

        if ( labs ) {
            # one more set of tests to add!
            # get values for this one distribution (r=0 is only choice)
            lmm_labs_0_vals <- tib_lmm_labs_0[[ metric ]]

            # because there are now multiple comparators, lets compare against best of rest only
            # let's narrow the choices to lmm r=0 vs pca best r
            # actually the best choice is already known, let's use that previous calculation
            best_vals <- if ( pca_b_lmm_0_best == 'PCA' ) pca_b_vals else lmm_0_vals
            # now perform test!
            obj <- test_distrs( best_vals, lmm_labs_0_vals, metric )
            best_lmm_labs_0_pval <- obj$pval
            best_lmm_labs_0_best <- if ( obj$reverse ) pca_b_lmm_0_best else 'LMM lab.'

            # add to report
            output_row$lmm_labs_0_calib <- test_calib( lmm_labs_0_vals, metric )
            output_row$best_lmm_labs_0_pval <- best_lmm_labs_0_pval
            output_row$best_lmm_labs_0_best <- best_lmm_labs_0_best
            output_row$best_lmm_labs_0_sig <- best_lmm_labs_0_pval < p_cut
            # ties are very messy here, let's do the extra work now because its easier at this point
            if ( output_row$best_lmm_labs_0_sig ) {
                # if difference is significant and "LMM lab." was the winner, stays as-is
                # else it could be a tie...
                # if not tie in first round, then report is correct already
                if ( best_lmm_labs_0_best != 'LMM lab.' && ! output_row$pca_b_lmm_0_sig )
                    output_row$best_lmm_labs_0_best <- 'PCA/LMM' # this is correct tie status
            } else {
                # here new test was not significant, but was previous eval a tie?
                if ( output_row$pca_b_lmm_0_sig ) {
                    # this is correct pairwise tie
                    output_row$best_lmm_labs_0_best <- paste0( pca_b_lmm_0_best, '/LMM lab.' )
                } else
                    # mark 3-way tie this way:
                    output_row$best_lmm_labs_0_best <- 'Tie'
            }
        }

        # append to output table
        output <- bind_rows( output, output_row )
    }
    
    return( output )
}


############
### DATA ###
############

# move to where the data is
setwd( '../data/' )

# read dataset paths and names for output
datasets <- read_tsv( 'datasets.txt', col_types = 'cccii' )

# the new evaluations (low herit and env) were not performed on real-sim (tree) datasets!
if ( m_causal_fac != 10 || herit != 0.8 || !is.na( env1 ) )
    datasets <- datasets[ datasets$type != 'Tree', ]

# number of tests in this table, for Bonferroni
n_tests <- nrow( datasets ) * length( methods$name ) * length( metrics ) * n_traits
# adjust p_cut accordingly
p_cut <- p_cut / n_tests

# to load real datasets and get back down easily
dir_orig <- getwd()

# path where data will be (for each dataset)
dir_out <- ''
if ( m_causal_fac != 10 )
    dir_out <- paste0( dir_out, 'm_causal_fac-', m_causal_fac, '/' )
if ( herit != 0.8 )
    dir_out <- paste0( dir_out, 'h', herit, '/' )
if ( !is.na( env1 ) )
    dir_out <- paste0( dir_out, 'env', env1, '-', env2, '/' )

# tibble for all datasets
output <- tibble( .rows = 0 )

for ( i in 1 : nrow( datasets ) ) {
    message( datasets$name_dir[ i ], ' RC' )
    setwd( datasets$name_dir[ i ] )
    dir_dataset <- getwd()
    if ( dir_out != '' )
        setwd( dir_out )
    name_paper <- datasets$name_paper[ i ]

    # load local data, calculate row, add to final output immediately
    output <- bind_rows(
        output, 
        process_dataset( name_paper, fes = FALSE )
    )
    
    # now do `fes = TRUE` case
    message( datasets$name_dir[ i ], ' FES' )
    setwd( dir_dataset )
    setwd( 'fes' )
    if ( dir_out != '' )
        setwd( dir_out )
    
    # load local data, calculate row, add to final output immediately
    output <- bind_rows(
        output, 
        process_dataset( name_paper, fes = TRUE )
    )
    
    # go back down
    setwd( dir_orig )
}

# validate our Bonferroni calculation
# (this table isn't complete until after we've applied all thresholds, that's why we calculate it first, then validate)
# only variable in columns rather than rows are the methods
stopifnot( nrow( output ) * length( methods$name ) == n_tests )
message( 'Number of tests (for Bonferroni): ', n_tests )

################
### CLEANUPS ###
################

# clean up tie situations (makes words easier to scan)
output$pca_b_lmm_0_best[ !output$pca_b_lmm_0_sig ] <- 'Tie'
output$pca_b_lmm_b_best[ !output$pca_b_lmm_b_sig ] <- 'Tie'
## # for labs case, tie is ambiguous because it's now 3 methods, so let's be more explicit
## if ( labs ) {
##     output$best_lmm_labs_0_best[ !output$best_lmm_labs_0_sig ] <- 'Tie'
## }

# more automatic reductions when LMM r=0 is best across the board
# NOTE: output$lmm_b_lmm_0_sig gets deleted later for all cases
# Update: let's always omit these cases, even in env they were largely the same as for LMM with r=0
#if ( all( !output$lmm_b_lmm_0_sig ) ) {
    # who cares if lmm_b is calibrated
    output$lmm_b_calib <- NULL
    # don't show lmm_b vs pca comparisons
    output$pca_b_lmm_b_best <- NULL
    output$pca_b_lmm_b_pval <- NULL
    output$pca_b_lmm_b_sig <- NULL
#}

########################
### PAPER FORMATTING ###
########################

# only function of next edits is to make table prettier, as it gets added to paper automatically

# reorder so all FES traits are listed first, then rmsd before auc
output <- arrange( output, trait, desc(metric) )

# round p-values and mark significant ones with asterisks
for ( colname in c('lmm_b_lmm_0', 'pca_b_lmm_0', 'pca_b_lmm_b', 'best_lmm_labs_0') ) {
    # get two column names to process
    colname_pval <- paste0( colname, '_pval' )
    colname_sig <- paste0( colname, '_sig' )
    # copy down values for ease
    pvals <- output[[ colname_pval ]]
    sigs <- output[[ colname_sig ]]
    # silently skip a case that might be missing already
    if ( is.null( pvals ) )
        next
    # round p-values
    pvals <- signif( pvals, 3 )
    # add asterisks, denoting significance, after rounding
    pvals[ sigs ] <- paste0( pvals[ sigs ], '*' )
    # copy pvals back
    output[[ colname_pval ]] <- pvals
    # delete indicator column now
    output[[ colname_sig ]] <- NULL
}

# metrics replace with LaTeX code
output$metric[ output$metric == 'rmsd' ] <- '$|\\rmsd|$'
output$metric[ output$metric == 'auc' ] <- '$\\auc$'

# uncapitalize booleans
# sadly very tedious
for ( colname in c('lmm_0_calib', 'lmm_b_calib', 'pca_b_calib', 'lmm_labs_0_calib') ) {
    # extract column
    x <- output[[ colname ]]
    # ignore columns that aren't present
    if ( is.null( x ) )
        next
    # edit
    x[ x == TRUE ] <- 'True'
    x[ x == 'FALSE' ] <- 'False' # cause now they're all strings, right?
    x[ is.na( x ) ] <- '' # looks better blank! (also NA are not stringified, stays proper NA)
    # write back
    output[[ colname ]] <- x
}

############
### SAVE ###
############

# go where output will be, path of data except RC version
if ( dir_out != '' )
    setwd( dir_out )

# write table to a file
write_tsv( output, file = 'stats.txt' )
