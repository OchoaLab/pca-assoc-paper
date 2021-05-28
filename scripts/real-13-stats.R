# this script reads big table of AUC and RMSD estimates, produces stats (whether anything is significantly different than the best, per metric).

library(readr)
library(ochoalabtools)
library(tibble)
library(dplyr)

# constants
# output file name (big table)
file_table <- 'sum.txt'
# P-value threshold for Wilcoxon 2-sample paired 1-tail comparisons
p_cut <- 0.01
# threshold at which SRMSDs are essentially good enough
srmsd_cut <- 0.01
# hardcode for all datasets
pcs <- 0 : 90

# final versions for papers
# pure plink version is only PCA version, add LMM too
methods <- tibble(
    code = c( 'pca-plink-pure', 'gcta' ),
    name = c( 'PCA', 'LMM' )
)
metrics <- c('rmsd', 'auc')


##################
### STAT TESTS ###
##################

get_best_pcs <- function(tib, metric, method) {
    # NOTE: input tibble should already be subset to be for a single method, but for several PCs and replicates per PC
    # extract this column we call very often
    tib_metric <- tib[[ metric ]]

    # get summary mean performance per PC
    obj_agg <- aggregate( tib_metric, list(tib$pc), mean)
    # make sure means are aligned to `pcs`
    stopifnot( obj_agg$pc == pcs )
    # get the only data we wanted out of this aggregation
    means <- obj_agg$x
    
    # an effect size assessment
    eff_size <- '' # for pc_best
    eff_size_min <- '' # for pc_best_min
    
    # this PC was closest to perfection
    if ( metric == 'rmsd' ) {
        # best RMSDs are near zero in absolute value
        pc_best <- pcs[ which.min( means ) ]
        alternative <- 'g' # alternative is greater than min
        # in this case we want to know if the mean was below a given threshold
        if ( min( means, na.rm = TRUE ) < srmsd_cut )
            eff_size <- '*'
    } else if ( metric == 'auc' ) {
        # best AUCs are largest
        pc_best <- pcs[ which.max( means ) ]
        alternative <- 'l' # alternative is lesser than max
    }
    # extract data corresponding to that PC only
    vals_best <- tib_metric[ tib$pc == pc_best ]
    
    # now we repeat tests again, but this time two-sample tests against the "best" one
    pvals <- vector( 'numeric', length(pcs) )
    for ( pc in pcs ) {
        if ( pc == pc_best ) {
            # self comparison should always carry a p-value of 1
            # actually best to avoid wilcox.test since it gives NaN instead because there's zeroes
            pvals[ pcs == pc ] <- 1
        } else {
            # subset table more
            # now this is just vector of values
            vals <- tib_metric[ tib$pc == pc ]
            # test if they are significantly different than the best or not
            # use paired version
            pvals[ pcs == pc ] <- wilcox.test( vals, vals_best, paired = TRUE, alternative = alternative )$p.value
        }
    }

    # this is the smallest PC that performed as well as the best, in the sense that it is not significantly different
    pc_best_min <- pcs[ min( which( pvals > p_cut ), na.rm = TRUE ) ]
    # figure out effect size mark for this case
    # do only for RMSD!
    if ( metric == 'rmsd' ) {
        if ( means[ pcs == pc_best_min ] < srmsd_cut )
            eff_size_min <- '*'
    }
    # extract data corresponding to pc_best_min only
    vals_best_min <- tib_metric[ tib$pc == pc_best_min ]

    # before returning, collapse pc_best and pc_best_min (often the same)
    # and add sRMSD effect size mark
    if ( pc_best != pc_best_min ) {
        pc_best <- paste0( pc_best, eff_size, ' (', pc_best_min, eff_size_min, ')' )
    } else {
        # just add effect size mark
        pc_best <- paste0( pc_best, eff_size )
    }

    # return that plus distribution of best values, for additional comparisons
    return(
        list(
            pc_best = pc_best,
            vals_best = vals_best, # in paper it makes more sense to talk about these
            vals_best_min = vals_best_min # in paper it makes more sense to talk about these
        )
    )
}

# final report compares best case across methods

report_cross_method <- function( method_to_vals, metric ) {
    # decide which performed best, by mean
    method_to_mean <- sapply( method_to_vals, mean )
    if ( metric == 'auc' ) {
        method_best <- methods$code[ which.max( method_to_mean ) ]
        alternative <- 'l' # alternative is lesser than max
    } else if ( metric == 'rmsd' ) {
        method_best <- methods$code[ which.min( method_to_mean ) ]
        alternative <- 'g' # alternative is greater than min
    }
    # to get direction right, need to identify worst method too
    method_worst <- setdiff( methods$code, method_best )
    
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

# gathers more loops essentially and massages data to look as needed for output table
process_dataset <- function( name, const_herit_loci ) {
    # read the big table!
    # load from wherever current location is
    tib <- read_tsv(
        file_table,
        col_types = 'ciiddd'
    )
    
    # RMSD data is signed, but in the comparisons of this entire script we don't care about the sign/direction
    # for simplicity, take absolute value now for all data (all PCs and all methods)!
    tib$rmsd <- abs( tib$rmsd )

    # for more processing per metric
    output_metric <- list()

    # begin processing each metric separately
    for ( metric in metrics ) {
        # other data we need to compare across methods
        method_to_vals_best <- list()
        method_to_vals_best_min <- list()
        method_to_pc <- list()
        
        # begin processing each method separately
        for ( method in methods$code ) {
            # subset big table
            tib_i <- tib[ tib$method == method, ]
            
            # apply function to data, RMSD version
            obj <- get_best_pcs(tib_i, metric, method)
            method_to_vals_best[[ method ]] <- obj$vals_best # save separately
            method_to_vals_best_min[[ method ]] <- obj$vals_best_min # save separately
            method_to_pc[[ method ]] <- obj$pc_best # save separately
        }
        
        # decide on best method!
        best_method <- report_cross_method( method_to_vals_best, metric )
        best_method_min <- report_cross_method( method_to_vals_best_min, metric )
        # combine into one scalar (most of the time these are identical)
        if ( best_method_min != best_method )
            best_method <- paste0( best_method, ' (', best_method_min, ')' )
        
        # reorganize data as we have it in the paper
        # NOTE: PCA is 1st, LMM is 2nd (use shorthand here to access methods)
        # put in a list per metric
        output_metric[[ metric ]] <- tibble(
            r_pca = method_to_pc[[ 1 ]],
            r_lmm = method_to_pc[[ 2 ]],
            method = best_method
        )
    }

    # massage table even more
    # right now we only have two rows, separate them
    output_rmsd <- output_metric[[ 'rmsd' ]]
    output_auc <- output_metric[[ 'auc' ]]

    # edit column names to label each metric
    names( output_rmsd ) <- paste0( 'rmsd_', names( output_rmsd ) )
    names( output_auc ) <- paste0( 'auc_', names( output_auc ) )

    # concatenate back!
    output <- bind_cols(
        tibble(
            name_paper = name,
            trait = if ( const_herit_loci ) 'FES' else 'RC'
        ),
        output_rmsd,
        output_auc
    )

    return( output )
}


############
### DATA ###
############

# move to where the data is
setwd( '../data/' )

# read dataset paths and names for output
datasets <- read_tsv( 'datasets.txt', col_types = 'cccii' )

# tibble for all datasets
output <- tibble( .rows = 0 )

for ( i in 1 : nrow( datasets ) ) {

    setwd( datasets$name_dir[ i ] )
    name_paper <- datasets$name_paper[ i ]

    # load local data, calculate row, add to final output immediately
    output <- bind_rows(
        output, 
        process_dataset( name_paper, const_herit_loci = FALSE )
    )
    
    # now do `const_herit_loci = TRUE` case
    setwd( 'const_herit_loci' )
    
    # load local data, calculate row, add to final output immediately
    output <- bind_rows(
        output, 
        process_dataset( name_paper, const_herit_loci = TRUE )
    )
    
    # go back down
    setwd( '../..' )
}

# replace method names for paper
for ( i in 1 : nrow( methods ) )
    output[ output == methods$code[ i ] ] <- methods$name[ i ]

# reorder so all Fixed Effect Size traits are listed first
output <- arrange( output, trait )

# write table to a file
write_tsv( output, file = 'stats.txt' )
