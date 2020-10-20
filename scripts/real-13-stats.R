# this script reads big table of AUC and RMSD estimates, produces stats (whether anything is significantly different than the best, per metric).

library(optparse)
library(readr)
library(ochoalabtools)
library(tibble)
library(dplyr)

# constants
# output file name (big table)
file_table <- 'sum.txt'

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--bfile", type = "character", default = NA, 
                help = "Base name for input plink files (.BED/BIM/FAM)", metavar = "character"),
    make_option("--complete", action = "store_true", default = FALSE, 
                help = "Plot only complete replicates (useful for partial runs)"),
    make_option("--cut", type = "double", default = 0.01, 
                help = "P-value threshold for Wilcoxon 2-sample comparisons", metavar = 'double')
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name <- opt$bfile
cut <- opt$cut

# stop if name is missing
if ( is.na(name) )
    stop('`--bfile` terminal option is required!')

# final versions for papers
# pure plink version is only PCA version, add LMM too
methods <- c(
    'pca-plink-pure',
    'gcta'
)
name_base <- 'sum-'

# move to where the data is
setwd( '../data/' )
setwd( name )

# read the big table!
tib <- read_tsv(
    file_table,
    col_types = 'ciiddd'
)

# extract methods from table itself
#methods <- names( method_to_label ) # not from table, but from hardcoded map, always lists PCA first!
# and PCs too, numerically
pcs <- sort( as.numeric( unique( tib$pc ) ) )

# a hack to only plot complete reps
if ( opt$complete ) {
    for ( rep in 1:50 ) {
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
### AUC/RMSD ###
################

get_best_pcs <- function(tib, metric, method) {
    # NOTE: input tibble should already be subset to be for a single method, but for several PCs and replicates
    
    # get summary mean performance per PC
    means <- aggregate( tib[[ metric ]], list(tib$pc), mean)$x

    # this PC was closest to perfection
    if ( metric == 'rmsd' ) {
        # best RMSDs are near zero in absolute value
        pc_best <- which.min( abs( means ) ) - 1
    } else if ( metric == 'auc' ) {
        # best AUCs are largest
        pc_best <- which.max( means ) - 1
    }
    # extract data corresponding to that PC only
    vals_best <- tib[[ metric ]][ tib$pc == pc_best ]
    
    # now we repeat tests again, but this time two-sample tests against the "best" one
    pvals <- vector( 'numeric', length(pcs) )
    for ( pc in pcs ) {
        if ( pc == pc_best ) {
            # self comparison should always carry a p-value of 1
            # actually best to avoid wilcox.test since it gives NaN instead because there's zeroes
            pvals[ pc + 1 ] <- 1
        } else {
            # subset table more
            # now this is just vector of values
            vals <- tib[[ metric ]][ tib$pc == pc ]
            # test if they are significantly different than the best or not
            # use paired version
            pvals[ pc + 1 ] <- wilcox.test( vals, vals_best, paired = TRUE )$p.value
        }
    }

    # this is the smallest PC that performed as well as the best, in the sense that it is not significantly different
    pc_best_min <- min( which( pvals > cut ), na.rm = TRUE ) - 1

    # create output tibble row
    data <- tibble(
        method = method,
        metric = metric,
        best = pc_best,
        min = pc_best_min
    )
    # return that plus distribution of best values, for additional comparisons
    return(
        list(
            data = data,
            vals = vals_best
        )
    )
}

# though this is kinda short, let's put in a tibble
output <- tibble( .rows = 0 )
# other data we need to compare across methods
method_to_vals_rmsd <- list()
method_to_vals_auc <- list()

# begin processing each method separately
for ( method in methods ) {
    # subset big table
    tib_i <- tib[ tib$method == method, ]
    
    # apply function to data, RMSD version
    obj <- get_best_pcs(tib_i, 'rmsd', method)
    method_to_vals_rmsd[[ method ]] <- obj$vals # save separately
    # concatenate output rows
    output <- bind_rows( output, obj$data )

    # repeat for AUCs
    obj <- get_best_pcs(tib_i, 'auc', method)
    method_to_vals_auc[[ method ]] <- obj$vals # save separately
    output <- bind_rows( output, obj$data )
}

# just print report when done
print( output )

# final report compares best case across methods

report_cross_method <- function( method_to_vals, metric ) {
    # decide which performed best, by mean
    method_to_mean <- sapply( method_to_vals, mean )
    if ( metric == 'auc' ) {
        method_best <- methods[ which.max( method_to_mean ) ]
    } else if ( metric == 'rmsd' ) {
        method_best <- methods[ which.min( abs(method_to_mean ) ) ]
    }

    # and decide if they're statistically significantly different or not
    sig <- wilcox.test(
        method_to_vals[[ methods[1] ]],
        method_to_vals[[ methods[2] ]],
        paired = TRUE
    )$p.value < cut

    message( 'best ', metric, ': ', method_best, ' ', if (sig) '(significant)' else '(tie)' )
}

report_cross_method( method_to_vals_rmsd, 'rmsd' )
report_cross_method( method_to_vals_auc, 'auc' )
