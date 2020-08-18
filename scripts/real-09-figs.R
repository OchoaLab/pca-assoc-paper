# this script reads big table of AUC and RMSD estimates

library(optparse)
library(scales) # for transparency
library(readr)
library(ochoalabtools)

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
    make_option("--pca", action = "store_true", default = FALSE, 
                help = "Compare PCA versions (PCs from R popkinsuppl::kinship_std vs pure plink)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name <- opt$bfile
pca_test <- opt$pca

# stop if name is missing
if ( is.na(name) )
    stop('`--bfile` terminal option is required!')

# a temporary boolean to compare PCA versions
# (PCs from R popkinsuppl::kinship_std vs pure plink)
if (pca_test) {
    method_to_label <- list(
        'pca-plink' = 'Fixed effects (PCA)',
        'pca-plink-pure' = 'Fixed effects (PCA) pure'
    )
    name_base <- 'sum-pca-'
} else {
    # final versions for papers
    # pure plink version is only PCA version, add LMM too
    method_to_label <- list(
        'pca-plink-pure' = 'Fixed effects (PCA)',
        gcta = 'Mixed effects (LMM+PCA)'
    )
    name_base <- 'sum-'
}
# hardcoded same order as method_to_label
# whether pca_test or not, there's only two things to look at
method_cols <- c(
    'red',
    'blue'
)

# move to where the data is
setwd( '../data/' )
setwd( name )

# read the big table!
tib <- read_tsv(
    file_table,
    col_types = 'ciiddd'
)

# extract methods from table itself
methods <- names( method_to_label ) # not from table, but from hardcoded map, always lists PCA first!
# and PCs too, numerically
pcs <- sort( as.numeric( unique( tib$pc ) ) )

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

################
### AUC/RMSD ###
################

# default legend position
legend_pos <- 'topright'
# hack for small sample size simulation only
if (name == 'sim-n100-k10-f0.1-s0.5-g1' && !pca_test)
    legend_pos <- 'bottomleft'

# gather data into lists, best for boxplots
data_rmsd <- list()
data_auc <- list()
for ( method in methods ) {
    # sublists we want
    data_rmsd_i <- list()
    data_auc_i <- list()
    # subset big table
    tib_i <- tib[ tib$method == method, ]
    for ( pc in pcs ) {
        # subset table more
        tib_ij <- tib_i[ tib_i$pc == pc, ]
        # transfer vectors to lists
        # (shift index, since R doesn't like zero index)
        data_rmsd_i[[ pc+1 ]] <- tib_ij$rmsd
        data_auc_i[[ pc+1 ]] <- tib_ij$auc
    }
    # copy sublists to big lists
    data_rmsd[[ method ]] <- data_rmsd_i
    data_auc[[ method ]] <- data_auc_i
}

# plotting labels
lab_rmsd <- expression( bold( SRMSD[p] ) )
lab_auc <- expression( bold( AUC[PR] ) )
lab_lambda <- expression( bold( paste("Inflation Factor (", lambda, ")") ) )

# global vars:
# - methods
# - method_to_label
# - method_cols
boxplots_rmsd_auc <- function(
                              name_out,
                              data_rmsd,
                              data_auc,
                              r_max,
                              legend_pos = 'topright',
                              alpha = 0.5
                              ) {
    # get common range of data plotted
    # always include zero in range for both
    range_rmsd <- range( 0, unlist(data_rmsd), na.rm = TRUE)
    range_auc <- range( 0, unlist(data_auc), na.rm = TRUE)
    
    # make labels
    # blank all non-multiples of 10, for plot
    rs <- 0 : r_max
    rs[ rs %% 10 != 0 ] <- NA
    # overlaid panels have all blanks (R doesn't handle this correctly otherwise)
    rs_blank <- rep.int( NA, r_max+1)
    
    # add transarency to colors!!!
    method_cols_alpha <- alpha(method_cols, alpha)
    
    # other shared params (same as standard boxplot)
    outline <- FALSE # no outliers plotted separately (so busy as it is)
    #  range <- 0 # make whiskers extend to full range
    range <- 1.5 # default whiskers range
    whisklty <- 1 # whisker line type (default 2?)
    
    # start PDF
    fig_start(name_out, width = 7, height = 5, mar_b = 1.5)
    # add lower margin, so inner margins can be even smaller
    par( oma = c(1.5, 0, 0, 0) )
    # two panels
    par( mfrow = c(2, 1) )
    # get this param after fig_start (which modifies lwd)
    medlwd <- par('lwd') # median line width (default is 3 times this)
    # boxplots!
    # top panel
    for ( method in methods ) {
        # get index to match colors up
        i <- which( methods == method )
        # first method has some special cases
        is_first <- i == 1
        # plot data!
        boxplot(
            data_rmsd[[ method ]],
            add = !is_first, # add all but first time
            names = if ( is_first ) rs else rs_blank,
            xlab = "",
            ylab = if ( is_first ) lab_rmsd else '',
            ylim = range_rmsd,
            border = method_cols_alpha[ i ],
            col = NA,
            outline = outline,
            range = range,
            whisklty = whisklty,
            medlwd = medlwd
        )
        # do this for first method of panel only
        if ( is_first ) {
            # mark zero line, significant in both metrics
            abline(
                h = 0,
                lty = 2,
                col = 'gray'
            )
            # add legend to top panel only
            legend(
                legend_pos,
                unlist( method_to_label ),
                text.col = method_cols,
                bty = 'n'
            )
        }
    }
    # bottom panel
    for ( method in methods ) {
        # get index to match colors up
        i <- which( methods == method )
        # first method has some special cases
        is_first <- i == 1
        # plot data!
        boxplot(
            data_auc[[ method ]],
            add = !is_first, # add all but first time
            names = if ( is_first ) rs else rs_blank,
            xlab = '',
            ylab = if ( is_first ) lab_auc else '',
            ylim = range_auc,
            border = method_cols_alpha[ i ],
            col = NA,
            outline = outline,
            range = range,
            whisklty = whisklty,
            medlwd = medlwd
        )
        # do this for first method of panel only
        if ( is_first ) {
            # mark zero line, significant in both metrics
            abline(
                h = 0,
                lty = 2,
                col = 'gray'
            )
        }
    }
    # add outer margin
    mtext(
        "Number of PCs (r)",
        side = 1,
        line = 0.5,
        adj = 0.55,
        outer = TRUE
    )
    fig_end()
}

# actually make plot
# for main figure output (goes on papers)
boxplots_rmsd_auc(
    name_out = paste0( name_base, 'rmsd-auc' ),
    data_rmsd = data_rmsd,
    data_auc = data_auc,
    r_max = max( pcs ),
    legend_pos = legend_pos
)

######################
### LAMBDA vs RMSD ###
######################

# kinda boring in this case, only make in main output
if ( !pca_test ) {
    # this was surprisingly cleaner than expected
    fig_start(
        paste0( name_base, 'rmsd-vs-lambda' )
    )
    plot(
        tib$rmsd,
        tib$lambda,
        xlab = lab_rmsd,
        ylab = lab_lambda,
        log = 'y',
        pch = '.'
    )
    fig_end()
}
