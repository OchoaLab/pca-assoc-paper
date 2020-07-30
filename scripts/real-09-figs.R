# this script reads big table of AUC and RMSD estimates

library(optparse)
library(scales) # for transparency
library(readr)
library(ochoalabtools)

# constants
# output file name (big table)
file_table <- 'sum.txt.gz'
# for main figure output (goes on papers)
name_out <- 'sum-rmsd-auc'

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--bfile", type = "character", default = NA, 
                help = "Base name for input plink files (.BED/BIM/FAM)", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name <- opt$bfile

# stop if name is missing
if ( is.na(name) )
    stop('`--bfile` terminal option is required!')

# move to where the data is
setwd( '../data/' )
setwd( name )

# read the big table!
tib <- read_tsv(
    file_table,
    col_types = 'ciiddd'
)

# extract methods from table itself
methods <- sort( unique(tib$method) )
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
    'sum-num-reps',
    width = 12
)
barplot(
    counts,
    beside = TRUE,
    xlab = '# PCs',
    ylab = 'replicates',
    legend.text = methods
)
fig_end()

################
### AUC/RMSD ###
################

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
lab_rmsd <- expression( bold( RMSD[p] ) )
lab_auc <- expression( bold( AUC[PR] ) )
lab_lambda <- expression( bold( paste("Inflation Factor (", lambda, ")") ) )

# fork of version Yiqi used, inputs are organized differently
boxplots_rmsd_auc <- function(
                              name_out,
                              data_rmsd,
                              data_auc,
                              legend_pos = 'topright',
                              name_fixed = 'pca-plink',
                              name_mixed = 'gcta',
                              col_fixed = 'red',
                              col_mixed = 'blue',
                              alpha = 0.5
                              ) {
    # extract fixed effects series
    rmsd_fixed <-  data_rmsd[[ name_fixed ]]
    auc_fixed <-  data_auc[[ name_fixed ]]
    # same for mixed effects
    rmsd_mixed <-  data_rmsd[[ name_mixed ]]
    auc_mixed <-  data_auc[[ name_mixed ]]
    
    # get common range of data plotted
    # always include zero in range for both
    range_rmsd <- range( 0, unlist(rmsd_fixed), unlist(rmsd_mixed), na.rm = TRUE)
    range_auc <- range( 0, unlist(auc_fixed), unlist(auc_mixed), na.rm = TRUE)
    r_max <- max(
        length( rmsd_fixed ), 
        length( auc_fixed ), 
        length( rmsd_mixed ), 
        length( auc_mixed )
    ) - 1 # turn indexes back to 0-based (first PC is 0)
    
    # make labels
    # blank all non-multiples of 10, for plot
    rs <- 0 : r_max
    rs[ rs %% 10 != 0 ] <- NA
    # overlaid panels have all blanks (R doesn't handle this correctly otherwise)
    rs_blank <- rep.int( NA, r_max+1)
    
    # add transarency to colors!!!
    col_fixed_alpha <- alpha(col_fixed, alpha)
    col_mixed_alpha <- alpha(col_mixed, alpha)
    
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
    boxplot(
        rmsd_fixed,
        names = rs,
        xlab = "",
        ylab = lab_rmsd,
        ylim = range_rmsd,
        border = col_fixed_alpha,
        col = NA,
        outline = outline,
        range = range,
        whisklty = whisklty,
        medlwd = medlwd
    )
    # mark zero line, significant in both metrics
    abline(
        h = 0,
        lty = 2,
        col = 'gray'
    )
    boxplot(
        rmsd_mixed,
        names = rs_blank,
        ylim = range_rmsd,
        add = TRUE,
        border = col_mixed_alpha,
        col = NA,
        outline = outline,
        range = range,
        whisklty = whisklty,
        medlwd = medlwd
    )
    # add legend to top panel only
    legend(
        legend_pos,
        c('Fixed effects (PCA)', 'Mixed effects (LMM+PCA)'),
        text.col = c( col_fixed, col_mixed ),
        bty = 'n'
    )
    # bottom panel
    boxplot(
        auc_fixed,
        names = rs,
        xlab = '',
        ylab = lab_auc,
        ylim = range_auc,
        border = col_fixed_alpha,
        col = NA,
        outline = outline,
        range = range,
        whisklty = whisklty,
        medlwd = medlwd
    )
    # mark zero line, significant in both metrics
    abline(
        h = 0,
        lty = 2,
        col = 'gray'
    )
    boxplot(
        auc_mixed,
        names = rs_blank,
        ylim = range_auc,
        add = TRUE,
        border = col_mixed_alpha,
        col = NA,
        outline = outline,
        range = range,
        whisklty = whisklty,
        medlwd = medlwd
    )
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
boxplots_rmsd_auc(
    name_out,
    data_rmsd,
    data_auc
)

######################
### LAMBDA vs RMSD ###
######################

# this was surprisingly cleaner than expected
fig_start(
    'sum-rmsd-vs-lambda'
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
