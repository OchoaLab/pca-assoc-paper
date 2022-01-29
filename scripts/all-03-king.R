# direct estimates of local/family structure

#library(popkin)
library(genio)
library(readr)
library(ochoalabtools)

# go where the data is
setwd( '../data/' )

# read datasets info (names for inputs and output, colors, line types)
datasets <- read_tsv( 'datasets.txt', col_types = 'cccii' )

# load all kinship matrices to a list
kinships <- list() 
for ( i in 1 : nrow( datasets ) ) {
    # gets reused often
    name_paper <- datasets$name_paper[ i ]
    # go where the data is
    setwd( datasets$name_dir[ i ] )

    # simulations (non-"Real") have data inside reps
    if ( datasets$type[ i ] != 'Real' )
        setwd( 'rep-1' )

    # parse precalculated binary matrix
    data <- read_grm( 'data', ext = 'king', shape = 'strict' )
    #fam <- data$fam
    kinships[[ name_paper ]] <- data$kinship
    
    # go back to base of data
    setwd( '..' )
    if ( datasets$type[ i ] != 'Real' )
        setwd( '..' )
}

# calculate cumulative curves
kin_cum <- function( kinship ) {
    # zero negative values
    kinship[ kinship < 0 ] <- 0
    # get unique values
    kinship_vec <- kinship[ lower.tri( kinship ) ]
    # table gives us order and counts per case!
    # values are ordered increasing *numerically* (was worried they'd be treated as text)
    x <- table( kinship_vec )
    # make decreasing now
    x <- rev( x )
    # this is desired cumulative counts
    y <- cumsum( x )
    # normalize counts so they're fraction
    # note max(y) == length( kinship_vec ) as desired
    y <- y / max( y )
    # recover values from table, as numeric
    x <- as.numeric( names( x ) )
    # return paired data
    return( list( x = x, y = y ) )
}
cumulatives <- lapply( kinships, kin_cum )

plot_king <- function( cumulatives, log = 'y' ) {
    # figure out shared ranges
    # x-axis is in linear scale, but y-axis in log potentially (don't include zero)
    x_max <- max( sapply( cumulatives, function( x ) max( x$x ) ) )
    y_min <- min( sapply( cumulatives, function( x ) min( x$y ) ) )

    width <- fig_width()
    fig_start(
        'king',
        width = width,
        height = width / 2,
        mar_t = 0.5
    )

    # create blank plot
    plot(
        NA,
        log = log,
        xlim = c( 0, x_max ),
        ylim = c( y_min, 1 ),
        xlab = 'Kinship estimate (KING-robust)',
        ylab = 'Cumulative distribution'
    )

    # now add lines for each dataset
    for ( i in 1 : nrow( datasets ) ) {
        cum_i <- cumulatives[[ i ]]
        lines(
            cum_i$x, 
            cum_i$y,
            col = datasets$col[i],
            lty = datasets$lty[i]
        )
    }

    # add legend
    legend(
        'topright',
        datasets$name_paper,
        col = datasets$col,
        lty = datasets$lty,
        cex = 0.8
    )
    fig_end()
}

plot_king( cumulatives )
#plot_king( cumulatives, log = '' )

# follow up of strange samples
# HO has the largest values, some exceeding sib and parent relatedness
kinship_ho <- kinships[[ 'Human Origins' ]]
xy <- arrayInd( which.max( kinship_ho ), dim( kinship_ho ) )
kinship_ho[ xy[1], xy[2] ]
# [1] 0.3989481 # matches plot
ids <- colnames( kinship_ho )
ids[ as.numeric( xy ) ]
# [1] "NAD96" "NAD95"

# on terminal I find this annotation:
## grep NAD95 data.fam 
## grep NAD96 data.fam 
## Chipewyan	NAD95	0	0	2	1
## Chipewyan	NAD96	0	0	2	1
