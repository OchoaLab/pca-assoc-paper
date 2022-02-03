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
    # special params if it's log-x (includes 'xy')
    is_logx <- grepl( 'x', log )
    
    # figure out shared ranges
    # x-axis is in linear scale, but y-axis in log potentially (don't include zero)
    # in linear y it's 1 that is boring... so take that out too
    x_max <- max( sapply( cumulatives, function( x ) max( x$x ) ) )
    y_min <- min( sapply( cumulatives, function( x ) min( x$y ) ) )
    y_max <- max( sapply( cumulatives, function( x ) max( x$y[ x$y < 1 ] ) ) )
    # if x-axis is in log scale, exclude zero from min
    x_min <- 0
    if ( is_logx )
        x_min <- 1e-4 # smaller values are observed but they are rare
        #x_min <- min( sapply( cumulatives, function( x ) min( x$x[ x$x > 0 ] ) ) )
    
    width <- fig_width()
    fig_start(
        paste0( 'king_log-', log ),
        width = width,
        height = width / 2,
        mar_t = 0.5,
        mar_r = if ( is_logx ) 0.5 else 0
    )

    # create blank plot
    plot(
        NA,
        log = log,
        xlim = c( x_min, x_max ),
        ylim = c( y_min, y_max ),
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
        if ( log == 'x' ) 'right' else if ( is_logx ) 'bottomleft' else 'topright',
        datasets$name_paper,
        col = datasets$col,
        lty = datasets$lty,
        cex = 0.8
    )

    # add reference line for 4th degree threshold
    x <- 2^(-11/2)
    abline( v = x, lty = 2, col = 'gray' )
    text( x = x, y = y_max, labels = '4th degree threshold', col = 'gray', adj = - 0.05 )
    
    fig_end()
}

# tried several versions
plot_king( cumulatives, log ='x' )
## plot_king( cumulatives, log ='y' )
## plot_king( cumulatives, log = 'xy' )
