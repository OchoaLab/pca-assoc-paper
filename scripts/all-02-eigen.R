# an attempt to talk about low/high dimensionality using eigenvalues, the old school way

library(popkin)
library(genio)
library(readr)
library(ochoalabtools)

# some other relevant file paths
smartpca_par <- '../../scripts/eigensoft-smartpca-par.txt' # relative to each data/$name/
file_popkin_tw <- 'popkin-tw.txt'
file_popkin_eval <- 'popkin.eval'
# pure eigensoft versions
file_eigensoft_tw <- 'eigensoft-tw.txt'
file_eigensoft_eval <- 'eigensoft.eval'
file_eigensoft_evec <- 'eigensoft.evec' # temporarily created
file_eigensoft_log <- 'eigensoft.log' # temporarily created
# p-value threshold (most things are either absurdly significant or very close to 1, so hopefully this doesn't matter that much)
tw_cut <- 0.01
twstats_tab <- '~/bin/EIG-7.2.1/POPGEN/twtable'

# uses globals twstats_tab, tw_cut
calc_tw <- function( file_eval, file_tw ) {
    # run twstats (from eigensoft)
    # don't overwrite for rsync sake
    if ( !file.exists( file_tw ) ) {
        system2(
            'twstats', # from path
            args = c(
                '-t', twstats_tab,
                '-i', file_eval,
                '-o', file_tw
            )
        )
    }

    # parse file and determine rank of matrices
    tw <- read_table( file_tw, col_types = 'iddddd' )
    # estimated rank!  Store for fig
    tw <- min( which( tw$`p-value` >= tw_cut ) ) - 1
    # NOTE: had done this (count cases below threshold) but p-values are not monotonic!
    # (the start getting bigger, then in the middle the get smaller again, and sometimes they're significant at the ends!!!)
    #sum( tw$`p-value` < tw_cut, na.rm = TRUE )

    # cleanup
    #unlink( file_eval ) # decided to keep these all!

    return( tw )
}

# go where the data is
setwd( '../data/' )

# read datasets info (names for inputs and output, colors, line types)
datasets <- read_tsv( 'datasets.txt', col_types = 'cccii' )
# add "density" for barplots (dashed lines mimic)
datasets$density <- ifelse( datasets$lty == 1, -1, 20 )

# load all kinship matrices to a list
#kinship <- list() 
evs_popkin <- list()
evs_eigensoft <- list()
tw_rank_popkin <- list()
tw_rank_eigensoft <- list()
for ( i in 1 : nrow( datasets ) ) {
    # gets reused often
    name_paper <- datasets$name_paper[ i ]
    # go where the data is
    setwd( datasets$name_dir[ i ] )

    # variable for this dataset (might get local edit)
    smartpca_par_local <- smartpca_par
    # simulations (non-"Real") have data inside reps
    if ( datasets$type[ i ] != 'Real' ) {
        setwd( 'rep-1' )
        smartpca_par_local <- paste0( '../', smartpca_par )
    }

    ### PURE EIGENSOFT
    
    # do it all the eigensoft way
    # this is slow, avoid repeating!
    if ( !file.exists( file_eigensoft_eval ) ) {
        system2(
            'smartpca', # from path
            args = c( '-p', smartpca_par_local ),
            stdout = 'eigensoft.log'
        )
        # remove other temporary files
        unlink( file_eigensoft_evec )
        unlink( file_eigensoft_log )
    }
    # load eigenvalues, store in list right away
    evs_eigensoft[[ name_paper ]] <- drop( read_matrix( file_eigensoft_eval ) )
    # calculate and save TW stats
    tw_rank_eigensoft[[ name_paper ]] <- calc_tw( file_eigensoft_eval, file_eigensoft_tw )

    ### POPKIN
    
    if ( !file.exists( file_popkin_eval ) ) {
        # load precomputed kinship matrix
        kinship_i <- read_grm( 'popkin', verbose = FALSE )$kinship
        # calculate eigenvalues as-is (no normalizations, coancestry transform or negative value caps)
        evs_popkin_i <- eigen( kinship_i )$values
        # save eigenvectors to a file, to pass to eigensoft
        write_matrix( file_popkin_eval, evs_popkin_i, ext = NA, verbose = FALSE )
    } else {
        # load precomputed data
        evs_popkin_i <- drop( read_matrix( file_popkin_eval ) )
    }
    
    # save for further analysis
    #kinship[[ name_paper ]] <- kinship_i
    evs_popkin[[ name_paper ]] <- evs_popkin_i

    # calculate and save TW stats too
    tw_rank_popkin[[ name_paper ]] <- calc_tw( file_popkin_eval, file_popkin_tw )

    # go back to base of data
    setwd( '..' )
    if ( datasets$type[ i ] != 'Real' )
        setwd( '..' )
}

# join into matrix
tw_rank <- rbind( unlist( tw_rank_eigensoft ), unlist( tw_rank_popkin ) )
# add row names for plots
rownames( tw_rank ) <- c('Standard', 'Popkin')

# report numbers too
print( tw_rank )
##          Large sample size sim. Small sample size sim. Family structure sim.
## Standard                      4                      3                    35
## Popkin                        5                      4                    21
##          Human Origins HGDP 1000 Genomes Human Origins sim. HGDP sim.
## Standard           574  150           29                534       151
## Popkin             106   46           67                101        30
##          1000 Genomes sim.
## Standard                12
## Popkin                  19

# normalized eigenvalues (variance explained)
varexp_popkin <- lapply( evs_popkin, function( x ) x / sum( x ) )
varexp_eigensoft <- lapply( evs_eigensoft, function( x ) x / sum( x ) )
# and cumulative version
varexpcum_popkin <- lapply( varexp_popkin, cumsum )
varexpcum_eigensoft <- lapply( varexp_eigensoft, cumsum )

## # what is the variance explained of the first 10 or 90 PCs?
## # failed experiment, no signal here either
## varexp_popkin90 <- sapply( varexp_popkin, function ( x ) sum( x[ 1:90 ] ) )
## varexp_popkin10 <- sapply( varexp_popkin, function ( x ) sum( x[ 1:10 ] ) )

#######################
### STANDALONE FIGS ###
#######################


# top PCs table across datasets, for simple barplots
eigenvalues_mat <- function( ev, r = 10 ) {
    sapply( ev, function(x) x[ 1:r ] )
}

# uses global `datasets`
plot_varexp10 <- function( varexp, r = 10, cex_leg = 1 ) {
    # make sure data is aligned!
    stopifnot( datasets$name_paper == names( varexp ) )
    # plot
    barplot(
        t( eigenvalues_mat( varexp, r = r ) ),
        beside = TRUE,
        legend.text = TRUE,
        xlab = 'Eigenvalue Rank',
        ylab = 'Variance Explained',
        names.arg = 1 : r,
        col = datasets$col,
        density = datasets$density,
        args.legend = list(
            cex = cex_leg,
            bty = 'n'
        )
    )
}

# uses global `datasets`
plot_varexpcum <- function( varexpcum, cex_leg = 1 ) {
    # make sure data is aligned!
    stopifnot( datasets$name_paper == names( varexpcum ) )
    # plot!
    plot(
        NA,
        xlim = c(0, 1),
        ylim = c(0, 1),
        xlab = 'Eigenvalue Rank Fraction',
        ylab = 'Cumulative Variance Explained'
    )
    for ( i in 1 : length( varexpcum ) ) {
        y <- varexpcum[[ i ]]
        x <- ( 1 : length( y ) ) / length( y )
        lines( x, y, col = datasets$col[i], lty = datasets$lty[i] )
    }
    legend(
        'bottomright',
        names( varexpcum ),
        col = datasets$col,
        lty = datasets$lty,
        cex = cex_leg,
        bty = 'n'
    )
}

# so far only eigensoft versions are standalone, but the other ones could be too!

# same for all legends
cex_legs <- 0.8
cex_legs2 <- 0.9 # actually last case looks better bigger

# get half width
width <- fig_width( )
# create square figure
fig_start(
    'eigensoft_varexpcum',
    width = width/2,
    height = width/2
)
plot_varexpcum( varexpcum_eigensoft, cex_leg = cex_legs * 3/4 )
fig_end()

# create wide figure
fig_start(
    'eigensoft_varexp10',
    width = width,
    height = width/2
)
plot_varexp10( varexp_eigensoft, cex_leg = cex_legs2 * 3/4 )
fig_end()



################
### MAIN FIG ###
################

# lower margin is bigger for first panel only
mar1 <- c(7, 3, 2, 0) + 0.2
mar2 <- c(3, 3, 2, 0) + 0.2

# truncating extreme rank estimates
rank_cut <- 160

# get max width
width <- fig_width()
# create square figure
fig_start(
    'eigen',
    width = width,
    height = width
)
# add a custom layout
layout(
    rbind(
        c(1, 2),
        c(3, 3)
    )
)

### TWSTATS panel
par( mar = mar1 )
barplot(
    unlist( tw_rank_popkin ),
    ylab = 'Kinship Rank Est.',
    las = 2,
    cex.names = cex_legs, # reduce bar labels too
    col = datasets$col,
    density = datasets$density
)
panel_letter('A', adj = -0.15 )

### VAR EXP total cumulative panel
par( mar = mar2 )
plot_varexpcum( varexpcum_popkin, cex_leg = cex_legs )
panel_letter('B', adj = -0.15 )

### VAR EXP top 10 panel
plot_varexp10( varexp_popkin, cex_leg = cex_legs2 )
panel_letter('C', adj = -0.07 )

fig_end()


# other similar plots that were not intesting and didn't make it to paper

#barplot( t( eigenvalues_mat( evs_popkin ) ), beside = TRUE, legend.text = TRUE, xlab = 'Eigenvalue rank', ylab = 'Eigenvalue', names.arg = 1:10 )
## barplot( t( eigenvalues_mat( varexp_popkin ) ), beside = TRUE, legend.text = TRUE, log = 'y' )
## barplot( t( eigenvalues_mat( varexpcum_popkin ) ), beside = TRUE, legend.text = TRUE )





## # apples-to-apples comparison between real datasets and their tree counterparts
## plot_eigenvals_two <- function( i, j, log = 'xy' ) {
##     plot(
##         evs_popkin[[i]],
##         evs_popkin[[j]],
##         xlab = datasets$name_paper[ i ],
##         ylab = datasets$name_paper[ j ],
##         log = log
##     )
##     abline(0,1, lty = 2, col = 'gray')
## }

## plot_eigenvals_two( 1, 3 )
## plot_eigenvals_two( 4, 7 )
## plot_eigenvals_two( 5, 8 )
## plot_eigenvals_two( 6, 9 )




## # compare to eigensoft
## # only did for first simulation
## # conclusion: top eigenvalues disagree because eigensoft removes some individuals treated as outliers, or maybe it's also the standard kinship bias

## setwd( datasets$name_dir[1] )
## setwd( 'rep-1' )
## eigensoft_evs <- as.numeric( read_lines( 'eigensoft.eval' ) )
## setwd( '../..' )
## # normalize
## eigensoft_evs <- eigensoft_evs / sum( eigensoft_evs )
## plot( evs[[1]][ 1 : length(eigensoft_evs) ], eigensoft_evs, log = 'xy' )
## abline(0,1, lty = 2, col = 'gray')


## # processing for a single kinship matrix
## # normalizes to variance explained, 
## eigenvalues <- function( kinship, coanc = FALSE, nonneg = FALSE ) {
##     # this may be better for accentuating dimensionality, not sure if it's appropriate for family case though as it ceases to be positive-definite
##     if ( coanc )
##         kinship <- inbr_diag( kinship )
##     # eigenvalues
##     ev <- eigen( kinship )$values
##     # treat negatives as zeroes
##     if ( nonneg )
##         ev[ ev < 0 ] <- 0
##     # turn to variance explained
##     ev <- ev / sum( ev )
##     return( ev )
## }

## # calculate eigenvalues for entire list
## evs <- lapply( kinship, eigenvalues ) # overwrites unnormalized eigenvectors, meh
## evc <- lapply( kinship, eigenvalues, coanc = TRUE )
## evsn <- lapply( kinship, eigenvalues, nonneg = TRUE )
## evcn <- lapply( kinship, eigenvalues, nonneg = TRUE, coanc = TRUE )

## # shows negatives are negligible here
## sapply( evs, function(x) sum(x[x<0]) )
## ## Large sample size sim. Small sample size sim.  Family structure sim. 
## ##           0.000000e+00           0.000000e+00           0.000000e+00 
## ##          Human Origins                   HGDP           1000 Genomes 
## ##           0.000000e+00           0.000000e+00          -1.117561e-03 
## ##     Human Origins sim.              HGDP sim.      1000 Genomes sim. 
## ##           0.000000e+00           0.000000e+00          -8.068318e-05

## # but they are non-negligible here!
## # way worst for family structure
## sapply( evc, function(x) sum(x[x<0]) )
## ## Large sample size sim. Small sample size sim.  Family structure sim. 
## ##            -0.19506873            -0.05797959            -1.15818374 
## ##          Human Origins                   HGDP           1000 Genomes 
## ##            -0.10670092            -0.29532426            -0.17944568 
## ##     Human Origins sim.              HGDP sim.      1000 Genomes sim. 
## ##            -0.12347837            -0.09465313            -0.17247519 

## # these show evs and evsn are practically identical
## # so we will stop considering it
## plot( unlist(evs), unlist(evsn) )
## abline(0,1, lty = 2, col = 'gray')
## plot( unlist(evs), unlist(evsn), log = 'xy' )
## abline(0,1, lty = 2, col = 'gray')

## # these show there is less agreement, but by far the largest differences are for the family simulation
## plot( unlist(evc), unlist(evcn) )
## abline(0,1, lty = 2, col = 'gray')
## plot( unlist(evc), unlist(evcn), log = 'xy' )
## abline(0,1, lty = 2, col = 'gray')

## # kinship had only two negatives, and only one of them is fairly big (TGP-real)
## sum( sapply( evs, min ) < 0 )
## #[1] 2
## min( sapply( evs, min ) )
## #[1] -0.001117561

## # coancestry all were negative, TGP-real was again most extreme
## sum( sapply( evc, min ) < 0 )
## #[1] 9
## min( sapply( evc, min ) )
## #[1] -0.00468201

## # validate non-neg cases
## sum( sapply( evcn, min ) < 0 ) # 0

## # plot comparing across datasets
## barplot( t( eigenvalues_mat( evs ) ), beside = TRUE, legend.text=TRUE )
## barplot( t( eigenvalues_mat( evs ) ), beside = TRUE, legend.text=TRUE, log = 'y' )
## barplot( t( eigenvalues_mat( evc ) ), beside = TRUE, legend.text=TRUE )
## barplot( t( eigenvalues_mat( evcn ) ), beside = TRUE, legend.text=TRUE )

## # coancestry seems harder to motivate and have big negative problems, let's ignore from here on

