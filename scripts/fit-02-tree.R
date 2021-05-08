library(optparse)
library(popkin)
library(genio)
library(ochoalabtools)
library(bnpsd)
library(ape)

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

file_data <- 'tree.RData'

# recalculate the tree data the first time
if ( !file.exists( file_data ) ) {

    # load precalculated matrix
    obj <- read_grm( 'popkin_subpops' )
    kinship_pop <- obj$kinship

    # fit tree!
    # NOTE: now best-match order is automatically performed here!
    tree <- fit_tree( kinship_pop )

    # NOTE: edge reordering messes up additive edges!
    # unfortunately coanc_tree just uses the previous ones, not recalculating, so we must do that now
    tree <- tree_additive( tree, force = TRUE )
    # and get coancestry of tree
    coanc_est <- coanc_tree( tree )

    # the coancestry has some subpops reorered, fix that here
    # this definitely affects tree fitting too!
    # NOTE: here we stick with the tree ordering, so keep that and order the other matrix instead
    ## indexes <- match( rownames( kinship_pop ), rownames( coanc_est ) )
    ## coanc_est <- coanc_est[ indexes, indexes ]
    indexes <- match( rownames( coanc_est ), rownames( kinship_pop ) )
    kinship_pop <- kinship_pop[ indexes, indexes ]

    # save all the data, particularly tree but also reordered raw data
    save( tree, coanc_est, kinship_pop, file = file_data )

} else {
    load( file_data ) # for replotting
}


## # score and decide how to reorder things in the tree, since the previous code wasn't perfect!
## names_exp <- rownames( kinship_pop )
## score_tip_order <- function( tree, names_exp ) {
##     # get actual names
##     names_obs <- tree$tip.label
##     # correspondence
##     indexes_obs <- match( names_obs, names_exp )
##     # score is departure from perfect order
##     indexes_exp <- 1 : length( names_obs )
##     # store first as vector (each position has a score)
##     scores <- indexes_obs - indexes_exp
##     # as departures can have positive and negative signs, overall score should be average of absolute values
##     # (avg or sum, doesn't matter much)
##     score_all <- mean( abs( scores ) )
##     return( score_all )
## }

## # after reordering
## score_tip_order( tree, names_exp )
## # [1] 28.04115
## # does it help to call this function again?
## tree2 <- rotateConstr( tree, names_exp )
## score_tip_order( tree2, names_exp )
## # [1] 28.04115 # NOPE, this algorithm is now stuck



## plot( tree )
## # add axis and label
## axisPhylo( backward = FALSE )
## mtext( 'Coancestry', side = 1, line = 2 )

# for trees
fig_start(
    'tree',
    width = 9 * 3,
    height = 4 * 3,
    mar_t = 2
)
#par( cex = 0.3 )
plot_popkin(
    list(
        tree,
        kinship_pop,
        coanc_est
    ),
    names = TRUE,
    names_cex = 0.5
)
fig_end()

## # verify rss/error is as expected
## # NOTE: RSS is a sum, so it probably scales with k_subpops^2
## tree$rss
## # [1] 0.02228358 # TGP
## # [1] 0.1705238  # HGDP
## # [1] 11.27201   # HO

## err <- kinship_pop - coanc_est
## err[ lower.tri(err) ] <- NA
## sum( err^2, na.rm=TRUE )
## # [1] 0.02228358 # TDP
## # [1] 0.1705238  # HGDP
## # [1] 11.27201   # HO
## mean( err^2, na.rm=TRUE )
## # [1] 6.348597e-05 # TGP
## # [1] 0.0001148308 # HGDP
## # [1] 0.0003802204 # HO
