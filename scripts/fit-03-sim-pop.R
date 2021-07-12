# set up admixture simulation with tree
# the tree was estimated in a previous step, here we just set up the admixture matrix and save/copy the necessary data into a new directory for replicates to read from

# a real-data/tree version of sim-00-sim-pop.R

# creates:
# - bnpsd.RData
# - fst.txt
# - kinship_mean.txt

library(optparse)
library(bnpsd)
library(genio)
library(readr)
library(popkin)

# the name is for dir only, actual file is just "data"
name_in <- 'data'

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

# load FAM table, for subpopulation labels
fam <- read_fam( name_in )

# read annotations
subpop_info <- read_tsv('pops-annot.txt', comment = '#')

# map subpopulations using sub-subpopulations
fam$superpop <- subpop_info$superpop[ match( fam$fam, subpop_info$pop ) ]

# load tree calculated last time
load( 'tree.RData' ) # loads: tree, coanc_est, kinship_pop

# save things into a new destination
name <- paste0( name, '_sim' )

# move to where the data will be
# let's not overwrite things, under the assumption that the simulations are very expensive to redo
# if the output directory exists, assume all the files we want are there too.  Only a non-existent output directory will work
setwd( '..' )
if ( dir.exists( name ) )
    stop('Output exists, will not overwrite: ', name)
# else create directory and move into there
dir.create( name )
setwd( name )


# labels come from fam file
# subpops are desired order, match tree order
admix_proportions <- admix_prop_indep_subpops( fam$fam, subpops = tree$tip.label )
# copy names of individuals, to have correspondence to real data (though these are completely simulated)
rownames( admix_proportions ) <- fam$id

# edit tree a bit
# save under this name
tree_subpops <- tree
# also remove root edge (otherwise `bnpsd` complains about it, doesn't make sense to simulate using it)
tree_subpops$root.edge <- NULL

# recalculate tree coancestry after root edge is removed!
coanc_est <- coanc_tree( tree_subpops )

# calculate FST from tree, ignoring individuals/admixture (trivial here).
# do want to weigh superpops evenly
# easiest to map onto tree, where tip.label is subpops, to ensure agreement
tree_subpops$superpop <- subpop_info$superpop[ match( tree_subpops$tip.label, subpop_info$pop ) ]
# calculate weights that balance everything!
weights_tree <- weights_subpops( tree_subpops$superpop )
# use coanc_est to get inbreeding values out of
# NOTE: must pass as inbreeding (vector), rather than self-kinship (matrix)!
fst_tree <- fst( diag( coanc_est ), weights_tree )
# save FST for table
write_lines( fst_tree, 'fst.txt' )

# compared manually, agreed!
## # as a check, calculate FST the regular way, with individuals
## # coanc_est is ancestral coancestry according to tree (came from 'tree.RData')
## coanc_tree_all <- coanc_admix( admix_proportions, coanc_est )
## # these are weights for individuals
## weights_ind <- weights_subpops( fam$superpop, fam$fam )
## # NOTE: must pass as inbreeding (vector), rather than self-kinship (matrix)!
## fst_ind <- fst( diag( coanc_tree_all ), weights_ind )

# save bnpsd data to an RData file
# here add copy of real `fam` table, for checks
save( admix_proportions, tree_subpops, fam, file = 'bnpsd.RData' )

## # try things with this info
## # compare to popkin
## # read existing popkin data
## data <- read_grm( 'popkin' )
## kinship <- data$kinship
## 
## # subpopulations are scrambled, but there's mostly agreement considering that
## library(popkin)
## plot_popkin(
##     list(
##         inbr_diag( kinship ),
##         coanc_tree_all
##     )
## )

# lastly, also need mean kinship to undifferentiate real MAF
# NOTES:
# - must be "unweighted", because MAF weighed individuals uniformly
# - but at individual level (not subpops)
# - and it must be kinship, not coancestry
# - (these are probably all small differences, but might as well get it right)
kinship_mean <- mean( coanc_to_kinship( coanc_admix( admix_proportions, coanc_est ) ) )
write_lines( kinship_mean, 'kinship_mean.txt' )
