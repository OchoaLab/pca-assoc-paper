# set up admixture simulation with tree
# the tree was estimated in a previous step, here we just set up the admixture matrix and save/copy the necessary data into a new directory for replicates to read from

# a real-data/tree version of sim-00-sim-pop.R

# creates:
# - bnpsd.RData

library(optparse)
library(bnpsd)
library(genio)

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

# load tree calculated last time
load( 'tree.RData' ) # load `tree` and other things we don't use here

# labels come from fam file
# subpops are desired order, match tree order
admix_proportions <- admix_prop_indep_subpops( fam$fam, subpops = tree$tip.label )
# copy names of individuals, to have correspondence to real data (though these are completely simulated)
rownames( admix_proportions ) <- fam$id

## # try things with this info
## # coanc_est is ancestral coancestry according to tree (came from 'tree.RData')
## coanc_tree_all <- coanc_admix( admix_proportions, coanc_est )
## 
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

# edit tree a bit
# save under this name
tree_subpops <- tree
# also remove root edge (otherwise `bnpsd` complains about it, doesn't make sense to simulate using it)
tree_subpops$root.edge <- NULL

# save bnpsd data to an RData file
# here add copy of real `fam` table, for checks
save( admix_proportions, tree_subpops, fam, file = 'bnpsd.RData' )
