library(optparse)
library(popkin)
library(genio)
library(readr)
library(ochoalabtools)
library(tibble)

## # load new functions from external scripts
## dir_orig <- getwd()
## setwd("../../scripts") # scripts of main GAS project
## source('gas_lmm_gcta.R')
## source('paths.R')
## setwd( dir_orig ) # go back to where we were

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

# read full fam data
fam <- read_fam( name_in )

# read existing popkin data
data <- read_grm( 'popkin' )
kinship <- data$kinship

# read annotations
subpop_info <- read_tsv('pops-annot.txt', comment = '#')

# map subpopulations using sub-subpopulations
fam$superpop <- subpop_info$superpop[ match( fam$fam, subpop_info$pop ) ]

# confirm matching orders here too
stopifnot( fam$id == colnames( kinship ) )
# reorder individuals so subpopulations come out in desired order:
indexes <- order( match( fam$fam, subpop_info$pop ) )
fam <- fam[ indexes, ]
kinship <- kinship[ indexes, indexes ]

# transform diagonal just before plotting
kinship <- inbr_diag( kinship )

##############
### SUBPOP ###
##############

# code adapted from Storey Lab (by Ochoa)

avgSubpopsKinship <- function(kinship, ind2pop, pops) {
    # make smaller kinship matrix averaging over populations
    # Note :
    # - ind2pop maps individuals in the order of kinship onto populations
    # - pops is unique populations in the desired order!

    # first a convenient cleanup: remove values in pops that are missing in ind2pop (that way we don't have to do it outside, we ensure agreement more generally)
    pops <- pops[ pops %in% ind2pop ]
    # initiate output matrix
    K <- length( pops )
    kinship_pops <- matrix( 0, nrow = K, ncol = K )
    
    # start averaging
    for (i in 1:K) {
        is <- ind2pop == pops[ i ] # booleans for individuals that are of this population
        # note i=j case is included 
        for (j in 1:i) {
            js <- ind2pop == pops[ j ] # booleans for individuals that are of this population
            phiij <- mean( kinship[ is, js ] ) # the mean value we want
            kinship_pops[ i, j ] <- phiij # store both ways
            kinship_pops[ j, i ] <- phiij # store both ways
        }
    }

    # store names of populations on matrix
    colnames(kinship_pops) <- pops
    rownames(kinship_pops) <- pops

    kinship_pops # return
}

kinship_pop <- avgSubpopsKinship( kinship, fam$fam, subpop_info$pop )

# save this matrix for later
write_grm( 'popkin_subpops', kinship_pop )

# means we can use subpop_info to annotate plot without changes (sorting or subsetting)
stopifnot( colnames( kinship_pop ) == subpop_info$pop )

############
### PLOT ###
############

# good kinship estimate has zero min (no capping necessary), only look at upper range and trust diagonal
alpha <- 0.01
capHi <- quantile( diag(kinship), probs = 1 - alpha )
kinship[kinship > capHi] <- capHi

# other parameters
line <- 3 # for outer labels
mar <- line + 1
# additional padding for legend top/bottom
leg_mar <- 0 # 8

# figure out dimensions
widthMax <- 7.87402
width <- widthMax
height <- width / 1.15 / 2

# start plot!
fig_start(
    'popkin_subpops',
    width = width,
    height = height
)

# main plot with all labeling bells and whistles
plot_popkin(
    list( kinship, kinship_pop ),
    labs = list(
        cbind(fam$fam, fam$superpop),
        cbind(subpop_info$pop, subpop_info$superpop)
    ),
    labs_cex = c(0.7, 1),
    labs_las = c(2, 0),
    labs_line = c(0.5, line),
    labs_lwd = c(0.1, 0.5),
    labs_sep = c(FALSE, TRUE),
    labs_even = c(TRUE, FALSE),
    leg_width = 0.2, # 0.12,
    mar = mar,
    leg_mar = c(mar + leg_mar, 0, leg_mar, 3),
    oma = 1
)

fig_end()
