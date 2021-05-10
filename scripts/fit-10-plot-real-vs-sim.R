library(optparse)
library(popkin)
library(genio)
library(readr)
library(ochoalabtools)
library(tibble)

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

#####################
### READ ALL DATA ###
#####################

# NOTE: for pipeline, it's most convenient to feed in _sim name
# here remove that suffix and read real data first
name_real <- sub( '_sim$', '', name )

# move to where the REAL data is
setwd( '../data/' )
setwd( name_real )

# read full fam data
fam <- read_fam( name_in )

# read existing popkin data
kinship_real <- read_grm( 'popkin' )$kinship

# read annotations
subpop_info <- read_tsv('pops-annot.txt', comment = '#')

# now switch dirs to simulated data
setwd( '..' )
setwd( name )
# data is always in first rep
setwd( 'rep-1' )

# read existing popkin data from simulation
kinship_sim <- read_grm( 'popkin' )$kinship

###############
### REORDER ###
###############

# map subpopulations using sub-subpopulations
fam$superpop <- subpop_info$superpop[ match( fam$fam, subpop_info$pop ) ]

# confirm matching orders here too
stopifnot( fam$id == rownames( kinship_real ) )
stopifnot( fam$id == rownames( kinship_sim ) )
# reorder individuals so subpopulations come out in desired order:
indexes <- order( match( fam$fam, subpop_info$pop ) )
fam <- fam[ indexes, ]
kinship_real <- kinship_real[ indexes, indexes ]
kinship_sim <- kinship_sim[ indexes, indexes ]

# transform diagonal just before plotting
kinship_real <- inbr_diag( kinship_real )
kinship_sim <- inbr_diag( kinship_sim )

############
### PLOT ###
############

# good kinship estimate has zero min (no capping necessary), only look at upper range and trust diagonal
alpha <- 0.01
capHi <- quantile( diag(kinship_real), probs = 1 - alpha )
kinship_real[kinship_real > capHi] <- capHi
kinship_sim[kinship_sim > capHi] <- capHi # just apply same cap as for real estimate

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
    'popkin_real-vs-sim',
    width = width,
    height = height
)

# main plot with all labeling bells and whistles
plot_popkin(
    list( kinship_real, kinship_sim ),
    labs = cbind(fam$fam, fam$superpop),
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
