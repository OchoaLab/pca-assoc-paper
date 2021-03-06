library(optparse)
library(popkin) # avg_kinship_subpops, etc
library(genio)
library(readr)
library(ochoalabtools)
library(tibble)

# calculates coancestry of subpopulations, and FST

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
subpop_info <- read_tsv( 'pops-annot.txt', comment = '#', show_col_types = FALSE )

# map subpopulations using sub-subpopulations
fam$superpop <- subpop_info$superpop[ match( fam$fam, subpop_info$pop ) ]

# confirm matching orders here too
stopifnot( fam$id == colnames( kinship ) )
# reorder individuals so subpopulations come out in desired order:
indexes <- order( match( fam$fam, subpop_info$pop ) )
fam <- fam[ indexes, ]
kinship <- kinship[ indexes, indexes ]

# calculate weights that balance everything!
weights <- weights_subpops( fam$superpop, fam$fam )
fst_est <- fst( kinship, weights )
# save FST for table
write_lines( fst_est, 'fst.txt' )

# transform diagonal just before plotting
kinship <- inbr_diag( kinship )

##############
### SUBPOP ###
##############

kinship_pop <- avg_kinship_subpops( kinship, fam$fam, subpop_info$pop )

# save this matrix for later
write_grm( 'popkin_subpops', kinship_pop )

# means we can use subpop_info to annotate plot without changes (sorting or subsetting)
# extra work is only necessary for Human Origins
names_pop <- colnames( kinship_pop )
if ( length( names_pop ) != nrow( subpop_info ) || all( names_pop != subpop_info$pop ) ) {
    # in Human Origins the problem is there are extra rows in info table, so filter that file
    indexes <- match( names_pop, subpop_info$pop )
    # subset now
    subpop_info <- subpop_info[ indexes, ]
    # check again for agreement
    stopifnot( subpop_info$pop == names_pop )
}

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
