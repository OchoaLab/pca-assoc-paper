# plotting code for F-fit tables

library(optparse)
library(readr)
library(ochoalabtools)

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--bfile", type = "character", default = NA, 
                help = "Base name for input plink files (.BED/BIM/FAM)", metavar = "character"),
    make_option(c("-m", "--m_loci"), type = "integer", default = 10000, 
                help = "number of loci", metavar = "int"),
    make_option("--maf_min", type = "double", default = 0.01, 
                help = "Minimum MAF, to simulate MAF-based ascertainment bias", metavar = "double")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name <- opt$bfile
m_loci <- opt$m_loci
maf_min <- opt$maf_min

# stop if name is missing
if ( is.na(name) )
    stop('`--bfile` terminal option is required!')

# directory above must already exist
setwd( '../data/' )
setwd( name )
setwd( 'fits' )

# define some paths (often shared input txt and output pdf)
name_F <- paste0( 'F-fit_m-', m_loci )
name_Fd <- paste0( 'Fd-fit_m-', m_loci )
name_Kd <- paste0( 'Kd-fit_m-', m_loci )
name_B <- paste0( 'B-fit_m-', m_loci )
name_all <- paste0( 'all-fits_m-', m_loci )

# load tables
file_F <- paste0( name_F, '.txt' )
file_Fd <- paste0( name_Fd, '.txt' )
file_Kd <- paste0( name_Kd, '.txt' )
file_B <- paste0( name_B, '.txt' )
# gracefully skip things that don't exist
dat_F <- if ( file.exists( file_F ) ) read_tsv( file_F ) else NULL
dat_Fd <- if ( file.exists( file_Fd ) ) read_tsv( file_Fd ) else NULL
dat_Kd <- if ( file.exists( file_Kd ) ) read_tsv( file_Kd ) else NULL
dat_B <- if ( file.exists( file_B ) ) read_tsv( file_B ) else NULL


### COMBINED PLOT F vs Fd vs B ###

# NOTE: NULL dat_* cases are ignored without errors!
range_kin <- range( 0, dat_F$kin, dat_Fd$kin, dat_Kd$kin, dat_B$kin )
range_maf <- range( 0, dat_F$maf, dat_Fd$maf, dat_Kd$maf, dat_B$maf )

fig_start( name_all )
plot(
    NA,
    xlab = 'kinship error',
    ylab = 'MAF error',
    xlim = range_kin,
    ylim = range_maf
)
# these also get skipped silently if dat_* are NULL
lines( dat_F$kin, dat_F$maf, col = 'red' )
lines( dat_Fd$kin, dat_Fd$maf, col = 'green' )
lines( dat_Kd$kin, dat_Kd$maf, col = 'orange' )
lines( dat_B$kin, dat_B$maf, col = 'blue' )
abline( h = 0, lty = 2, col = 'gray' )
abline( v = 0, lty = 2, col = 'gray' )
# leave them all there, whatevs
legend(
    'topright',
    c('undiff_af fixed', 'undiff_af distr', 'kin distr', 'rbeta'),
    col = c('red', 'green', 'orange', 'blue'),
    pch = 1,
    cex = 0.7
)
fig_end()


### SINGLE DATASET ###

# these we will skip completely if there isn't data to plot

if ( !is.null( dat_F ) ) {
    fig_start( name_F, width = 6 )
    par( mfrow = c(1, 2) )
    plot(
        dat_F$fst,
        dat_F$kin,
        ylim = range( 0, dat_F$kin ),
        xlab = 'FST for undiff_af',
        ylab = 'kinship error'
    )
    plot(
        dat_F$fst,
        dat_F$maf,
        ylim = range( 0, dat_F$maf ),
        xlab = 'FST for undiff_af',
        ylab = 'MAF error'
    )
    fig_end()
}

if ( !is.null( dat_Fd ) ) {
    fig_start( name_Fd, width = 6 )
    par( mfrow = c(1, 2) )
    plot(
        dat_Fd$fst,
        dat_Fd$kin,
        ylim = range( 0, dat_Fd$kin ),
        xlab = 'FST for undiff_af',
        ylab = 'kinship error'
    )
    plot(
        dat_Fd$fst,
        dat_Fd$maf,
        ylim = range( 0, dat_Fd$maf ),
        xlab = 'FST for undiff_af',
        ylab = 'MAF error'
    )
    fig_end()
}

if ( !is.null( dat_Kd ) ) {
    fig_start( name_Kd, width = 6 )
    par( mfrow = c(1, 2) )
    plot(
        dat_Kd$fst,
        dat_Kd$kin,
        ylim = range( 0, dat_Kd$kin ),
        xlab = 'kinship_mean for undiff_af',
        ylab = 'kinship error'
    )
    plot(
        dat_Kd$fst,
        dat_Kd$maf,
        ylim = range( 0, dat_Kd$maf ),
        xlab = 'kinship_mean for undiff_af',
        ylab = 'MAF error'
    )
    fig_end()
}

if ( !is.null( dat_B ) ) {
    fig_start( name_B, width = 6 )
    par( mfrow = c(1, 2) )
    plot(
        dat_B$alpha,
        dat_B$kin,
        ylim = range( 0, dat_B$kin ),
        xlab = 'alpha for rbeta',
        ylab = 'kinship error',
        log = 'x'
    )
    plot(
        dat_B$alpha,
        dat_B$maf,
        ylim = range( 0, dat_B$maf ),
        xlab = 'alpha for rbeta',
        ylab = 'MAF error',
        log = 'x'
    )
    fig_end()
}
