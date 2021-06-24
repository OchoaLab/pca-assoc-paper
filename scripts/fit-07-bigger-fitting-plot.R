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

# load tables
name_F <- paste0( 'F-fit_m-', m_loci )
name_Fd <- paste0( 'Fd-fit_m-', m_loci )
name_B <- paste0( 'B-fit_m-', m_loci )
name_FB <- paste0( 'F-fit_m-', m_loci, '_vs-B' )

dat_F <- read_tsv( paste0( name_F, '.txt' ) )
dat_Fd <- read_tsv( paste0( name_Fd, '.txt' ) )
dat_B <- read_tsv( paste0( name_B, '.txt' ) )


### COMBINED PLOT F vs Fd vs B ###

range_kin <- range( 0, dat_F$kin, dat_Fd$kin, dat_B$kin )
range_maf <- range( 0, dat_F$maf, dat_Fd$maf, dat_B$maf )

fig_start( name_FB )
plot(
    NA,
    xlab = 'kinship error',
    ylab = 'MAF error',
    xlim = range_kin,
    ylim = range_maf
)
lines( dat_F$kin, dat_F$maf, col = 'red' )
lines( dat_Fd$kin, dat_Fd$maf, col = 'green' )
lines( dat_B$kin, dat_B$maf, col = 'blue' )
abline( h = 0, lty = 2, col = 'gray' )
abline( v = 0, lty = 2, col = 'gray' )
legend(
    'topright',
    c('undiff_af fixed', 'undiff_af distr', 'rbeta'),
    col = c('red', 'green', 'blue'),
    pch = 1
)
fig_end()


### SINGLE DATASET ###

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
