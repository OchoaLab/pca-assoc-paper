# makes sure all p-value vectors have the right length (number of loci in data)

library(optparse)
library(genio)

# constants
methods <- c('pca-plink-pure', 'gcta') # 'pca-plink', 
# the name is for dir only, actual file is just "data"
name_in <- 'data'

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--bfile", type = "character", default = NA, 
                help = "Base name for input plink files (.BED/BIM/FAM)", metavar = "character"),
    make_option("--n_pcs", type = "integer", default = 90,
                help = "Max number of PCs", metavar = "int"),
    make_option(c("-r", "--rep"), type = "integer", default = 50,
                help = "Max replicates", metavar = "int"),
    make_option("--sim", action = "store_true", default = FALSE, 
                help = "Genotypes are simulated (rather than real; alters location only)"),
    make_option("--final", action = "store_true", default = FALSE, 
                help = "Dies if files are missing (otherwise they are skipped silently)"),
    make_option("--const_herit_loci", action = "store_true", default = FALSE, 
                help = "Causal coefficients constructed to result in constant per-locus heritability (saved in diff path)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name <- opt$bfile
rep_max <- opt$rep
n_pcs_max <- opt$n_pcs
const_herit_loci <- opt$const_herit_loci

# stop if name is missing
if ( is.na(name) )
    stop('`--bfile` terminal option is required!')

# move to where the data is
setwd( '../data/' )
setwd( name )

# in simulations, each rep has its own bim file, but they all have the same number of loci, so let's just read it from rep-1
if (opt$sim)
    name_in <- paste0( 'rep-1/', name_in )
# get number of loci
m_loci <- count_lines( name_in, 'bim' )

for ( rep in 1 : rep_max ) {
    # move higher to the "reps" location
    # this is so GCTA's temporary files don't overwrite files from other parallel runs
    dir_out <- paste0( 'rep-', rep )
    # skip reps that we haven't calculated at all
    if ( !dir.exists( dir_out ) )
        next
    setwd( dir_out )
    
    # move in an additional level in this case
    if ( const_herit_loci ) {
        dir_phen <- 'const_herit_loci/'
        # directory can be missing, skip in that case
        if ( !dir.exists( dir_phen ) ) {
            # just move down from rep-* case
            setwd( '..' )
            # skip rest
            next
        }
        # else move in
        setwd( dir_phen )
    }
    
    # start a big loop
    for ( method in methods ) {
        for ( n_pcs in 0 : n_pcs_max ) {
            # file to read
            file_pvals <- paste0( 'pvals_', method, '_', n_pcs, '.txt.gz' )

            if ( !file.exists( file_pvals ) ) {
                if ( opt$final ) {
                    # fatal if file is missing and we set `--final` option
                    stop(
                        'rep: ', rep,
                        ', method: ', method,
                        ', pcs: ', n_pcs,
                        ', MISSING!'
                    )
                } else 
                    next
            }

            # count lines using linux system terminal (fastest)
            # run wc on terminal
            n_lines <- system(
                paste0(
                    'zcat ',
                    file_pvals,
                    ' | wc -l'
                ),
                intern = TRUE
            )
            # returns string, turn into numeric
            n_lines <- as.numeric( n_lines )

            # data is either complete (m_loci) or result was NULL (1 line), nothing else is acceptable
            if ( n_lines != m_loci && n_lines != 1 ) {
                # this is what we wanted to find!
                stop(
                    'rep: ', rep,
                    ', method: ', method,
                    ', pcs: ', n_pcs,
                    ', lines: ', n_lines,
                    ', expected: ', m_loci
                )
            }
        }
    }
    
    # move back down when done with this rep
    setwd( '..' )
    # move back an additional level in this case
    if ( const_herit_loci )
        setwd( '..' )
}
