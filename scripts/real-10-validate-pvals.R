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
    make_option("--herit", type = "double", default = 0.8, 
                help = "heritability", metavar = "double"),
    make_option("--m_causal_fac", type = "double", default = 10,
                help = "Proportion of individuals to causal loci", metavar = "double"),
    make_option("--fes", action = "store_true", default = FALSE, 
                help = "Use FES instead of RC trait model")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name <- opt$bfile
rep_max <- opt$rep
n_pcs_max <- opt$n_pcs
fes <- opt$fes
m_causal_fac <- opt$m_causal_fac
herit <- opt$herit

# stop if name is missing
if ( is.na(name) )
    stop('`--bfile` terminal option is required!')

# move to where the data is
setwd( '../data/' )
setwd( name )

# remember this location, using global path, for easily returning across multiple levels when done
dir_base <- getwd()

# in simulations, each rep has its own bim file, but they all have the same number of loci, so let's just read it from rep-1
if (opt$sim)
    name_in <- paste0( 'rep-1/', name_in )
# get number of loci
m_loci <- count_lines( name_in, 'bim' )

for ( rep in 1 : rep_max ) {
    # specify location of files to process, as many levels as needed
    dir_out <- paste0( 'rep-', rep, '/' )
    if ( fes )
        dir_out <- paste0( dir_out, 'fes/' )
    if ( m_causal_fac != 10 )
        dir_out <- paste0( dir_out, 'm_causal_fac-', m_causal_fac, '/' )
    if ( herit != 0.8 )
        dir_out <- paste0( dir_out, 'h', herit, '/' )
    
    # skip reps that we haven't calculated at all
    if ( !dir.exists( dir_out ) ) {
        if ( opt$final ) {
            # fatal if file is missing and we set `--final` option
            stop( 'whole rep ', rep, ' MISSING!' )
        } else 
            next
    }
    setwd( dir_out )
    
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
    
    # return to base when done with this rep
    setwd( dir_base )
}
