library(ochoalabtools)

# slurm job submission script for DCC

# shared items for all runs
plink <- TRUE # FALSE
fes <- FALSE # TRUE
#bfile <- 'sim-n100-k10-f0.1-s0.5-g1'; short <- 's'
#bfile <- 'sim-n1000-k10-f0.1-s0.5-g1'; short <- 'l'
#bfile <- 'sim-n1000-k10-f0.1-s0.5-g20'; short <- 'f'
#bfile <- 'HoPacAll_ld_prune_1000kb_0.3'; short <- 'h'
#bfile <- 'hgdp_wgs_autosomes_ld_prune_1000kb_0.3_maf-0.01'; short <- 'd'
#bfile <- 'all_phase3_filt-minimal_ld_prune_1000kb_0.3_maf-0.01'; short <- 'k'
#bfile <- 'HoPacAll_ld_prune_1000kb_0.3_sim'; short <- 'H'
bfile <- 'hgdp_wgs_autosomes_ld_prune_1000kb_0.3_maf-0.01_sim'; short <- 'D'
#bfile <- 'all_phase3_filt-minimal_ld_prune_1000kb_0.3_maf-0.01_sim'; short <- 'K'

# needed to make sure each process gets enough memory if all my jobs saturate the machines
# (didn't have to change last time, but before TGP was thinned I needed 4 threads for GCTA runs)
threads <- 1

#############################

# set some settings/etc automatically fron here down

# add trait type marker to output
short <- paste0( short, if ( fes ) 'f' else 'r' )

# GCTA uses more memory, this works for the largest cases
mem <- if ( plink ) '4G' else '16G'

# select script automatically from boolean
script <- if ( plink ) 'real-06-pca-plink.R' else 'real-05-gcta.R'

# so names don't overlap between plink and gcta runs
short <- paste0( short, if (plink) 'p' else 'g' )

# main submission steps
# global vars: script, bfile, short, mem, plink
# new version uses array jobs, arraying over PCs
submit_rep_pcs <- function( rep ) {
    # pass params to command
    # first do core set of parameters, more to be added below as needed
    commands <- paste0(
        'time Rscript ', script,
        ' --bfile ', bfile,
        ' -r ', rep,
        ' --n_pcs ', '$SLURM_ARRAY_TASK_ID',
        ' -t ', threads,
        ' --dcc'
    )

    # use desired trait type from global var (indicates trait file locations)
    if ( fes )
        commands <- paste0( commands, ' --fes' )

    # infer if this is a simulation or not (indicates genotype file locations)
    if ( grepl( 'sim', bfile ) )
        commands <- paste0( commands, ' --sim' )
    
    # load plink module if needed
    if ( plink ) {
        commands <- c(
            'module load Plink/2.00a3LM',
            commands,
            'module unload Plink/2.00a3LM'
        )
    }
    
    # give name clear params too
    name <- paste0( short, rep, '-%a' )
    
    # create submission file
    batch_writer(
        commands,
        name,
        mem = mem,
   	threads = threads,
        array = '0-90' # number of PCs hardcoded here
    )

    # submit job!
    batch_submit( name )

    # remove script when done (not needed after submission)
    batch_cleanup( name )
}

# setup loop
reps_max <- 50

# I
# test one rep
rep <- 1
submit_rep_pcs( rep )

## # II
## # finish rest of reps
## for ( rep in 2 : reps_max ) {
##     submit_rep_pcs( rep )
## }

########################

## # start loop
## for ( rep in 1 : reps_max ) {
##     submit_rep_pcs( rep )
## }

