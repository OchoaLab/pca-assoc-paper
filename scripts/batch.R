library(ochoalabtools)

# slurm job submission script for DCC

# using biostats server
#partition <- 'biostat'
# using ochoalab server
#partition <- 'ochoalab'
# draw partitions probabilistically!
partitions <- c( 'ochoalab', 'biostat' )
# note this works without normalizing!
partitions_probs <- c(2, 5)

# shared items for all runs
herit_low <- FALSE
env <- TRUE
king_cutoff <- TRUE
#bfile <- 'sim-n100-k10-f0.1-s0.5-g1'; short_orig <- 's'
#bfile <- 'sim-n1000-k10-f0.1-s0.5-g1'; short_orig <- 'l'
#bfile <- 'sim-n1000-k10-f0.1-s0.5-g20'; short_orig <- 'f'
#bfile <- 'HoPacAll_ld_prune_1000kb_0.3_maf-0.01'; short_orig <- 'h'
#bfile <- 'hgdp_wgs_autosomes_ld_prune_1000kb_0.3_maf-0.01_geno-0.1'; short_orig <- 'd'
bfile <- 'tgp-nygc-autosomes_ld_prune_1000kb_0.3_maf-0.01'; short_orig <- 'k'
#bfile <- 'HoPacAll_ld_prune_1000kb_0.3_maf-0.01_sim'; short_orig <- 'H'
#bfile <- 'hgdp_wgs_autosomes_ld_prune_1000kb_0.3_maf-0.01_geno-0.1_sim'; short_orig <- 'D'
#bfile <- 'tgp-nygc-autosomes_ld_prune_1000kb_0.3_maf-0.01_sim'; short_orig <- 'K'

# needed to make sure each process gets enough memory if all my jobs saturate the machines
# (didn't have to change last time, but before TGP was thinned I needed 4 threads for GCTA runs)
threads <- 1L
# setup loop
reps_max <- 50L

plink_cases <- c(TRUE, FALSE)
fes_cases <- c(FALSE, TRUE)

if ( king_cutoff ) {
    # automatically fix dataset name for this case
    bfile <- paste0( bfile, '_king-cutoff-4' )
    # and only run FES
    fes_cases <- TRUE
}

#############################

# main submission steps
# uses array jobs, arraying over PCs
submit_rep_pcs <- function( rep, bfile, threads, fes, herit_low, env, plink, partition, labs = FALSE ) {
    # GCTA uses more memory, this works for the largest cases
    mem <- if ( plink ) '4G' else '16G'

    # select script automatically from boolean
    script <- if ( plink ) 'real-06-pca-plink.R' else 'real-05-gcta.R'

    # default number of PCs hardcoded here
    pcs <- '0-90'
    # only king_cutoff has obvious special cases (hack of "array" notation for a single or two values)
    if ( king_cutoff )
        pcs <- if ( plink ) 20 else '0,10'
    
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

    # a pair of parameters for low heritability simulations
    if ( herit_low )
        commands <- paste0( commands, ' --herit 0.3 --m_causal_fac 27' )
    
    # a pair of parameters for env simulations
    if ( env ) {
        commands <- paste0( commands, ' --env1 0.3 --env2 0.2' )
        # labs only makes sens for env and gcta
        if ( labs && !plink )
            commands <- paste0( commands, ' -l' )
    }
    
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
        account = partition, # account == partition in both cases that apply to me
        partition = partition,
        array = pcs
    )

    # submit job!
    batch_submit( name )

    # remove script when done (not needed after submission)
    batch_cleanup( name )
}


# automatically loop through main sims
for ( plink in plink_cases ) {
    # decide labs situation
    # most of the time there's no labs option
    labs_cases <- FALSE
    # labs needs GCTA and env (global)
    if ( !plink && env )
        labs_cases <- c(FALSE, TRUE)

    for ( labs in labs_cases ) {
        
        for ( fes in fes_cases ) {
            # set some settings/etc automatically fron here down

            ############
            ### NAME ###
            ############

            # reset `short` name after every loop, otherwise it gets messed up!
            short <- short_orig
            
            # set one of the levels of type of simulation in output name, so they don't overlap
            if ( env ) {
                # env only makes sense in combination with low heritability
                herit_low <- TRUE
                short <- paste0( short, 'e' ) # env (+ low-herit implied)
            } else if ( herit_low ) {
                short <- paste0( short, 'l' ) # low herit
            } else {
                short <- paste0( short, 'h' ) # high herit
            }
            
            # add trait type marker to output
            short <- paste0( short, if ( fes ) 'f' else 'r' )

            # so names don't overlap between plink and gcta runs
            short <- paste0( short, if ( plink ) 'p' else if ( labs ) 'l' else 'g' )

            ###########
            ### ETC ###
            ###########
            
            ## # I
            ## # test one rep
            ## rep <- 1L
            ## partition <- sample( partitions, 1L, prob = partitions_probs )
            ## submit_rep_pcs( rep, bfile, threads, fes, herit_low, env, plink, partition, labs )

            ## # II
            ## # finish rest of reps
            ## for ( rep in 2L : reps_max ) {
            ##     partition <- sample( partitions, 1L, prob = partitions_probs )
            ##     submit_rep_pcs( rep, bfile, threads, fes, herit_low, env, plink, partition, labs )
            ## }

            ########################

            # start loop
            for ( rep in 1L : reps_max ) {
                partition <- sample( partitions, 1L, prob = partitions_probs )
                submit_rep_pcs( rep, bfile, threads, fes, herit_low, env, plink, partition, labs )
            }
        }
    }
}
