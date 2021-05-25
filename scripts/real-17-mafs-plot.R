# plots MAF distributions gathered in previous script

library(tibble)
library(ochoalabtools)

# names of datasets (paired list)
# data matches rmsd-vs-lambdas fig, including short names and colors
datasets <- tibble(
    name_short = c(
        'Large sample size sim.',
        'Small sample size sim.',
        'Family structure sim.',
        'Human Origins sim',
        'Human Origins',
        'HGDP sim.',
        'HGDP',
        '1000 Genomes sim.',
        '1000 Genomes'
    ),
    name_long = c(
        'sim-n1000-k10-f0.1-s0.5-g1',
        'sim-n100-k10-f0.1-s0.5-g1',
        'sim-n1000-k10-f0.1-s0.5-g20',
        'HoPacAll_ld_prune_1000kb_0.3_sim',
        'HoPacAll_ld_prune_1000kb_0.3',
        #'hgdp_wgs_autosomes_ld_prune_1000kb_0.3',
        'hgdp_wgs_autosomes_ld_prune_1000kb_0.3_maf-0.01_sim',
        'hgdp_wgs_autosomes_ld_prune_1000kb_0.3_maf-0.01',
        #'all_phase3_filt-minimal_ld_prune_1000kb_0.3_thinned-0.1'
        'all_phase3_filt-minimal_ld_prune_1000kb_0.3_maf-0.01_sim',
        'all_phase3_filt-minimal_ld_prune_1000kb_0.3_maf-0.01'
    ),
    col = c(1, 2, 3, 4, 4, 5, 5, 6, 6),
    lty = c(2, 2, 2, 2, 1, 2, 1, 2, 1)
)

# move to where the data is
setwd( '../data/' )

# loads individual `maf` vectors into named list for all datasets
mafs <- list()
for ( name in datasets$name_long ) {
    setwd( name )
    load( 'maf.RData' )
    mafs[[ name ]] <- maf
    setwd( '..' )
}

plot_mafs_panel <- function( datasets, mafs, leg_cex = 0.5, leg_seg_len = 3 ) {
    # set up plot area
    plot(
        NA,
        # axis labels for first time
        xlab = 'Minor Allele Frequency',
        ylab = 'Cumulative Probability',
        main = '',
        xlim = c(0, 0.5),
        ylim = c(0, 1)
    )

    # navigate datasets
    for ( i in 1 : nrow( datasets ) ) {
        # get name used in key
        name_long <- datasets$name_long[ i ]
        
        # a nicer alternative to `ecdf`, more flexibility for me
        # https://stat.ethz.ch/pipermail/r-help/2016-June/439416.html
        x <- sort( mafs[[ name_long ]] )
        y <- ppoints( x )
        
        # plot, with add as needed
        lines(
            x,
            y,
            col = datasets$col[ i ],
            lty = datasets$lty[ i ]
        )
    }

    # legend, (almost the) same as RMSD-vs-lambda fig
    # reverse order because it matches curve appearance better
    legend(
        'bottomright',
        # these two are reversed together
        rev( datasets$name_short ),
        text.col = rev( datasets$col ),
        title = 'Dataset',
        title.col = 'black',
        pch = NA,
        col = rev( datasets$col ),
        lty = rev( datasets$lty ),
        bty = 'n',
        cex = leg_cex,
        seg.len = leg_seg_len
    )
}

# make independent plot
fig_start( 'mafs' )
plot_mafs_panel( datasets, mafs )
fig_end()

# save data and panel function for replotting as part of bigger, multipanel figure
save( plot_mafs_panel, datasets, mafs, file = 'mafs.RData' )
