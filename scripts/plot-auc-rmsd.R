# global constants

# plotting labels
lab_rmsd <- expression( bold( SRMSD[p] ) )
lab_auc <- expression( bold( AUC[PR] ) )

# transparency values for quartiles and extrema
alpha_q <- 0.4
alpha_e <- 0.2

# rule-of-thumb sRMSD cut
srmsd_cut <- 0.01

# reorganizes tibble into lists
# globals:
# - methods, pcs
tibble_to_lists_rmsd_auc <- function( tib ) {
    # initialize lists
    data_rmsd <- list()
    data_auc <- list()
    for ( method in methods ) {
        # sublists we want
        data_rmsd_i <- list()
        data_auc_i <- list()
        # subset big table
        tib_i <- tib[ tib$method == method, ]
        for ( pc in pcs ) {
            # subset table more
            tib_ij <- tib_i[ tib_i$pc == pc, ]
            # transfer vectors to lists
            # (shift index, since R doesn't like zero index)
            data_rmsd_i[[ pc + 1 ]] <- tib_ij$rmsd
            data_auc_i[[ pc + 1 ]] <- tib_ij$auc
        }
        # copy sublists to big lists
        data_rmsd[[ method ]] <- data_rmsd_i
        data_auc[[ method ]] <- data_auc_i
    }
    # return data in yet another list
    return(
        list(
            rmsd = data_rmsd,
            auc = data_auc
        )
    )
}

# make quantiles matrix
# expects single data type (either RMSD or AUC) for a single method (either PCA or LMM)
get_mean_quants_from_list <- function( data ) {
    # get number of columns essentially
    n <- length( data )
    # output matrix
    # 5 rows: standard quantiles (min, LQ, M, UQ, max)
    stats <- matrix( NA, nrow = 5, ncol = n )
    for ( i in 1 : n ) {
        # extract vector of data
        data_i <- data[[ i ]]
        # add stats to column i
        # NOTE: data can have NAs, let's always ignore those
        stats[ , i ] <- quantile( data_i, na.rm = TRUE )
    }
    # all done, return matrix!
    return( stats )
}

# assumes plot has been started, adds lines
plot_mean_quarts_lines <- function( stats, x, col, alpha_q, alpha_e ) {
    
    # plot median as simple line
    l <- 3 # median is 3rd stat row
    # there can be NAs when a PC/method is completely missing
    # this doesn't work with lines and polygons
    # identify by missing means only
    indexes_keep <- !is.na( stats[ l, ] )
    # subset entire data now, for simplicity
    x <- x[ indexes_keep ]
    stats <- stats[ , indexes_keep, drop = FALSE ]
    # actually plot lines now
    lines(
        x,
        stats[ l, ],
        lty = 1,
        col = col
    )
    
    # build polygon for quartiles
    # essentially, first walk through one line, then walk the other line backwards
    xq <- c( x, rev(x) )
    yq <- c( stats[ 2, ], rev( stats[ 4, ] ) )
    polygon(
        xq,
        yq,
        col = alpha( col, alpha_q ),
        border = FALSE # no border lines
    )

    # repeat with greater transparency for the extrema
    # reuse xq above
    yq <- c( stats[ 1, ], rev( stats[ 5, ] ) )
    polygon(
        xq,
        yq,
        col = alpha( col, alpha_e ),
        border = FALSE # no border lines
    )
}

# generic legend for quantiles data
leg_mean_quarts <- function( alpha_q, alpha_e, x = 'topright', col = 'black', cex = 1, title = 'Distribution' ) {
    # because the quartile area is always inside the extrema area, the combined alpha is this:
    # https://en.wikipedia.org/wiki/Alpha_compositing
    alpha_qe <- alpha_q + alpha_e * ( 1 - alpha_q )
    # not sure where to put it yet, will be data dependent
    # idea for pch for shaded area: (using `fill` didn't work!)
    # https://stat.ethz.ch/pipermail/r-help/2007-September/140415.html
    legend(
        x,
        c( 'median', 'quartiles', 'extrema' ),
        lty = c( 1, NA, NA ),
        pch = c( NA, 15, 15 ),
        pt.cex = 2 * cex,
        col = c(
            col,
            alpha( col, alpha_qe ),
            alpha( col, alpha_e )
        ),
        cex = cex,
        bty = 'n',
        title = title
    )
}

# globals:
# - method_to_label
# - method_cols
leg_methods <- function( inset_x = 0.15, cex = 1 ) {
    legend(
        'topright',
        unlist( method_to_label ),
        text.col = method_cols,
        lty = 1,
        col = method_cols,
        cex = cex,
        bty = 'n',
        title = 'Assoc. Model',
        title.col = 'black',
        inset = c( inset_x, 0 )
    )
}

# globals:
# - methods
# - method_cols
# - alpha_q, alpha_e
lineplots_rmsd_auc_one_panel <- function( data, lab, r_max, guide_max = FALSE, main = '' ) {
    # assumes PDF or whatever device has been created, and panels laid out
    # does not plot legends (that's external, for more control)
    
    # get common range of data plotted
    # always include zero in range for all cases (RMSD and AUC)
    ylim <- range( 0, unlist(data), na.rm = TRUE)

    # data locations
    x <- 0 : r_max
    # set a tighter xlim (default leaves lots of space on the sides)
    xlim <- range( x ) + 2 * c( 1, -1 )
    
    # start blank plot
    plot(
        NA,
        xlim = xlim,
        ylim = ylim,
        xlab = "", # will use outer labels here
        ylab = lab,
        main = main
    )
    
    # background stuff
    # mark zero line, significant in both metrics
    abline(
        h = 0,
        lty = 2,
        col = 'gray'
    )

    # method curves/areas go last
    # first pass gets data (in case we want a guide line in the background)
    stats <- list()
    for ( method in methods ) {
        # get quantile matrix for this method
        stats[[ method ]] <- get_mean_quants_from_list( data[[ method ]] )
    }
    
    # plot guideline if needed
    if ( guide_max ) {
        # gave it a short name, but this is the maximum of both medians
        # extract it from the data
        median_max <- 0 # always larger than this in the AUC case we're interested in
        for ( method in methods ) {
            # median is third row, get its max
            median_max <- max( stats[[ method ]][ 3, ], median_max )
        }
        # now plot line
        abline(
            h = median_max,
            lty = 2,
            col = 'gray'
        )
    } else {
        # assumptions: `guide_max = TRUE` only for AUC plots, so it's FALSE for RMSD plots (here)
        # add area to mark acceptably low SRMSDs
        # NOTE: `x` overflows from `c(0, 90)` (from -10 to 100) because we want area to extend like it's just in the background (looks better, contrasts with data)
        polygon(
            c(-10, 100, 100, -10),
            c( srmsd_cut, srmsd_cut, -srmsd_cut, -srmsd_cut ),
            col = alpha( 'gray', alpha_q ), # the darker of the two areas to see it more easily
            border = FALSE # no border lines
        )
    }
    
    # second pass actually plots
    for ( method in methods ) {
        # get index to match colors up
        i <- which( methods == method )

        # plot lines now
        plot_mean_quarts_lines( stats[[ method ]], x = x, col = method_cols[ i ], alpha_q, alpha_e )
    }
}

# global vars:
# - methods
# - method_to_label
# - method_cols
lineplots_rmsd_auc <- function(
                               name_out,
                               data_rmsd,
                               data_auc,
                               r_max,
                               legend_pos = 'topright'
                               ) {
    # start PDF
    fig_start(name_out, width = 7, height = 4, mar_b = 1.5)
    # add lower margin, so inner margins can be even smaller
    par( oma = c(1.5, 0, 0, 0) )
    # two panels
    par( mfrow = c(2, 1) )
    # change tick mar frequency on x axis
    par( lab = c(10, 3, 7) )

    # top panel: RMSD
    lineplots_rmsd_auc_one_panel( data_rmsd, lab_rmsd, r_max )
    
    # add all legends to top panel
    # in single-dataset version, legends just need to be a bit smaller and be shifted a bit
    leg_methods( inset_x = 0.2, cex = 0.8 )
    leg_mean_quarts( alpha_q, alpha_e, cex = 0.8 )
    
    # bottom panel: AUC
    lineplots_rmsd_auc_one_panel( data_auc, lab_auc, r_max, guide_max = TRUE )
    
    # add outer margin label
    mtext(
        "Number of PCs (r)",
        side = 1,
        line = 0.5,
        adj = 0.55,
        outer = TRUE
    )
    fig_end()
}

