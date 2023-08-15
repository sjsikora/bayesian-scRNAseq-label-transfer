assign_layer_labels <- function(cds_ref_D, cds_qry_D, monoProb, basoProb) {

    ontogeny <- list(
        c("L1.1", "L2.1", "DC", "DC"),
        c("L1.1", "L2.2", "NK", "NK"),
        c("L1.1", "L2.2", "B", "B"),
        c("L1.2", "L2.3", "L3.1", "Neu"),
        c("L1.2", "L2.3", "L3.1", "Baso"),
        c("L1.2", "L2.3", "Mono", "Mono"),
        c("L1.2", "Ery", "Ery", "Ery"),
        c("Undifferentiated", "Undifferentiated", "Undifferentiated", "Undifferentiated"),
        c("Prog", "Prog", "Prog", "Prog"),
        c("T", "T", "T", "T"),
        c("L1.1", "L2.1", "Mono", "Mono"),
        c("L1.2", "Ery", "Ery", "Baso")
    )

    names(ontogeny) <- c("DC", "NK", "B", "Neu", "Baso", "Mono", "Ery", "Undifferentiated", "Prog", "T", "MonoBreak", "BasoBreak")

    set.seed(005)

    list_of_cds <- list(cds_ref_D, cds_qry_D)

    for(i in 1:2) {
        cds <- list_of_cds[[i]]

        cell_labels <- colData(cds)$Cell.type.annotation

        df_list <- vector("list", length(cell_labels))
        
        for(j in 1:length(cell_labels)) {

            cell_label <- cell_labels[j]

            if(cell_label == 'Mono') {
                if(runif(1) < monoProb) {
                    cell_label <- 'MonoBreak'
                }
            }

            if(cell_label == 'Baso') {
                if(runif(1) < basoProb) {
                    cell_label <- 'BasoBreak'
                }
            }

            df_list[[j]] <- ontogeny[[cell_label]]
            
        }

        df <- do.call(rbind, df_list)

        colData(cds)[["L1"]] <- NULL
        colData(cds)[["L2"]] <- NULL
        colData(cds)[["L3"]] <- NULL
        colData(cds)[["L4"]] <- NULL

        colData(cds)[["L1"]] <- df[, 1]
        colData(cds)[["L2"]] <- df[, 2]
        colData(cds)[["L3"]] <- df[, 3]
        colData(cds)[["L4"]] <- df[, 4]

        list_of_cds[[i]] <- cds
    }

    return(list_of_cds)

}

check_ontogeny <- function(
    cds_nn
) {

    ontogeny <- list(
        c("L1.1", "L2.1", "DC", "DC"),
        c("L1.1", "L2.2", "NK", "NK"),
        c("L1.1", "L2.2", "B", "B"),
        c("L1.2", "L2.3", "L3.1", "Neu"),
        c("L1.2", "L2.3", "L3.1", "Baso"),
        c("L1.2", "L2.3", "Mono", "Mono"),
        c("L1.2", "Ery", "Ery", "Ery"),
        c("Undifferentiated", "Undifferentiated", "Undifferentiated", "Undifferentiated"),
        c("Prog", "Prog", "Prog", "Prog"),
        c("T", "T", "T", "T")
    )

    ontogenyTrack <- rep("Linear", dim(cds_nn)[1])

    names(ontogeny) <- c("DC", "NK", "B", "Neu", "Baso", "Mono", "Ery", "Undifferentiated", "Prog", "T")

    number_of_rows <- dim(cds_nn)[1]

    number_of_mono_cells <- 0
    number_of_brokenmono_cells <- 0

    number_of_baso_cells <- 0
    number_of_brokenbaso_cells <- 0

    for(i in 1:number_of_rows) {

        strict_ontogeny <- ontogeny[[cds_nn[i, 4]]]
        assigned_ontogeny <- cds_nn[i ,]


        if(assigned_ontogeny[4] == 'Mono') number_of_mono_cells <- number_of_mono_cells + 1
        if(assigned_ontogeny[4] == 'Baso') number_of_baso_cells <- number_of_baso_cells + 1


        if (!all(strict_ontogeny == assigned_ontogeny)) {
            print(paste0('Cell ', i, ' has a mismatch'))
            print(paste0('Strict: ', strict_ontogeny))
            print(paste0('Assigned: ', assigned_ontogeny))

            if(assigned_ontogeny[4] == 'Mono') {
                number_of_brokenmono_cells <- number_of_brokenmono_cells + 1
                ontogenyTrack[i] <- 'MonoBreak'
            }
            if(assigned_ontogeny[4] == 'Baso') {
                number_of_brokenbaso_cells <- number_of_brokenbaso_cells + 1
                ontogenyTrack[i] <- 'BasoBreak'
            }
        }
    }

    print(paste0('Number of mono cells: ', number_of_mono_cells))
    print(paste0('Number of broken mono cells: ', number_of_brokenmono_cells))
    
    print(paste0('Number of baso cells: ', number_of_baso_cells))
    print(paste0('Number of broken baso cells: ', number_of_brokenbaso_cells))

    return(ontogenyTrack)
}