shuffle <- function(
    cds_ref,
    cds_qry
) {

    list <- list(cds_ref, cds_qry)

    for(cds in list) {

        cell_labels <- colData(cds)$cell_type_sub

        cell_labels[cell_labels == "enveloping layer (EVL)"] <- "enveloping layer (EVL)/periderm 1"
        cell_labels[cell_labels == "periderm 1"] <- "enveloping layer (EVL)/periderm 1"

        cell_labels[cell_labels == "hatching gland 1"] <- "hatching gland 1/2/periderm 5"
        cell_labels[cell_labels == "hatching gland 2"] <- "hatching gland 1/2/periderm 5"
        cell_labels[cell_labels == "periderm 5"] <- "hatching gland 1/2/periderm 5"

        cell_labels[cell_labels == "intestinal (mid)"] <- "intestinal (mid)/macrophage"
        cell_labels[cell_labels == "macrophage"] <- "intestinal (mid)/macrophage"

        cell_labels[cell_labels == "slow-committed myocyte" ] <- "slow-committed myocyte/mature slow muscle 1"
        cell_labels[cell_labels == "mature slow muscle 1"] <- "slow-committed myocyte/mature slow muscle 1"

        cell_labels[cell_labels == "satellite cells"] <- "satellite cells/mature fast muscle 1"
        cell_labels[cell_labels == "mature fast muscle 1"] <- "satellite cells/mature fast muscle 1"

        colData(cds)$cell_type_sub <- cell_labels
    }
    
    return(list)
}


shuffle_all <- function(
    cds_ref,
    cds_qry
) {

    set.seed(001)
    list <- list(cds_ref, cds_qry)

    for(i in 1:2) {
        cds <- list[[i]]
        for(col in c("L1", "L2", "L3", "L4")) {
            colData(cds)[, col] <- sample(colData(cds)[, col])
        }
        list[[i]] <- cds
    }

    return(list)
}