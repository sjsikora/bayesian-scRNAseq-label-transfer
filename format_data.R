
load_split_downloaded_data <- function() {

    expression_matrix <- Matrix::readMM('lineage_pap/downloaded_10x/GSM4185643_stateFate_inVivo_normed_counts.mtx')
    expression_matrix <- Matrix::t(expression_matrix)

    cell_metadata <- read.delim('lineage_pap/downloaded_10x/GSM4185643_stateFate_inVivo_metadata.txt')
    
    
    cds <- new_cell_data_set(expression_matrix,
                    cell_metadata = cell_metadata)

    gene_names <- readLines('lineage_pap/downloaded_10x/GSM4185643_stateFate_inVivo_gene_names.txt')
    rownames(cds) <- gene_names

    cell_barcodes <- readLines('lineage_pap/downloaded_10x/GSM4185643_stateFate_inVivo_cell_barcodes.txt')
    colnames(cds) <- cell_barcodes

    number_of_cells <- dim(colData(cds_qry))[1]
    rownames(colData(cds_qry)) <- paste0('cell_qry_', 1:number_of_cells)

    #-----Split-------
    
    librarys <- colData(cds)$Library
    unique_librarys <- unique(librarys)

    set.seed(001)
    qry_librarys <- sample(unique_librarys, 25)
    ref_librarys <- unique_librarys[!unique_librarys %in% qry_librarys]

    #Subset cds based on qry_libtarys

    #cds_qry_D makes up about %36 of the total data
    cds_qry_D <- cds[, !is.na(colData(cds)$Library) & colData(cds)$Library %in% qry_librarys]
    cds_ref_D <- cds[, !is.na(colData(cds)$Library) & colData(cds)$Library %in% ref_librarys]

    cds_ref_D <- estimate_size_factors(cds_ref_D)
    cds_qry_D <- estimate_size_factors(cds_qry_D)

    reducedDim(cds_ref_D) <- NULL
    reducedDim(cds_ref_D) <- NULL
    reducedDim(cds_ref_D) <- NULL

    reducedDim(cds_qry_D) <- NULL
    reducedDim(cds_qry_D) <- NULL
    reducedDim(cds_qry_D) <- NULL

    cds_ref_D <- preprocess_cds(cds_ref_D)
    cds_ref_D <- reduce_dimension(cds_ref_D, build_nn_index=TRUE)

    save_transform_models(cds_ref_D, 'cds_ref_R_D_models')

    cds_qry_D <- load_transform_models(cds_qry_D, 'cds_ref_R_D_models')

    cds_qry_D <- preprocess_transform(cds_qry_D)
    cds_qry_D <- reduce_dimension_transform(cds_qry_D)

    colData(cds_ref_D)[['data_set']] <- 'reference'
    colData(cds_qry_D)[['data_set']] <- 'query'

    return(list(cds_ref_D, cds_qry_D))
}

assign_layer_labels <- function(cds_ref_D, cds_qry_D, monoProb = 0.3, basoProb = 0.3) {

    ontology <- list(
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

    names(ontology) <- c("DC", "NK", "B", "Neu", "Baso", "Mono", "Ery", "Undifferentiated", "Prog", "T", "MonoBreak", "BasoBreak")

    set.seed(005)

    list_of_cds <- list(cds_ref_D, cds_qry_D)

    for(i in 1:2) {
        cds <- list_of_cds[[i]]

        cell_labels <- colData(cds)$Cell.type.annotation

        df_list <- vector("list", length(cell_labels))
        
        for(j in 1:length(cell_labels)) {

            cell_label <- cell_labels[j]

            if(cell_label == 'Mono') {
                #10% to switch
                if(runif(1) < monoProb) {
                    cell_label <- 'MonoBreak'
                    print("monoBreak")
                }
            }

            if(cell_label == 'Baso') {
                #5% to switch
                if(runif(1) < basoProb) {
                    cell_label <- 'BasoBreak'
                    print("BasoBreak")
                }
            }

            df_list[[j]] <- ontology[[cell_label]]
            
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




check_ontology <- function(
    cds_nn
) {

    ontology <- list(
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

    ontologyTrack <- rep("No", dim(cds_nn)[1])

    names(ontology) <- c("DC", "NK", "B", "Neu", "Baso", "Mono", "Ery", "Undifferentiated", "Prog", "T")

    number_of_rows <- dim(cds_nn)[1]

    number_of_mono_cells <- 0
    number_of_brokenmono_cells <- 0

    number_of_baso_cells <- 0
    number_of_brokenbaso_cells <- 0

    for(i in 1:number_of_rows) {

        strict_ontology <- ontology[[cds_nn[i, 4]]]
        assigned_ontology <- cds_nn[i ,]


        if(assigned_ontology[4] == 'Mono') number_of_mono_cells <- number_of_mono_cells + 1
        if(assigned_ontology[4] == 'Baso') number_of_baso_cells <- number_of_baso_cells + 1


        if (!all(strict_ontology == assigned_ontology)) {
            print(paste0('Cell ', i, ' has a mismatch'))
            print(paste0('Strict: ', strict_ontology))
            print(paste0('Assigned: ', assigned_ontology))

            if(assigned_ontology[4] == 'Mono') {
                number_of_brokenmono_cells <- number_of_brokenmono_cells + 1
                ontologyTrack[i] <- 'MonoBreak'
            }
            if(assigned_ontology[4] == 'Baso') {
                number_of_brokenbaso_cells <- number_of_brokenbaso_cells + 1
                ontologyTrack[i] <- 'BasoBreak'
            }
        }
    }

    print(paste0('Number of mono cells: ', number_of_mono_cells))
    print(paste0('Number of broken mono cells: ', number_of_brokenmono_cells))
    
    print(paste0('Number of baso cells: ', number_of_baso_cells))
    print(paste0('Number of broken baso cells: ', number_of_brokenbaso_cells))

    return(cbind(cds_nn, ontologyTrack))
}

newFormatData <- function(
    cds_ref,
    cds_qry
) {
    #Delete rows if they have cell type T, Prog, or Undifferentiated

    cds_ref <- cds_ref[, !(colData(cds_ref)$Cell.type.annotation %in% c('T', 'Prog', 'Undifferentiated'))]
    cds_qry <- cds_qry[, !(colData(cds_qry)$Cell.type.annotation %in% c('T', 'Prog', 'Undifferentiated'))]
}