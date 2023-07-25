shuffle_labes <- function(
    cds_ref,
    cds_qry
) {

    cell_labels <- colData(cds_ref)$cell_type_sub


    cell_labels[cell_labels == "enveloping layer (EVL)"] <- "enveloping layer (EVL)/periderm 1"
    cell_labels[cell_labels == "periderm 1"] <- "enveloping layer (EVL)/periderm 1"

    cell_labels[cell_labels == "hatching gland 1"] <- "hatching gland 1/2/periderm 5"
    cell_labels[cell_labels == "hatching gland 2"] <- "hatching gland 1/2/periderm 5"
    cell_labels[cell_labels == "periderm 5"] <- "hatching gland 1/2/periderm 5"

    colData(cds_ref)$cell_type_sub <- cell_labels
}

# CHAT-GPT CODE
# TEST BEFORE DEPLOYING

calculate_cds_nn <- function(ref_coldata, query_search, priors) {
    # Number of cells in the reference data
    number_of_reference_cells <- nrow(ref_coldata)
  
    # Number of cells in the query data
    number_of_cells <- nrow(query_search[['nn.idx']])

    # Number of labels in the reference data
    number_of_labels <- ncol(ref_coldata) - length(ref_column_names)

    # Extract the reference ontology without duplicates
    ref_ontology <- unique(ref_coldata[, ref_column_names])

    # Initialize the output data frame
    cds_nn <- data.frame(matrix(NA, nrow = number_of_cells, ncol = number_of_labels))

    # Iterate over the cells in the query data
    for (i in (number_of_reference_cells + 1):number_of_cells) {
        # Get the nearest neighbors' indices for the current query cell
        nn_indices <- query_search[['nn.idx']][i, ]
      
        # Extract the reference ontology data for the nearest neighbors
        nn_table <- ref_coldata[nn_indices[number_of_reference_cells + 1:number_of_cells], ref_column_names]

        # Calculate the ratio of paths for each label
        ratio_of_paths <- sapply(1:nrow(ref_ontology), function(j) {
            path <- ref_ontology[j, ]
            ratios <- sapply(1:number_of_labels, function(h) {
                sum(nn_table[[colnames(nn_table)[h]]] == path[[colnames(path)[h]]]) / nrow(nn_table)
            })
            return(ratios)
        })

        # Calculate the posteriors for the current query cell
        posteriors <- priors %*% ratio_of_paths

        # Get the index of the maximum posterior probability
        index_of_max <- which.max(posteriors)
        
        # Get the final path using the index
        final_path <- ref_ontology[index_of_max, ]
      
        # Assign the final path to the corresponding row in the output data frame
        cds_nn[i, ] <- final_path
    }

    return(cds_nn)
}
