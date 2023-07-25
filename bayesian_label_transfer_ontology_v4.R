# Sam Sikora Summer Student 2023
#
# !!! A large amount of this file was built on top of Maddy Duran's
# work in label_transfer.R in monocle3 1.3.1 !!!!

# This function is used in train_priors_on_reference's optim
# function. The purpose of this function is to calculate the
# log-likelyhood of the reference data set given the priors.
# This function assumes that the correct path for the cell is
# in the first col of the data matrix. And matrixs with all 1
# or matrixs where the path is not in the matrix are not included.
expectation <- function(
    data,
    par,
    measured_index
) {

    #Ensure parameters are positive, between 0-1, and add up to one
    par <- abs(par)

    par <- par / sum(par)

    # This vector will hold the likelihood of the correct path for
    # each cell in the reference data set.
    expectation_vector <- vector("numeric", length = length(data))

    #Mutiply priors and normalize
    for(i in 1:length(data)) {
        matrix <- data[[i]]
        prob_of_paths <- par %*% matrix
        prob_of_paths <- prob_of_paths / sum(prob_of_paths)
        expectation_vector[i] <- prob_of_paths[[measured_index[i]]]

        #expectation_vector[i] <- prob_of_paths[[1]]
    }

    likelyhood <- sum(log(expectation_vector))

    #Optim wants to minimize, so we return the negative of the likelyhood
    return(-likelyhood)
}

# This function is called by train_priors_on_reference
# the purpose of this function is to create a matrix
# where the columns report how many k labels are reporting
# that path. Every value is (Number of times that label 
# was present in k-NN)/(k). Also, this ensures that the
# correct path is in the first column of the matrix.

#orthology_paths_to_matrix <- function(
#    orthology_paths_,
#    number_of_labels_,
#    measured_,
#    nn_table_
#) {
#
#    #Declare the final matrix
#    ratio_of_paths <- matrix(NA, nrow=number_of_labels_, ncol=nrow(orthology_paths_))
#
#    for(j in 1:nrow(orthology_paths_)) { 
#
#        path <- orthology_paths_[j, ]
#
#        #loop through each layer and calculate the ratio of labels reporting that path.
#        ratios <- vector("numeric", length = number_of_labels_)
#        for(h in 1:number_of_labels_) ratios[h] <- sum(nn_table_[[colnames(nn_table_)[h]]] == path[[colnames(measured_)[h]]]) / nrow(nn_table_)
#
#        #Ensure that the correct path is in the first column of the matrix
#        if(all(measured_ == path) && j != 1) {
#            temp <- ratio_of_paths[, 1]
#            ratio_of_paths[, 1] <- ratios
#            ratio_of_paths[, j] <- temp
#        } else {
#            ratio_of_paths[, j] <- ratios
#        }
#    }
#
#    return(ratio_of_paths)
#}

nn_table_to_matrix <- function(
    ref_ontology,
    nn_table,
    number_of_labels
) {

    list_of_nn_table_colnames <- colnames(nn_table)
    number_of_paths <- nrow(ref_ontology)
    number_of_k_neighbors <- nrow(nn_table)

    ratio_of_paths <- matrix(0, nrow=number_of_labels, ncol=number_of_paths)

    #This chunk will calculate the ratio of paths for unique label in cds_ref
    ratio_of_paths <- sapply(1:number_of_paths, function(j) {
        
        path <- ref_ontology[j, ]
        
        ratios <- sapply(1:number_of_labels, function(h) {
            sum(nn_table[[list_of_nn_table_colnames[h]]] == path[[colnames(path)[h]]]) / number_of_k_neighbors
        })
        
        return(ratios)
    })
    
    return(ratio_of_paths)
}



# This function has the end goal of maximizing the priors_
# The way it does this is that it first calcuates the likelyhood
# that we observe the dataset we have.
# Then it uses the optim function to find the priors_ that maximize
# the likelyhood of the dataset we have.
train_priors_on_reference <- function(
    priors_,
    query_search_,
    ref_coldata_,
    ref_column_names_,
    ref_ontology_,
    number_of_reference_cells_,
    number_of_labels_
) {

    #List of matrixs that report the ratio of labels reporting that path
    list_of_ref_cells_paths <- vector("list", length = number_of_reference_cells_)
    vector_of_measured_index <- vector("numeric", length = number_of_reference_cells_)

    current_index_in_list <- 0

    for(i in 1:number_of_reference_cells_) {
        
        ref_neighbors <- query_search_[['nn.idx']][i,]
        nn_table <- ref_coldata_[ref_neighbors, ref_column_names_]

        measured <- nn_table[1,]
        nn_table <- nn_table[-1,]

        ratio_of_paths <- nn_table_to_matrix(ref_ontology_, nn_table, number_of_labels_)

        #Find index in ontology that matches the measured
        measured_index <- which(apply(ref_ontology_, 1, function(row) all(row == measured)))

        #If all 1 or 0, skip
        measured_ratio_of_paths <- ratio_of_paths[, measured_index]
        if(all(measured_ratio_of_paths == 0) | all(measured_ratio_of_paths == 1) ) next

        #Add the ratio of paths to the list
        current_index_in_list <- current_index_in_list + 1
        
        vector_of_measured_index[current_index_in_list] <- measured_index
        list_of_ref_cells_paths[[current_index_in_list]] <- ratio_of_paths
    }

    #Trunciate the list of ref cells paths to avoid 0's
    list_of_ref_cells_paths <- list_of_ref_cells_paths[1:current_index_in_list]
    vector_of_measured_index <- vector_of_measured_index[1:current_index_in_list]

    #Optimize the priors
    optim_result <- optim(par = priors_, 
                        fn = expectation, 
                        data = list_of_ref_cells_paths, 
                        measured_index = vector_of_measured_index, 
                        method="BFGS")

    #Normalize the priors
    priors_ <- optim_result$par
    priors_ <- abs(priors_)
    priors_ <- priors_ / sum(priors_)

    return(priors_)

}


# Calculate the posteriors. Using the priors that we have optimized
# mutiply that by the ratio of paths for each cell. and report the max
# as the label for that cell.
calculate_posteriors_and_label <- function(
    priors,
    query_search,
    ref_column_names,
    ref_coldata,
    ref_ontology,
    number_of_labels,
    number_of_query_cells,
    number_of_reference_cells,
    number_of_cells
) {

    cds_nn <- data.frame(matrix(NA, nrow=number_of_query_cells, ncol=number_of_labels))
    list_of_final_paths <- vector("list", length = number_of_query_cells)

    for(i in 1:number_of_query_cells) {

        nn_table <- ref_coldata[query_search[['nn.idx']][i + number_of_reference_cells,], ref_column_names]

        #Calculate the ratio of paths for each label
        ratio_of_paths <- nn_table_to_matrix(ref_ontology, nn_table, number_of_labels)

        posteriors <- priors %*% ratio_of_paths

        #Get path name
        index_of_max <- which.max(posteriors)
        final_path <- ref_ontology[index_of_max, ]

        #Add it to the list
        list_of_final_paths[[i]] <- final_path
    }

    #Convert list to dataframe
    cds_nn <- Reduce(rbind, list_of_final_paths)

    return(cds_nn)
}

# k-NN table -> dataframe of labels.
# Train priors on the reference cds
# Calculate the posteriors and label the query cds
get_nn_ontology_cell_labels <- function(
    query_data,
    query_search,
    ref_coldata,
    ref_column_names
) {

    #Get the number of labels and cells
    number_of_cells <- nrow(query_data)
    number_of_reference_cells <- nrow(ref_coldata)
    number_of_query_cells <- number_of_cells - number_of_reference_cells
    number_of_labels <- length(ref_column_names)

    ref_ontology <- ref_coldata[, ref_column_names]
    ref_ontology <- unique(ref_ontology)

    priors <- rep((1/number_of_labels), number_of_labels)

    # Report back an optimized prior by training it on 
    # a k-NN search of the reference data set.
    priors <- train_priors_on_reference(
        priors, 
        query_search, 
        ref_coldata, 
        ref_column_names, 
        ref_ontology,
        number_of_reference_cells, 
        number_of_labels
    )

    # TEMP ----------------------------------------------------------------
    print(paste0("Priors: ", priors))
    # TEMP ----------------------------------------------------------------

    #Use priors to calcuate the posteriors and then find label
    cds_nn <- calculate_posteriors_and_label(
        priors, 
        query_search,
        ref_column_names,
        ref_coldata,
        ref_ontology,
        number_of_labels,
        number_of_query_cells,
        number_of_reference_cells,
        number_of_cells
    )

    return(cds_nn)
}


bayesian_ontology_label_transferv3 <- function(
    cds_query,
    cds_ref,

    reduction_method = c("UMAP", "PCA", "LSI"), 
    ref_column_names,
    query_column_names = ref_column_names,
    transform_models_dir = NULL,
    k = 10,
    nn_control = list(),
    verbose = FALSE

) {

    assertthat::assert_that(methods::is(cds_query, 'cell_data_set'),
                          msg= paste0('cds_query parameter is not a cell_data_set'))
                           
    assertthat::assert_that(methods::is(cds_ref, 'cell_data_set'),
                        msg= paste0('cds_ref parameter is not a cell_data_set'))

    assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
                error = function(e) FALSE),
    msg = "reduction_method must be 'UMAP', 'PCA', or 'LSI'")

    assertthat::assert_that(assertthat::is.count(k))

    reduction_method <- match.arg(reduction_method)
    ref_coldata <- colData(cds_ref)

    if(!is.data.frame(ref_coldata)) ref_coldata <- as.data.frame(ref_coldata)

    assertthat::assert_that(all(ref_column_names %in% colnames(ref_coldata)),
                            msg= paste0('ref_column_name \'', ref_column_names, '\' is not in the ref_col_data'))

    #Added checks:

    assertthat::assert_that(length(ref_column_names) > 1, 
                            msg = "Length of ref_column_names must be greater than 1.")
    
    assertthat::assert_that(length(ref_column_names) == length(unique(ref_column_names)),
                            msg = "ref_column_names must be unique.")

    assertthat::assert_that(length(query_column_names) == length(unique(query_column_names)),
                            msg = "query_column_names must be unique.")

    assertthat::assert_that(length(ref_column_names) == length(query_column_names),
                            msg = "ref_column_names and query_column_names must be the same length.")

    for(column_name in ref_column_names) {
        assertthat::assert_that(is.character(ref_coldata[[column_name]]),
                                msg = paste0('ref_coldata column \'', column_name, '\' is not a character vector'))
    }

    if(reduction_method == 'UMAP') {
        nn_control_default <- list(method='annoy', metric='euclidean', n_trees=50, M=48, ef_construction=200, ef=150, grain_size=1, cores=1)
    } else {
        nn_control_default <- list(method='annoy', metric='cosine', n_trees=50, M=48, ef_construction=200, ef=150, grain_size=1, cores=1)
    }


  # Use set_nn_control to find nn method, which we need in order to select the correct index,
  # and we may need the index to set the nn_control[['search_k']].
    nn_control_tmp <- set_nn_control(mode=2,
                                   nn_control=nn_control,
                                   nn_control_default=nn_control_default,
                                   nn_index=NULL,
                                   k=k,
                                   verbose=verbose)

    #To make a nn index with both the reference and query data, we need to combine the cds

    cds_ref_temp <- cds_ref
    cds_query_temp <- cds_query
    
    colData(cds_ref_temp)[['data_set']] <- 'reference'
    colData(cds_query_temp)[['data_set']] <- 'query'
    cds_com <- combine_cds(list(cds_ref_temp, cds_query_temp), keep_all_genes=TRUE, cell_names_unique=TRUE, keep_reduced_dims=TRUE)
    cds_com <- load_transform_models(cds_com, directory_path=transform_models_dir)
    cds_com <- preprocess_cds(cds_com)

    cds_nn_index <- get_cds_nn_index(cds=cds_com, reduction_method=reduction_method, nn_control_tmp[['method']], verbose=verbose) 

    cds_reduced_dims <- SingleCellExperiment::reducedDims(cds_com)[[reduction_method]]

    if(ncol(cds_reduced_dims) != cds_nn_index[['ncol']]) {
        stop('transfer_cell_labels: reduced dimension matrix and nearest neighbor index dimensions do not match')
    }

    checksum_matrix_rownames <- cds_nn_index[['checksum_rownames']]
    if(!is.na(checksum_matrix_rownames)) {
        checksum_coldata_rownames <- digest::digest(sort(rownames(ref_coldata)))
        if(checksum_matrix_rownames != checksum_coldata_rownames) {
            # In earlier versions (<2022-08-29), I did not sort the rownames. Preserve compatibility.
            checksum_coldata_rownames <- digest::digest(rownames(ref_coldata))
            if(checksum_matrix_rownames != checksum_coldata_rownames) {
            stop('transfer_ontology_cell_labels: matrix and colData rownames do not match')
            }
        }
    } else if(!is.na(cds_nn_index[['nrow']]) && (nrow(ref_coldata) != cds_nn_index[['nrow']])) {
        stop('transfer_ontology_cell_labels: matrix and colData row counts do not match')
    }

    nn_control <- set_nn_control(mode=2,
                               nn_control=nn_control,
                               nn_control_default=nn_control_default,
                               nn_index=cds_nn_index,
                               k=k,
                               verbose=verbose)

    # Load the reference projection models and nn indexes
    # into the query cds.
    if(!is.null(transform_models_dir)) {
        #cds_query <- load_transform_models(cds=cds_query, directory_path=transform_models_dir)
        cds_com <- load_transform_models(cds=cds_com, directory_path=transform_models_dir)
    }
  
    assertthat::assert_that(!is.null(cds_query@reduce_dim_aux[[reduction_method]]),
                            msg=paste0("Reduction Method '", reduction_method, "' is not in the",
                                    "loaded model object."))


    # Search the reference reduction_method space for nearest neighbors
    # to the query cells.
    # The cds_reduced_dims contains the query cell coordinates
    # after projection into the reference space.
    # The cds@reduce_dim_aux[[reduction_method]] contains the reduction_method
    # coordinates for the reference data set, which were
    # loaded using load_transform_models() above.
    cds_res <- search_nn_index(query_matrix=cds_reduced_dims, nn_index=cds_nn_index,
                             k=k, nn_control=nn_control, verbose=verbose)


    
    cds_nn <- get_nn_ontology_cell_labels(
        query_data=cds_reduced_dims,
        query_search=cds_res,
        ref_coldata=ref_coldata,
        ref_column_names=ref_column_names
    )


    colnames(cds_nn) <- query_column_names
    
    colData(cds_query) <- cbind(colData(cds_query), cds_nn)

    return(cds_query)
}