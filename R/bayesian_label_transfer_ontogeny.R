  # Sam Sikora Summer Student 2023
  #
  # !!! This is an extension of Maddy Duran's
  # work in label_transfer.R in monocle3 1.3.1 !!!!

  # This function takes in a list of matrices. Each matrix
  # is a ratio of paths for a cell. The function will
  # return the summaiton hinge loss of the priors and the ratio of
  # paths for each cell. Called by train_priors_on_reference.
  hinge_loss <- function(
      x0,
      data,
      measured_index
  ) {

    #Ensure priors are positive and add up to one
    par <- abs(x0)
    par <- par / sum(par)
    
    #Create vector to store the likelyhood of each path
    expectation_vector <- vector("numeric", length = length(data))
    
    for(i in 1:length(data)) {

      #Mutiply priors and normalize
      matrix <- data[[i]]
      prob_of_paths <- par %*% matrix

      #Soft Max the resulting variables
      prob_of_paths <- exp(prob_of_paths) / sum(exp(prob_of_paths))

      # Loss is calcuated by the distance between the maxium path and
      # the measured path.
      expectation <- 1 - (max(prob_of_paths) - prob_of_paths[[measured_index[i]]])

      #If there is a tie, return 0.5. We do not want ties.
      if(length(prob_of_paths[prob_of_paths == max(prob_of_paths)]) > 1) expectation <- 0.5

      expectation_vector[i] <- expectation
    }

    return(-sum(expectation_vector))
  }


  # This function takes in a k-NN dataframe. It then returns
  # a matrix where each coloumn is a path down the ontogeny,
  # each row is a layer, and every entry is how many times
  # that label was seen in the k-NN dataframe. Called by
  # train_priors_on_reference and calculate_posteriors_and_label.

  nn_table_to_matrix <- function(
      ref_ontogeny,
      nn_table,
      NUMBER_OF_LABELS
  ) {
    
    list_of_nn_table_colnames <- colnames(nn_table)
    NUMBER_OF_ONTOGENY_PATHS <- nrow(ref_ontogeny)
    NUMBER_OF_NEAREST_NEIGHBORS <- nrow(nn_table)
    
    ratio_of_paths <- matrix(0, nrow=NUMBER_OF_LABELS, ncol=NUMBER_OF_ONTOGENY_PATHS)

    #This chunk will calculate the ratio of paths for unique label in cds_ref
    ratio_of_paths <- sapply(1:NUMBER_OF_ONTOGENY_PATHS, function(j) {
      path <- ref_ontogeny[j, ]

      ratios <- sapply(1:NUMBER_OF_LABELS, function(h) {
        sum(nn_table[[list_of_nn_table_colnames[h]]] == path[[colnames(path)[h]]]) / NUMBER_OF_NEAREST_NEIGHBORS
      })
      return(ratios)
    })
    
    return(ratio_of_paths)
  }


  # This function has the end goal of maximizing the priors
  # on the reference data set. The first part of the function
  # will setup the data so it can used to evalute the priors.
  # The second part of the function will first globally
  # optimize the priors, then use a local optmizer to fine
  # tune the priors. The function will return the optimized
  # priors.
  train_priors_on_reference <- function(
      priors,
      query_search,
      ref_coldata,
      ref_column_names,
      ref_ontogeny,
      maxeval,
      NUMBER_OF_REFERENCE_CELLS,
      NUMBER_OF_LABELS
  ) {
    
    #List of matrixs that report the ratio of labels reporting that path
    list_of_ref_cells_paths <- vector("list", length = NUMBER_OF_REFERENCE_CELLS)
    vector_of_measured_index <- vector("numeric", length = NUMBER_OF_REFERENCE_CELLS)
    
    current_index_in_list <- 0
    
    for(i in 1:NUMBER_OF_REFERENCE_CELLS) {
      
      #Search for k-NN
      ref_neighbors <- query_search[['nn.idx']][i,]
      nn_table <- ref_coldata[ref_neighbors, ref_column_names]

      #If all the neighbors are the same, skip (arent important in optmizing)
      if(dim(unique(nn_table))[1] == 1) next
      
      #Get the measured path and the index of that path in the ontogeny
      measured <- nn_table[1,]
      nn_table <- nn_table[-1,]

      #Make sure measured path isnt all zeros no reason to optimize
      if(nrow(merge(measured, nn_table)) == 0) next

      #Get index of path in ontogeny
      measured_index <- which(apply(ref_ontogeny, 1, function(row) all(row == measured)))

      #Turn the nn_table to a matrix of ratios of paths
      ratio_of_paths <- nn_table_to_matrix(ref_ontogeny, nn_table, NUMBER_OF_LABELS)
      
      #Add the ratio of paths to the list
      current_index_in_list <- current_index_in_list + 1

      vector_of_measured_index[current_index_in_list] <- measured_index
      list_of_ref_cells_paths[[current_index_in_list]] <- ratio_of_paths
    }

    # If there were no cells that were not all 0's or all 1's just transfer
    # at the most speficic label.
    if(current_index_in_list == 0) return(c(rep(0, NUMBER_OF_LABELS - 1), 1))

    #Trunciate the list of ref cells paths to avoid NA's
    list_of_ref_cells_paths <- list_of_ref_cells_paths[1:current_index_in_list]
    vector_of_measured_index <- vector_of_measured_index[1:current_index_in_list]

    #Globally optimize the priors
    optim_result_GN <- nloptr(
      opts = list("algorithm"="NLOPT_GN_ESCH", "xtol_rel"=1.0e-4, "maxeval"=maxeval),

      x0 = priors,
      lb = rep(0, NUMBER_OF_LABELS),
      ub = rep(1, NUMBER_OF_LABELS),

      eval_f = hinge_loss,
      data = list_of_ref_cells_paths, 
      measured_index = vector_of_measured_index
    )

    #Normalize the priors
    priors <- optim_result_GN$solution
    priors <- abs(priors)
    priors <- priors / sum(priors)

    #Locally optimize the priors
    optim_result_LN <- nloptr(
      opts = list("algorithm"="NLOPT_LN_BOBYQA", "xtol_rel"=1.0e-8, "maxeval"=maxeval),

      x0 = priors,
      lb = rep(0, NUMBER_OF_LABELS),
      ub = rep(1, NUMBER_OF_LABELS),

      eval_f = hinge_loss,
      data = list_of_ref_cells_paths, 
      measured_index = vector_of_measured_index
    )

    #Normalize the priors
    priors <- optim_result_LN$solution
    priors <- abs(priors)
    priors <- priors / sum(priors)

    return(priors) 
  }

  # Calculate the posteriors. Using the priors that we have optimized
  # mutiply that by the ratio of paths for each cell. and report the max
  # as the label for that cell.
  calculate_posteriors_and_label <- function(
      priors,
      query_search,
      ref_column_names,
      ref_coldata,
      ref_ontogeny,
      NUMBER_OF_LABELS,
      NUMBER_OF_QUERY_CELLS,
      NUMBER_OF_REFERENCE_CELLS,
      NUMBER_OF_CELLS
  ) {
    
    cds_nn <- data.frame(matrix(NA, nrow=NUMBER_OF_QUERY_CELLS, ncol=NUMBER_OF_LABELS))
    list_of_final_paths <- vector("list", length = NUMBER_OF_QUERY_CELLS)
    
    for(i in 1:NUMBER_OF_QUERY_CELLS) {
      
      nn_table <- ref_coldata[query_search[['nn.idx']][i + NUMBER_OF_REFERENCE_CELLS,], ref_column_names]
      
      #Calculate the ratio of paths for each label
      ratio_of_paths <- nn_table_to_matrix(ref_ontogeny, nn_table, NUMBER_OF_LABELS)
      
      posteriors <- priors %*% ratio_of_paths
      
      #Get path name
      index_of_max <- which.max(posteriors)
      final_path <- ref_ontogeny[index_of_max, ]
      
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
  get_nn_ontogeny_cell_labels <- function(
      query_data,
      query_search,
      ref_coldata,
      ref_column_names,
      maxeval
  ) {
    
    #Get the number of labels and cells
    NUMBER_OF_CELLS <- nrow(query_data)
    NUMBER_OF_REFERENCE_CELLS <- nrow(ref_coldata)
    NUMBER_OF_QUERY_CELLS <- NUMBER_OF_CELLS - NUMBER_OF_REFERENCE_CELLS
    NUMBER_OF_LABELS <- length(ref_column_names)
    
    ref_ontogeny <- ref_coldata[, ref_column_names]
    ref_ontogeny <- unique(ref_ontogeny)
    ref_ontogeny <- as.data.frame(ref_ontogeny)
    
    priors <- rep((1/NUMBER_OF_LABELS), NUMBER_OF_LABELS)
    
    # Report back an optimized prior by training it on 
    # a k-NN search of the reference data set.
    priors <- train_priors_on_reference(
      priors, 
      query_search, 
      ref_coldata, 
      ref_column_names, 
      ref_ontogeny,
      maxeval,
      NUMBER_OF_REFERENCE_CELLS, 
      NUMBER_OF_LABELS
    )

    for(i in 1:NUMBER_OF_LABELS) {
      print(paste0("Prior for ", ref_column_names[i], ": ", priors[i]))
    }
    
    #Use priors to calcuate the posteriors and then find label
    cds_nn <- calculate_posteriors_and_label(
      priors, 
      query_search,
      ref_column_names,
      ref_coldata,
      ref_ontogeny,
      NUMBER_OF_LABELS,
      NUMBER_OF_QUERY_CELLS,
      NUMBER_OF_REFERENCE_CELLS,
      NUMBER_OF_CELLS
    )
    
    return(cds_nn)
  }


  #' @title Transfer ontogeny labels
  #'
  #' @description Transfer ontogeny labels from
  #' a reference to a query dataset. To do this,
  #' first, cds_qry and cds_ref are combined into
  #' a combo cds. Then, the combo cds is reduced
  #' using the reduction_method. Using k-NN the combo
  #' cds is compared to cds_ref. 
  #' 
  #' Then, priors for each layer of the ontogeny is
  #' optimized on the reference data set. The priors
  #' are then used to calculate the posteriors for each
  #' cell in the query data set. Finally, the cell is
  #' assigned a label based on the max posterior.
  #' 
  #' @param cds_query A cell_data_set object. An unkown
  #' set of cells to label
  #' @param cds_ref A cell_data_set object. Know set of
  #' cells to base the labels off of.
  #' @param reduction_method The method used to reduce
  #' the dimensionality of the combo cds. Must be one of
  #' 'UMAP', 'PCA', or 'LSI'.
  #' @param ref_column_names The column names of colData(cds_ref)
  #' that contain the ontogeny. The order of the column names
  #' must be from most broad to most speficic.
  #' @param query_column_names (Optional) The column names of
  #' colData(cds_qry) that will be affixed to the cell_data_set.
  #' @param transform_models_dir (Optional) The directory path
  #' to the transform models. If NULL, then the transform models
  #' will be loaded from the cds_ref.
  #' @param k The number of nearest neighbors to use in the
  #' k-NN search.
  #' @param nn_control (Optional) A list of parameters to control
  #' the k-NN search.
  #' @param maxeval (Optional) The maximum number of iterations
  #' to use in the optimization of the priors. Higher means more accurate
  #' priors, but slower function.
  #' @param verbose (Optional) A boolean. If TRUE, then the function
  #' will print out progress messages.
  #' @return A cell_data_set object with the query_column_names
  #' affixed to the colData.
  #' 
  bayesian_ontogeny_label_transfer <- function(
    cds_query,
    cds_ref,
    
    reduction_method = c("UMAP", "PCA", "LSI"), 
    ref_column_names,
    query_column_names = ref_column_names,
    transform_models_dir = NULL,
    k = 10,
    maxeval = 500,
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
    ref_coldata <- as.data.frame(colData(cds_ref))
    
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
    
    # To make a nn index with both the reference and query data, we need to combine the cds
    # However, the are two important things to ensure. First, the column names of the
    # reference and query data must be the same. Second, every ID of a cell in that
    # data is unique.

    cds_ref_temp <- cds_ref
    cds_query_temp <- cds_query

    colData(cds_ref_temp) <- colData(cds_ref_temp)[, ref_column_names]
    for(col in ref_column_names) {
      colData(cds_query_temp)[[col]] <- NULL
      colData(cds_query_temp)[[col]] <- rep(NA, ncol(cds_query_temp))
    }
    colData(cds_query_temp) <- colData(cds_query_temp)[, ref_column_names]

    rownames(colData(cds_query_temp)) <- paste0("cds_qry", 1:dim(cds_query_temp)[2])
    rownames(colData(cds_ref_temp)) <- paste0("cell_ref", 1:dim(cds_ref_temp)[2])
    colnames(cds_query_temp) <- paste0("cds_qry", 1:dim(cds_query_temp)[2])
    colnames(cds_ref_temp) <- paste0("cell_ref", 1:dim(cds_ref_temp)[2])

    cds_com <- combine_cds(list(cds_ref_temp, cds_query_temp), keep_reduced_dims = TRUE, cell_names_unique = TRUE)
    cds_com <- load_transform_models(cds_com, directory_path=transform_models_dir)
    cds_com <- preprocess_cds(cds_com)
    
    cds_nn_index <- get_cds_nn_index(cds=cds_com, reduction_method=reduction_method, nn_control_tmp[['method']], verbose=verbose) 
    
    cds_reduced_dims <- SingleCellExperiment::reducedDims(cds_com)[[reduction_method]]
    
    if(ncol(cds_reduced_dims) != cds_nn_index[['ncol']]) {
      stop('transfer_cell_labels: reduced dimension matrix and nearest neighbor index dimensions do not match')
    }
    
    nn_control <- set_nn_control(mode=2,
                                nn_control=nn_control,
                                nn_control_default=nn_control_default,
                                nn_index=cds_nn_index,
                                k=k,
                                verbose=verbose)
    

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
    
    cds_nn <- get_nn_ontogeny_cell_labels(
      query_data=cds_reduced_dims,
      query_search=cds_res,
      ref_coldata=ref_coldata,
      ref_column_names=ref_column_names,
      maxeval=maxeval
    )

    #cds_nn <- check_ontogeny(cds_nn)
    #colnames(cds_nn) <- c(query_column_names, "break")
    #rownames(colData(cds_query)) <- 1:54277
    #colData(cds_query) <- cbind(colData(cds_query), cds_nn)
    #plot_cells(cds_query, color_cells_by="break")

    colnames(cds_nn) <- query_column_names
    
    colData(cds_query) <- cbind(colData(cds_query), cds_nn)
    
    return(cds_query)
  }