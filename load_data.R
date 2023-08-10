library("dplyr")
library("monocle3")
library("Matrix")
library("assertthat")
library("nloptr")

source("nearest_neighbors.R")
source("bayesian_label_transfer_ontology_v4.R")
source("shuffle.R")
source("format_data.R")

doItAll <- function() {
    load_data()
    benchmark()
}


L <- function() {
    cds_ref <<- load_monocle_objects("cds_ref")
    cds_qry <<- load_monocle_objects("cds_qry")
}

LD <- function(monoProb = 0.3, basoPorb= 0.3) {
    cds_ref <- load_monocle_objects("lineage_pap/cds_ref_R_D")
    cds_qry <- load_monocle_objects("lineage_pap/cds_qry_R_D")
    list_of_cds <- assign_layer_labels(cds_ref, cds_qry, monoProb, basoPorb)
    cds_ref <<- list_of_cds[[1]]
    cds_qry <<- list_of_cds[[2]]
}


benchmark <- function() {

    source("bayesian_label_transfer_ontology_v4.R")

    cds_nn2 <- bayesian_ontology_label_transferv3(
        cds_query = cds_qry,
        cds_ref = cds_ref,
        reduction_method="UMAP",
        ref_column_names=c("germ_layer", "tissue", "cell_type_broad", "cell_type_sub"),
        query_column_names = c("bay_germ", "bay_tissue", "bay_broad", "bay_sub"),
        transform_models_dir = 'cds_ref_test_models',
        k = 25,
        nn_control = list(),
        verbose = FALSE
    )
        
    cds_nn_non_ortho <- transfer_cell_labels(
        cds_qry, 
        reduction_method='UMAP', 
        ref_coldata=colData(cds_ref), 
        ref_column_name="cell_type_sub",
        query_column_name ='temp_single', 
        transform_models_dir='cds_ref_test_models',
        k = 20,
        top_frac_threshold = 0.7,
        top_next_ratio_threshold = 1.7
    )

    cds_qry_lab_fix <- fix_missing_cell_labels(cds_nn_non_ortho, reduction_method='UMAP', from_column_name='temp_single', to_column_name='single_transfer')

    cds_nn$single_transfer <- cds_qry_lab_fix$single_transfer

    return(cds_nn)
}

run_main_bay <- function() {

    cds_nn <- bayesian_transfer_cell_labels(
        cds_qry,
        cds_ref,
        "UMAP",
        colData(cds_ref),
        "cell_type_sub",
        transform_models_dir = 'cds_ref_test_models',
        k = 20
    )
}

load_data_bay <- function() {

    gapData <- readRDS("./full_gap_hf_ctrl_ref_mito-filt_1.25M_model-update_anno_cds.RDS")
    cds_ref <- gapData[, !is.na(colData(gapData)$gene_target) & !is.na(colData(gapData)$expt) & colData(gapData)$expt %in% c("GAP13", "GAP14", "GAP18") & colData(gapData)$gene_target %in% c("ctrl-uninj" , "ctrl-inj")]
    cds_qry <- gapData[, !is.na(colData(gapData)$gene_target) & !is.na(colData(gapData)$expt) & colData(gapData)$expt %in% c("GAP16") & colData(gapData)$gene_target %in% c("ctrl-uninj" , "ctrl-inj")]

    cds <- cds_ref

    assigned_type_marker_test_res <- top_markers(cds,
                                             group_cells_by="cell_type_broad",
                                             reference_cells=1000,
                                             cores=8)
    
    garnett_markers <<- assigned_type_marker_test_res %>%
                        filter(marker_test_q_value < 0.01 & specificity >= 0.5) %>%
                        dplyr::group_by(cell_group) %>%
                        dplyr::top_n(5, marker_score)

    garnett_markers <- garnett_markers %>% 
                        dplyr::group_by(gene_short_name) %>%
                        filter(n() == 1)
}


load_data <- function() {

    gapData <<- readRDS("./full_gap_hf_ctrl_ref_mito-filt_1.25M_model-update_anno_cds.RDS")

    #For testing, the cell data set reference is the 72hpf and the query is the 48hp

    #cds_ref <<- gapData[, !is.na(colData(gapData)$timepoint) & colData(gapData)$timepoint == 72 ]
    #cds_qry <<- gapData[, !is.na(colData(gapData)$timepoint) & colData(gapData)$timepoint == 48 ]

    cds_ref <<- gapData[, !is.na(colData(gapData)$gene_target) & !is.na(colData(gapData)$expt) & colData(gapData)$expt %in% c("GAP13", "GAP14", "GAP18") & colData(gapData)$gene_target %in% c("ctrl-uninj" , "ctrl-inj")]
    cds_qry <<- gapData[, !is.na(colData(gapData)$gene_target) & !is.na(colData(gapData)$expt) & colData(gapData)$expt %in% c("GAP16") & colData(gapData)$gene_target %in% c("ctrl-uninj" , "ctrl-inj")]
    cds_com <<- gapData[, !is.na(colData(gapData)$gene_target) & !is.na(colData(gapData)$expt) & colData(gapData)$expt %in% c("GAP13", "GAP14", "GAP16", "GAP18") & colData(gapData)$gene_target %in% c("ctrl-uninj" , "ctrl-inj")]

    genes_ref <- row.names(cds_ref)
    genes_qry <- row.names(cds_qry)

    genes_shared <- intersect(genes_ref, genes_qry)
    cds_ref <<- cds_ref[genes_shared,]
    cds_qry <<- cds_qry[genes_shared,]

    cds_qry <<- cds_qry[, colData(cds_qry)$n.umi >= 1000]
    cds_ref <<- cds_ref[, colData(cds_ref)$n.umi >= 1000]

    cds_ref <<- estimate_size_factors(cds_ref)
    cds_qry <<- estimate_size_factors(cds_qry)

    reducedDim(cds_ref) <<- NULL
    reducedDim(cds_ref) <<- NULL
    reducedDim(cds_ref) <<- NULL

    reducedDim(cds_qry) <<- NULL
    reducedDim(cds_qry) <<- NULL
    reducedDim(cds_qry) <<- NULL

    cds_ref <<- preprocess_cds(cds_ref, num_dim=100)
    cds_ref <<- reduce_dimension(cds_ref, build_nn_index=TRUE)

    save_transform_models(cds_ref, 'cds_ref_test_models')
    cds_qry <<- load_transform_models(cds_qry, 'cds_ref_test_models')

    cds_qry <<- preprocess_transform(cds_qry)
    cds_qry <<- reduce_dimension_transform(cds_qry)

    colData(cds_ref)[['data_set']] <<- 'reference'
    colData(cds_qry)[['data_set']] <<- 'query'
}



print_pdf <- function() {

    gapData <<- readRDS("./full_gap_hf_ctrl_ref_mito-filt_1.25M_model-update_anno_cds.RDS")

    cds_com <<- gapData[, !is.na(colData(gapData)$gene_target) & !is.na(colData(gapData)$expt) & colData(gapData)$expt %in% c("GAP13", "GAP14", "GAP16", "GAP18") & colData(gapData)$gene_target %in% c("ctrl-uninj" , "ctrl-inj")]

    cds_com <<- cds_com[, colData(cds_com)$n.umi >= 1000]
    cds_com <<- estimate_size_factors(cds_com)

    reducedDim(cds_com) <<- NULL
    reducedDim(cds_com) <<- NULL
    reducedDim(cds_com) <<- NULL

    cds_com <<- preprocess_cds(cds_com, num_dim=100)
    cds_com <<- reduce_dimension(cds_com)

    cds_com <<- cluster_cells(cds_com, resolution=1e-5)

    colData(cds_com)[['data_set']] <<- ifelse((cds_com$expt %in% c("GAP13", "GAP14", "GAP18")), 'reference', 'query')

    pdf(file="./figures/colored_by_cluster.pdf")
    plot_cells(cds_com)
    dev.off()

    pdf(file="./figures/colored_by_data_set.pdf")
    plot_cells(cds_com, color_cells_by='data_set')
    dev.off()

    pdf(file="./figures/colored_by_cell_type_only_reference.pdf")
    plot_cells(cds_com[, cds_com$data_set == 'reference'], color_cells_by='cell_type_broad')
    dev.off()


    cds_ref <- cds_com[, !is.na(colData(cds_com)$gene_target) & !is.na(colData(cds_com)$expt) & colData(cds_com)$expt %in% c("GAP13", "GAP14", "GAP18") & colData(cds_com)$gene_target %in% c("ctrl-uninj" , "ctrl-inj")]
    cds_qry <- cds_com[, !is.na(colData(cds_com)$gene_target) & !is.na(colData(cds_com)$expt) & colData(cds_com)$expt %in% c("GAP16") & colData(cds_com)$gene_target %in% c("ctrl-uninj" , "ctrl-inj")]
    
}

V <- function() {
    cds_query <<- cds_qry
    reduction_method <<- "UMAP"
    ref_column_names <<- c("germ_layer", "tissue", "cell_type_broad", "cell_type_sub")
    query_column_names <<- c("bay_germ", "bay_tissue", "bay_broad", "bay_sub")
    transform_models_dir <<- 'cds_ref_test_models'
    ref_column_ontology = list()
    k <<- 30
    verbose <<- FALSE
    nn_control <<- list()
    ref_column_ontology <<- list()
  
}

VD <- function() {
    cds_query <<- cds_qry
    reduction_method <<- "UMAP"
    ref_column_names <<- c("L1", "L2", "L3", "L4")
    query_column_names <<- c("bay_L1", "bay_L2", "bay_L3", "bay_L4")
    transform_models_dir <<- 'lineage_pap/cds_ref_R_D_models'
    maxeval <<- 500
    k <<- 30
    verbose <<- FALSE
    nn_control <<- list()
    ref_column_ontology <<- list()
}


givemeanum <- function(x) {
    x <- x + 15479
    ref_neighbors <- query_search[['nn.idx']][x,]
    ref_labels <- ref_coldata[ref_neighbors, ref_column_name]
    return(ref_labels)
}


R <- function() {
    cds_qry <- bayesian_ontology_label_transferv3(
        cds_query = cds_qry,
        cds_ref = cds_ref,
        reduction_method="UMAP",
        ref_column_names=c("L1", "L2", "L3", "L4"),
        query_column_names = c("bay_L1", "bay_L2", "bay_L3", "bay_L4"),
        transform_models_dir = 'lineage_pap/cds_ref_R_D_models',
        k = 30
    )
}





cds_to_csv <- function(cds, everything = TRUE) {
    if (everything) write.csv(as.data.frame(colData(cds)), "cds_colData.csv")
    else {
        write.csv(as.data.frame(colData(cds)[c("germ_layer", "tissue", "cell_type_broad", "cell_type_sub", "bay_germ", "bay_tissue", "bay_broad", "bay_sub", "single_transfer")]), "cds_colData.csv")
    }
}


load_tut_data <- function() {
    library(monocle3)
    library(Matrix)

    # Load the reference data set.
    matrix_ref <- readMM(gzcon(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/mouse/data/cao.mouse_embryo.sample.mtx.gz")))
    cell_ann_ref <- read.csv(gzcon(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/mouse/data/cao.mouse_embryo.sample.coldata.txt.gz"), text=TRUE), sep='\t')
    gene_ann_ref <- read.csv(gzcon(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/mouse/data/cao.mouse_embryo.sample.rowdata.txt.gz"), text=TRUE), sep='\t')

    cds_ref <<- new_cell_data_set(matrix_ref,
                                cell_metadata = cell_ann_ref,
                                gene_metadata = gene_ann_ref)

    # Load the query data set.
    matrix_qry <- readMM(gzcon(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/mouse/data/srivatsan.mouse_embryo_scispace.sample.mtx.gz")))
    cell_ann_qry <- read.csv(gzcon(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/mouse/data/srivatsan.mouse_embryo_scispace.sample.coldata.txt.gz"), text=TRUE), sep='\t')
    gene_ann_qry <- read.csv(gzcon(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/mouse/data/srivatsan.mouse_embryo_scispace.sample.rowdata.txt.gz"), text=TRUE), sep='\t')

    cds_qry <<- new_cell_data_set(matrix_qry,
                                cell_metadata = cell_ann_qry,
                                gene_metadata = gene_ann_qry)

    # Genes in reference.
    genes_ref <- row.names(cds_ref)

    # Genes in query.
    genes_qry <- row.names(cds_qry)

    # Shared genes.
    genes_shared <- intersect(genes_ref, genes_qry)

    # Remove non-shared genes.
    cds_ref <<- cds_ref[genes_shared,]
    cds_qry <<- cds_qry[genes_shared,]

    # Reference data set UMI cutoff.
    numi_ref <- min(colData(cds_ref)[['Total_mRNAs']])
    # numi_ref is 1001
    numi_qry <- min(colData(cds_qry)[['n.umi']])
    # numi_qry is 1000

    cds_ref <<- estimate_size_factors(cds_ref)
    cds_qry <<- estimate_size_factors(cds_qry)

    cds_ref <<- preprocess_cds(cds_ref, num_dim=100)
    cds_ref <<- reduce_dimension(cds_ref, build_nn_index=TRUE)
    # Save the PCA and UMAP transform models for use with projection.
    save_transform_models(cds_ref, 'cds_ref_test_models')

    # Load the reference transform models into the query cds.
    cds_qry <<- load_transform_models(cds_qry, 'cds_ref_test_models')
    # Apply the reference transform models to the query cds.
    cds_qry <<- preprocess_transform(cds_qry)
    cds_qry <<- reduce_dimension_transform(cds_qry)

    colData(cds_ref)[['data_set']] <<- 'reference'
    colData(cds_qry)[['data_set']] <<- 'query'
}


T <- function() {

    source('bayesian_label_transfer_ontology_v4.R')

  assertthat::assert_that(methods::is(cds_query, 'cell_data_set'),
                          msg= paste0('cds_query parameter is not a cell_data_set'))
  
  assertthat::assert_that(methods::is(cds_ref, 'cell_data_set'),
                          msg= paste0('cds_ref parameter is not a cell_data_set'))
  
  assertthat::assert_that(assertthat::is.count(k))
  
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
  
  #To make a nn index with both the reference and query data, we need to combine the cds
  cds_com <- combine_cds(list(cds_ref, cds_query), keep_reduced_dims = TRUE)
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

    query_data <<- cds_reduced_dims
    query_search <<- cds_res
    ref_coldata <<- ref_coldata
    ref_column_names <<- ref_column_names
    ref_column_ontology <<-ref_column_ontology
    cds_reduced_dims <<- cds_reduced_dims
}

expectationNumberRight <- function(
  data,
  par,
  measured_index
) {
  #This is a benchmark to detemrine how many ref cells were got right with priors

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

    expectation_vector[i] <- ifelse(max(prob_of_paths) == prob_of_paths[[measured_index[i]]], 1, 0)
  }

  return(sum(expectation_vector))
  
}