library("dplyr")
library("monocle3")
library("Matrix")
library("assertthat")

source("../R/bayesian_label_transfer_ontogeny.R")
source("../R/nearest_neighbors.R")
source("tutorial_functions.R")

# Load the scRNAseq datasets into monolce3 objects:
expression_matrix <- Matrix::readMM('/Users/samsikora/Desktop/Trapnell_Lab/sclassification/bayesian-scRNAseq-label-transfer/cds_data/lineage_pap/downloaded_10x/GSM4185643_stateFate_inVivo_normed_counts.mtx')
expression_matrix <- Matrix::t(expression_matrix)
cell_metadata <- read.delim('/Users/samsikora/Desktop/Trapnell_Lab/sclassification/bayesian-scRNAseq-label-transfer/cds_data/lineage_pap/downloaded_10x/GSM4185643_stateFate_inVivo_metadata.txt')

cds <- new_cell_data_set(expression_matrix, cell_metadata = cell_metadata)

gene_names <- readLines('/Users/samsikora/Desktop/Trapnell_Lab/sclassification/bayesian-scRNAseq-label-transfer/cds_data/lineage_pap/downloaded_10x/GSM4185643_stateFate_inVivo_gene_names.txt')
rownames(cds) <- gene_names

cell_barcodes <- readLines('/Users/samsikora/Desktop/Trapnell_Lab/sclassification/bayesian-scRNAseq-label-transfer/cds_data/lineage_pap/downloaded_10x/GSM4185643_stateFate_inVivo_cell_barcodes.txt')
colnames(cds) <- cell_barcodes

# Split the data into a query and reference dataset by spliting by librarys
librarys <- colData(cds)$Library
unique_librarys <- unique(librarys)

set.seed(001)
qry_librarys <- sample(unique_librarys, 25)
ref_librarys <- unique_librarys[!unique_librarys %in% qry_librarys]

cds_qry <- cds[, !is.na(colData(cds)$Library) & colData(cds)$Library %in% qry_librarys]
cds_ref <- cds[, !is.na(colData(cds)$Library) & colData(cds)$Library %in% ref_librarys]

cds_ref <- estimate_size_factors(cds_ref)
cds_qry <- estimate_size_factors(cds_qry)

# Reduce the dimension of the reference cell data set
# then use that model to reduce the dimension of the query

cds_ref <- preprocess_cds(cds_ref)
cds_ref <- reduce_dimension(cds_ref, build_nn_index=TRUE)
save_transform_models(cds_ref, 'cds_ref_tutorial_model')


cds_qry <- load_transform_models(cds_qry, 'cds_ref_tutorial_model')
cds_qry <- preprocess_transform(cds_qry)
cds_qry <- reduce_dimension_transform(cds_qry)

list_of_cds <- assign_layer_labels(cds_ref, cds_qry, monoProb = 0.3, basoProb = 0.3)

cds_ref <- list_of_cds[[1]]
cds_qry <- list_of_cds[[2]]

# Run the main function to transfer the labels
cds_qry <- bayesian_ontogeny_label_transfer(
    cds_query = cds_qry,
    cds_ref = cds_ref,

    reduction_method = "UMAP",
    ref_column_names = c("L1", "L2", "L3", "L4"),
    query_column_names = c("bay_L1", "bay_L2", "bay_L3", "bay_L4"),
    transform_models_dir = 'cds_ref_tutorial_model',
    k = 30,
    maxeval = 500,
    nn_control = list(),
    verbose = FALSE
)

# Determine what cells have broken ontogeny in both datasets
cds_qry_transfered_labels <- as.data.frame(colData(cds_qry)[, c("bay_L1", "bay_L2", "bay_L3", "bay_L4")])
colData(cds_query) <- cbind(colData(cds_query), check_ontogeny(cds_transfered_labels))

cds_ref_measured_labels <- as.data.frame(colData(cds_ref)[, c("L1", "L2", "L3", "L4")])
colData(cds_ref) <- cbind(colData(cds_ref), check_ontogeny(cds_measured_labels))

#To plot, ensure the rownames are unqiue
rownames(colData(cds_query)) <- paste0(rownames(colData(cds_query)), "_", 1:54277)

#Recreate plots:

#Plot measured cell_types in cds_ref
plot_cells(cds_ref, color_cells_by = "Cell.type.annotation")

#Plot simulated cell_type breakage in cds_ref
plot_cells(cds_ref, color_cells_by = "break")

#Plot estimated cell_types in cds_qry
plot_cells(cds_qry, color_cells_by = "bay_L4")

#Plot simulated cell_type breakage in cds_qry
plot_cells(cds_qry, color_cells_by = "break")