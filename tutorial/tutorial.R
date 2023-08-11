library("dplyr")
library("monocle3")
library("Matrix")
library("assertthat")

source("../R/bayesian_label_transfer_ontology.R")
source("../R/nearest_neighbors.R")

# Theses steps have been done for you, but are important in new analyses
# see monocle3 documentation for more details.
# 
# cds_ref <- preprocess_cds(cds_ref)
# cds_ref <- reduce_dimension(cds_ref, build_nn_index=TRUE)
# save_transform_models(cds_ref, 'cds_ref_R_D_models')
# 
# cds_qry <- load_transform_models(cds_qry, 'cds_ref_R_D_models')
# cds_qry <- preprocess_transform(cds_qry)
# cds_qry <- reduce_dimension_transform(cds_qry)
#
# save_monocle_object(cds_qry, 'cds_qry_R_D')
# save_monocle_object(cds_ref, 'cds_ref_R_D')


# Load the data:
cds_ref <- load_monocle_objects("cds_ref_R_D")
cds_qry <- load_monocle_objects("cds_qry_R_D")

# Run the main function:

cds_qry <- bayesian_ontology_label_transfer(
    cds_query = cds_qry,
    cds_reference = cds_ref,

    reduction_method = "UMAP",
    ref_column_names = c("L1", "L2", "L3", "L4"),
    query_column_names = c("bay_L1", "bay_L2", "bay_L3", "bay_L4"),
    transform_models_dir = 'cds_ref_R_D_models',
    k = 30,
    max_eval = 500,
    nn_control = list(),
    verbose = FALSE
)

# Plot the newly transfer cells:
plot_cells(cds_qry, color_cells_by = "bay_L4")
