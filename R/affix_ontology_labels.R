
#' @title Affix full layer labels to a cell_data_set.
#'
#' @description For every cell label in a cell_data_set,
#' this function will find the full path of that in label
#' within ontology and affix it to the cell_data_set.
#' 
#' The ontology is found within ref_coldata. 
#' 
#' *** The function assumes that there is only one path to 
#' get to the most speficic cell label and that every path 
#' is given in the ref_coldata. ***
#' 
#' @param cds_qry A cell_data_set object. To pull cells
#' to label from
#' @param ref_coldata colData(reference) containing 
#' the ontology of the cell labels.
#' @param ref_column_names The column names of ref_coldata
#' that contain the ontology. The order of the column names
#' must be from most broad to most speficic.
#' @param query_column_name The column name of colData(cds_qry)
#' that contains the cell labels to be affixed.
#' @param query_column_names (Optional) The column names of
#' colData(cds_qry) that will be affixed to the cell_data_set.
#' @return a cell data object
#' @examples
#' \dontrun{
#'   cds_qry <- affix_ontology_labels(
#'     cds_qry,
#'     colData(cds_ref),
#'     c("germ_layer", "tissue", "cell_type_broad", "cell_type_sub"),
#'     "cell_type_sub",
#'     c("transfered_germ_layer", "transfered_tissue", 
#'          "transfered_cell_type_broad", "transfered_cell_type_sub")
#' }
#'

affix_ontology_labels <- function(
    cds_qry,
    ref_coldata,
    ref_column_names,
    query_column_name,

    query_column_names = ref_column_names
) {

    assertthat::assert_that(
        methods::is(cds_qry, 'cell_data_set'),
        is.character(ref_column_names),
        is.character(query_column_name),
        is.character(query_column_names)
    , msg = "affix_ontology_labels: Invalid input")

    assertthat::assert_that(
        length(ref_column_names) == length(query_column_names)
    , msg = "affix_ontology_labels: ref_column_names and query_column_names must be the same length")

    assertthat::assert_that(
        all(ref_column_names %in% colnames(ref_coldata))
    , msg = "affix_ontology_labels: ref_column_names must be a subset of the column names of ref_coldata")

    assertthat::assert_that(
        query_column_name %in% colnames(colData(cds_qry))
    , msg = "affix_ontology_labels: query_column_name must be a column name of cds_qry")

    
    query_cells <- colData(cds_qry)[, query_column_name]
    ref_ontology <- as.data.frame(unique(ref_coldata[, ref_column_names]))

    NUMBER_OF_QUERY_CELLS <- length(query_cells)
    NUMBER_OF_LAYERS <- length(ref_column_names)

    #Turn ref_ontology into a list where the index is the most speficic label and the value is the row
    ref_ontology <- split(ref_ontology, ref_ontology[, NUMBER_OF_LAYERS])
    list_of_full_cell_paths <- list(length=NUMBER_OF_QUERY_CELLS)


    list_of_full_cell_paths <- lapply(query_cells, function(label) {
        if(is.na(label)) return(rep(NA, NUMBER_OF_LAYERS))
        return(as.vector(unlist(ref_ontology[[label]])))
    })

    df_to_bind <- Reduce(rbind, list_of_full_cell_paths)

    colData(cds_qry)[, query_column_names] <- df_to_bind

    return(cds_qry)
}