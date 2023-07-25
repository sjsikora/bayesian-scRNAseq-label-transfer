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
    ref_ontology <- unique(ref_coldata[, ref_column_names])

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

affix_ontology_labels(cds_qry, colData(cds_ref), ref_column_names, query_column_name, c("G", "T", "B", "S"))