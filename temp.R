temp <- function() {

    list_of_priors <- list()


    for(i in 1:10) {

        priors <- c(0.25, 0.25, 0.25, 0.25)

        priors <- train_priors_on_reference(
            priors, 
            query_search, 
            ref_coldata, 
            ref_column_names, 
            ref_ontology,
            NUMBER_OF_REFERENCE_CELLS, 
            NUMBER_OF_LABELS
        )

        list_of_priors[[i]] <- priors
    }

    return(list_of_priors)
}