a <- function(
    orthology_paths,
    nn_table,
    number_of_labels
) {
  
  
  ratio_of_paths <- matrix(NA, nrow=number_of_labels, ncol=nrow(orthology_paths))
  
  for(j in 1:nrow(orthology_paths)) {
    path <- orthology_paths[j, ]
    ratios <- vector("numeric", length = length(ref_column_names))
    
    for(h in 1:number_of_labels) ratios[h] <- sum(nn_table[[colnames(nn_table)[h]]] == path[[colnames(path)[h]]]) / nrow(nn_table)
    
    ratio_of_paths[, j] <- ratios
  }
  
  return(ratio_of_paths)
  
}

b <- function(
    ref_ontology,
    nn_table,
    number_of_labels
) {
  
  list_of_nn_table_colnames <- colnames(nn_table)
  number_of_paths <- nrow(ref_ontology)
  number_of_k_neighbors <- nrow(nn_table)
  
  ratio_of_paths <- matrix(0, nrow=number_of_labels, ncol=number_of_paths)
  
  #This chunk willcalculate the ratio of paths for unique label in cds_ref
  ratio_of_paths <- sapply(1:number_of_paths, function(j) {
    
    path <- ref_ontology[j, ]
    
    ratios <- sapply(1:number_of_labels, function(h) {
      sum(nn_table[[list_of_nn_table_colnames[h]]] == path[[colnames(path)[h]]]) / number_of_k_neighbors
    })
    
    return(ratios)
  })
  
  return(ratio_of_paths)
}