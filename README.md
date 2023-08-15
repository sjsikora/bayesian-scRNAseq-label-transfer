# scRNA-seq Bayesian Label Transfer

As a part of my summer 2023 intership at the Trapnell Lab, I developed a novel way to transfer cell-labels from a reference monocle3 cell data set to a query cell data set respecting ontology.

Single-cell RNA sequencing (scRNA-seq) assays are powerful tools for analyzing the gene expression space of single cells. However, they can generate large datasets with many cells, making it laborious to label each cell manually. **Label transfer algorithms** are used to quickly annotate cells in scRNA-seq datasets. Yet, to my knowledge, none of the current algorithms incorporate cell-ontology information. In addition, these label transfer enforce a strict ontology to label cells, which in some cases does not reflect real biology.

The key difference between this repo's label transfer and all the others is the addition of **priors** for each layer of a cell ontology. The introduction of these priors ensure that this information is inculded in the classification of a query cell and have the capiabilty of transfering cells that "break" the strict ontology.

## Algorithm

- The reference and query cell data sets are reduced into coordinate space. 
- A k-nearest-neighbor framework is built to compare both reference and query data sets to reference data set. 
- Priors are initialized and trained on reference:
    - For every reference cell, retrieve k-NN of that reference cell to k other reference cells.
    - Transform k-NN table to a ontology matrix where every column is a path down cell ontology, every row is a cell, and every entry is the (number of k-NN reporting that label)/k
    - To train the priors, calculate the hinge-loss* of a given prior and refine.
- Once priors are optimized on the reference cell data set, retrieve k-NN of the query cell data to reference cell data set.
- Transform k-NN table to a ontology matrix.
- Multiply priors and ontology matrix and report back the maxium path.

## *Note on Hinge Loss:

A big point of diffculty in this project was the loss function. I settled on hinge loss because of its demphasis on increasing the margin of correct anwsers. While this may seem counter intutive, the result of methods that emphasize increasing margins are priors of 0 0 0 1.  