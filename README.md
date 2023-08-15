# scRNA-seq Bayesian Label Transfer

As a part of my summer 2023 intership at the Trapnell Lab, I developed a novel way to transfer cell-labels from a reference monocle3 cell data set to a query cell data set respecting ontology.

Single-cell RNA sequencing (scRNA-seq) assays are powerful tools for analyzing the gene expression space of single cells. However, they can generate large datasets with many cells, making it laborious to label each cell manually. **Label transfer algorithms** are used to quickly annotate cells in scRNA-seq datasets. Yet, to my knowledge, none of the current algorithms incorporate cell-ontology information. In addition, these label transfer enforce a strict ontology to label cells, which in some cases does not reflect real biology.

The key difference between this repo's label transfer and all the others is the addition of **priors** for each layer of a cell ontology. The introduction of these priors ensure that this information is inculded in the classification of a query cell and have the capiabilty of transfering cells that "break" the strict ontology.

## Algorithm

First, the reference and query cell data sets are reduced into coordinate space. 