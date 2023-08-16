# scRNA-seq Bayesian Label Transfer

During my summer 2023 internship at the Trapnell Lab, I developed a novel method to transfer cell labels from a reference monocle3 cell dataset to a query cell dataset while respecting ontogeny.

Single-cell RNA sequencing (scRNA-seq) assays are powerful tools for analyzing the gene expression profiles of individual cells. However, they often produce extensive datasets with numerous cells, making manual cell labeling a labor-intensive task. **Label transfer algorithms** are employed to rapidly annotate cells in scRNA-seq datasets. However, to the best of my knowledge, none of the existing algorithms incorporate cell ontogeny information. Furthermore, these label transfer methods enforce a strict ontogeny for cell labeling, which may not always accurately represent actual biological.

The primary distinction between the label transfer technique in this repository and others lies in the integration of **priors** for each layer of a cell ontogeny. Introducing these priors ensures that this information is taken into account during the classification of a query cell and allows for the transfer of cells that may not conform to the strict ontogeny.

## Benchmark

## Algorithm

- The reference and query cell datasets are reduced into a coordinate space.
- A k-nearest-neighbor framework is established to compare both the reference and query datasets to the reference dataset.
- Priors are initialized and trained on the reference dataset:
    - For every reference cell, retrieve the k-nearest neighbors of that cell among k other reference cells.
    - Transform the k-nearest neighbor table into an ontogeny matrix, where each column represents a path within the cell ontogeny, each row represents a cell, and each entry indicates the proportion of k-nearest neighbors that report that label divided by k.
    - To train the priors, calculate the hinge loss* of a given prior and refine accordingly.
- Once the priors are optimized on the reference cell dataset, retrieve the k-nearest neighbors of the query cell dataset with respect to the reference cell dataset.
- Transform the k-nearest neighbor table into an ontogeny matrix.
- Multiply the priors matrix by the ontogeny matrix and report the path with the highest value.

## *Note on Hinge Loss:

A major challenge in this project was designing an effective loss function. I settled on the hinge loss due to its emphasis on increasing the margin between correct answers. While this might seem counterintuitive, methods that prioritize increasing margins often result in priors like [0, 0, 0, 1].

## Citations and Acknowledgements

Monocle3:

Trapnell, C., Cacchiarelli, D., Grimsby, J. et al. The dynamics and regulators of cell fate decisions are revealed by pseudotemporal ordering of single cells. Nat Biotechnol 32, 381–386 (2014). https://doi.org/10.1038/nbt.2859

Qiu, X., Mao, Q., Tang, Y. et al. Reversed graph embedding resolves complex single-cell trajectories. Nat Methods 14, 979–982 (2017). https://doi.org/10.1038/nmeth.4402

Cao, J., Spielmann, M., Qiu, X. et al. The single-cell transcriptional landscape of mammalian organogenesis. Nature 566, 496–502 (2019). https://doi.org/10.1038/s41586-019-0969-x

UMAP:

McInnes et al., (2018). UMAP: Uniform Manifold Approximation and Projection. Journal of Open Source Software, 3(29), 861, https://doi.org/10.21105/joss.00861

Tutorial Data:

Weinreb, C., Rodriguez-Fraticelli, A., Camargo, F. D., & Klein, A. M. (2020). Lineage tracing on transcriptional landscapes links state to fate during differentiation. Science (New York, N.Y.), 367(6479), eaaw3381. https://doi.org/10.1126/science.aaw3381

Weinreb, C., & Klein, A. M. (2020). Lineage reconstruction from clonal correlations. Proceedings of the National Academy of Sciences - PNAS, 117(29), 17041-17048. https://doi.org/10.1073/pnas.2000238117