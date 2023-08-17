# scRNA-seq Bayesian Label Transfer

During my summer 2023 internship at the Trapnell Lab, I developed a novel method to transfer cell labels from a reference [monocle3](https://github.com/cole-trapnell-lab/monocle3) cell dataset to a query cell dataset while respecting ontogeny.

Single-cell RNA sequencing (scRNA-seq) assays are powerful tools for analyzing the gene expression profiles of individual cells. However, they often produce extensive datasets with numerous cells, making manual cell labeling a labor-intensive task. **Label transfer algorithms** are employed to rapidly annotate cells in scRNA-seq datasets. However, to the best of my knowledge, none of the existing algorithms incorporate cell ontogeny information. Moreover, these methods enforce strict ontogeny constraints for cell labeling, which may not accurately represent the true biological context.

The primary distinction between the label transfer technique in this repository and others lies in the integration of **priors** for each layer of a cell ontogeny. Introducing these priors ensures that this information is taken into account during the classification of a query cell and allows for the transfer of cells that may not conform to the strict ontogeny.

## Benchmark
In a paper by Weinreb and Klein (2020), the group developed a computational tool, CLiNC, to learn cell fate choices, or a cell ontogeny, from single-cell clonal barcoding. In their paper, they utilized CLiNC on a published hematopoiesis scRNAseq dataset provided by Weinreb, Rodriguez-Fraticelli, Camargo, and Klein (2020).

CLiNC, not only reconstructs a cell ontogeny but also detected cross-tree transition:
<img width="400" alt="image" src="https://github.com/sjsikora/bayesian-scRNAseq-label-transfer/assets/20007305/d012c3cc-42fb-456a-a9f6-7bc1bfed9e98">


To assess this algorithm's ability to transfer cross-tree transitions in a query, the same scRNAseq dataset was imported into monocle3. Undifferentiated and progenitor cells were removed. The dataset was split by libraries into a reference and query data set. Since the scRNAseq only annotated the cells by cell type, a cell ontogeny was applied to the reference cell dataset with cross-tree transitions simulated by random chance (30%). Finally, the reference and the query were run through the main function. 

<details>
    <summary>Reference cell data set plotted by measured cell type:</summary>
        <img width="600" alt="cds_ref_cell_type" src="https://github.com/sjsikora/bayesian-scRNAseq-label-transfer/assets/20007305/f4105040-2ec5-4016-b956-d10260a85748">

</details>
<details>
    <summary>Reference cell data set plotted by simulated cross-tree transitions:</summary>
        <img width="600" alt="cds_ref_breakage" src="https://github.com/sjsikora/bayesian-scRNAseq-label-transfer/assets/20007305/ad049d85-c5b8-4d11-9215-b9f8e1336dbf">
</details>
<details>
    <summary>Query cell data set plotted by predicted cell type:</summary>
        <img width="600" alt="cds_qry_cell_type" src="https://github.com/sjsikora/bayesian-scRNAseq-label-transfer/assets/20007305/2e98fa0f-167c-488a-a86f-24839a9ce90c">
</details>
<details>
    <summary>Query cell data set plotted by predicted cross-tree transitions:</summary>
        <img width="600" alt="cds_qry_breakage" src="https://github.com/sjsikora/bayesian-scRNAseq-label-transfer/assets/20007305/f76754dc-4b7b-4163-a2e2-bab0187c8083">
</details>

Results:
- The algorithm was able to correctly identify 99.12% of the query cell types.
- It detected 215 / 4766 Mono cell and 9 / 1188 Baso cell cross-tree transitions.

## Algorithm

- The reference and query cell datasets are reduced into a coordinate space.
- A k-nearest-neighbor framework is established to compare both the reference and query datasets to the reference dataset.
- Priors are initialized and trained on the reference dataset:
    - For every reference cell, retrieve the k-nearest neighbors of that cell among k other reference cells.
    - Transform the k-nearest neighbor table into an ontogeny matrix, where each column represents a path within the cell ontogeny, each row represents a cell, and each entry indicates the proportion of k-nearest neighbors that report that label divided by k.
    - To train the priors, calculate the loss of a given prior and refine accordingly.
- Once the priors are optimized on the reference cell dataset, retrieve the k-nearest neighbors of the query cell dataset with respect to the reference cell dataset.
- Transform the k-nearest neighbor table into an ontogeny matrix.
- Multiply the priors matrix by the ontogeny matrix and report the path with the highest value.

## Citations and Acknowledgements

monocle3:

Trapnell, C., Cacchiarelli, D., Grimsby, J. et al. The dynamics and regulators of cell fate decisions are revealed by pseudotemporal ordering of single cells. Nat Biotechnol 32, 381–386 (2014). https://doi.org/10.1038/nbt.2859

Qiu, X., Mao, Q., Tang, Y. et al. Reversed graph embedding resolves complex single-cell trajectories. Nat Methods 14, 979–982 (2017). https://doi.org/10.1038/nmeth.4402

Cao, J., Spielmann, M., Qiu, X. et al. The single-cell transcriptional landscape of mammalian organogenesis. Nature 566, 496–502 (2019). https://doi.org/10.1038/s41586-019-0969-x

UMAP:

McInnes et al., (2018). UMAP: Uniform Manifold Approximation and Projection. Journal of Open Source Software, 3(29), 861, https://doi.org/10.21105/joss.00861

Tutorial Data:

Weinreb, C., Rodriguez-Fraticelli, A., Camargo, F. D., & Klein, A. M. (2020). Lineage tracing on transcriptional landscapes links state to fate during differentiation. Science (New York, N.Y.), 367(6479), eaaw3381. https://doi.org/10.1126/science.aaw3381

Weinreb, C., & Klein, A. M. (2020). Lineage reconstruction from clonal correlations. Proceedings of the National Academy of Sciences - PNAS, 117(29), 17041-17048. https://doi.org/10.1073/pnas.2000238117
