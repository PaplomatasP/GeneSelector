# GeneSelector - Advanced Biomarker Discovery Tool

GeneSelector is a sophisticated tool designed to enhance the exploration of biomarkers in single-cell RNA sequencing (scRNA-seq) data. This innovative software combines two crucial aspects of biomarker discovery: Differentially Expressed Genes (DEGs) and Highly Variable Genes (HVGs) within a single framework. GeneSelector leverages the power of Genetic Algorithms (GAs) to perform feature selection, enabling the identification of genes that possess both high variance and biological significance.

## Introduction

GeneSelector simplifies the complex task of biomarker discovery in scRNA-seq data. By merging DEGs, which signify genes with differential expression across different conditions or cell types, and HVGs, which highlight genes with substantial variability within a dataset, GeneSelector ensures a comprehensive analysis of gene expression patterns. This integration allows researchers to pinpoint genes that not only exhibit significant biological relevance but also demonstrate a high degree of variability, making them prime candidates for further investigation as potential biomarkers.

## Usage

GeneSelector is a powerful tool for researchers seeking to uncover genes that play a crucial role in cellular processes and disease mechanisms. With its GA-based feature selection approach, GeneSelector empowers scientists to make informed decisions and advance our understanding of gene expression dynamics in single-cell RNA sequencing studies.

### Case Study: Glioblastoma Multiforme (GBM)

We used GeneSelector to study Glioblastoma multiforme (GBM), a brain cancer known for its aggressive nature and poor prognosis, with distinct subtypes: proneural and mesenchymal. Targeted therapy development hinges on a deep understanding of these subtypes. Reliably pinpointing gene biomarkers specific to each subtype is challenging but critical.

In our study, we present a genetic optimization algorithm tailored to select a precise set of genes capable of clearly differentiating between the proneural and mesenchymal subtypes of GBM. This algorithm innovatively combines differential gene expression and gene variability, employing a dual-criterion strategy. Our approach guarantees the selection of genes that are not just differentially expressed among subtypes but also exhibit a consistent variability pattern, thereby ensuring the biological pertinence of our results. 

Utilizing actual scRNA-seq experimental data, this method aims to discover gene biomarkers that effectively distinguish between the proneural and mesenchymal GBM subtypes. Our findings reveal significant genes with a strong connection to GBM’s critical aspects.

## Results

This study's results highlight the effectiveness of the GeneSelector framework in identifying critical biomarkers in glioblastoma multiforme (GBM) through a novel genetic algorithm (GA)-based feature selection. The integration of differentially expressed genes (DEGs) and highly variable genes (HVGs) within this framework facilitates the discovery of genes that are both biologically significant and exhibit high variance, a necessary step in understanding GBM's complexity.

### Key Observations

- Clear demarcation between the two GBM subtypes (see Fig. 1) indicates the promise of gene-centric analysis.
- Limited overlap between DEGs and HVGs (Figure 2) emphasizes the need for a dual-method approach.
- Integration of DEGs and HVGs within a GA enriches feature selection, yielding robust candidate biomarkers.
- The GA's performance trajectory (Figure 3) demonstrates its adeptness in narrowing down biologically pertinent features.

## Contributing

Contributions to GeneSelector are welcome! If you'd like to contribute to the project or report issues, please refer to our [Contribution Guidelines](CONTRIBUTING.md).

## License

GeneSelector is licensed under the GNU AFFERO GENERAL PUBLIC LICENSE

## Acknowledgments

We would like to thank all contributors to the GeneSelector project for their valuable input and support.
