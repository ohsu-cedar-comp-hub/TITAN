# TITAN (Topic Inferance Of Transcriptionally Associated Networks)

### Summary

TITAN (Topic Inference of Transcriptionally Associated Networks is an unsupervised Bayesian topic modelling based approach utilizing LDA (Latent Dirichlet Allocation). In this approach, latent transcriptional networks are considered as ‘topics’, linking gene patterns to cells. TITAN uses a similar approach to cisTopic in that it combines LDA with collapsed Gibbs-sampling to output topics that are linked together by a distribution of scores that relate a) genes to topics and b) topics to cells. These outputs can be used to find underlying transcriptional patterns by analysing gene network contributions to a topic pattern and the expression levels of different topics in each single cell.

In this sense, it is an alternative to traditional scRNA-seq analysis methods for datasets that are herterogenious or where the differences between groups of cells may not be as distinct. TITAN uses topic modeling to create n number of topics, each with its own unique gene expression profile. The cells are then scored based on how much they express that profile. One way TITAN seperates itself from traditional clustering is that the topic-cell association is not descrete, the topic scoring system allows for a cell to score highly in multiple topics based upon the gene expression profile. One can then extract the gene expression profile from a topic to attempt to connect the topic to a known transcriptional network. TITAN allows you to identify multiple different transcriptional networks that are present in the dataset and determine which cells express that network and how strongly they are expressing it.



Tutorial on hormone treatments in breast cancer can be found [here](https://github.com/JuliusCampbell/TITAN/blob/master/vignettes/TITAN_vignette.md)
