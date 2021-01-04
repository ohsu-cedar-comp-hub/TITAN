# TITAN (Topic Inferance Of Transcriptionally Associated Networks)

### Summary

We developed TITAN (Topic Inference of Transcriptionally Associated Networks), an unsupervised Bayesian topic modelling based approach utilizing Latent Dirichlet Allocation (LDA). LDA has the ability to infer sparse data in a bag-of-words model which often contains dropouts and is therefore directly applicable to single cell data. Here, latent transcriptional networks are considered as ‘topics’, with genes as words and cells as documents. TITAN uses an approach similar to cisTopic in that it combines LDA with collapsed Gibbs-sampling to output topics that are linked together by a distribution of scores that relate genes to topics and topics to cells. These outputs can be used to find underlying transcriptional patterns by analyzing gene set contributions to a topic pattern and the expression levels of different topics in each single cell. From a transcriptional perspective, TITAN utilizes the gene-topic and cell-topic scores to visualize underlying transcriptional patterns that can be seen outside the traditional clustering methods, which normally cluster cells into discrete definitions. TITAN distinguishes itself in its ability to assay distinct cell states that can be visualized against traditional clusters and UMAP based projections to infer granular details of scRNA-seq data.



Tutorial on hormone treatments in breast cancer can be found [here](https://github.com/JuliusCampbell/TITAN/blob/master/vignettes/TITAN_vignette.md)
