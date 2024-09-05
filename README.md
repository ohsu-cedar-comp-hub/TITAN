# TITANv2 - in development (Topic Inferance Of Transcriptionally Associated Networks)

### Disclaimer

This is the development branch of TITANv2 which has the capability to incorporate spatial coordinates into the model, creating topics containing both RNA and spatial information. To run TITANv2, first run the create_invDistMat.py script to build the inverse distance matrix. Then run the script spatial_TITANv3.R to create the model. Once the model is created, it can be analyzed using the built-in functions of TITANv1 which can be installed from the main branch.

### Summary

We developed TITAN (Topic Inference of Transcriptionally Associated Networks), an unsupervised Bayesian topic modelling based approach utilizing Latent Dirichlet Allocation (LDA). LDA has the ability to infer sparse data in a bag-of-words model which often contains dropouts and is therefore directly applicable to single cell data. Here, latent transcriptional networks are considered as ‘topics’, with genes as words and cells as documents. TITAN uses an approach similar to cisTopic in that it combines LDA with collapsed Gibbs-sampling to output topics that are linked together by a distribution of scores that relate genes to topics and topics to cells. These outputs can be used to find underlying transcriptional patterns by analyzing gene set contributions to a topic pattern and the expression levels of different topics in each single cell. From a transcriptional perspective, TITAN utilizes the gene-topic and cell-topic scores to visualize underlying transcriptional patterns that can be seen outside the traditional clustering methods, which normally cluster cells into discrete definitions. TITAN distinguishes itself in its ability to assay distinct cell states that can be visualized against traditional clusters and UMAP based projections to infer granular details of scRNA-seq data.



Tutorial on hormone treatments in breast cancer can be found [here](https://github.com/JuliusCampbell/TITAN/blob/master/vignettes/TITAN_vignette.md)
