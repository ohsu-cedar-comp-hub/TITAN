# TITAN (Topic Inferance Of Transcriptionally Associated Networks)

### Summary

TITAN is a tool used to analyze single-cell RNA sequencing data. It is an alternative to the traditional clustering methods for datasets where the differences between groups of cells may not be as distinct. TITAN uses topic modeling to create a number of topics, each with its own unique gene expression profile. The cells are then scored based on how much they express that profile. This is what separates TITAN from traditional clustering. In traditional clustering, the clusters are binary, a cell belongs to one cluster and that cluster only. This topic scoring system allows for a cell to score highly in multiple topics. One can then extract the gene expression profile from a topic to attempt to connect the topic to some sort of transcriptional network. Now you can identify multiple different transcriptional networks that are present in the dataset and determine which cells express that network and how strongly they are expressing it



Tutorial on hormone treatments in breast cancer can be found [here](https://github.com/AlexChitsazan/TITAN/blob/master/vignettes/TITAN_vignette.md)
