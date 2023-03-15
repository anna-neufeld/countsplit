What is countsplit?
-----

The ``countsplit`` R package splits an integer-valued matrix into a training matrix and a test matrix using binomial thinning. Under a Poisson assumption, the training and test matrices are independent. 

The motivation for this method is described in Neufeld et al., 2022 [(link to paper)](http://arxiv.org/abs/2207.00554) in the context of inference after latent variable estimation for single cell RNA sequencing data. Briefly, count splitting allows users to perform differential expression analysis to see which genes vary across estimated cell types (such as those obtained via clustering) or along an estimated cellular trajectory (pseudotime). 

Recent package updates
-----

We have improved the ability of the package to work with sparse matrices. We have added a negative binomial count splitting function, but the tutorials for this function are still a work in progress. This function implements the decomposition of the negative binomial described in Neufeld et al., 2022 [(link to preprint)](https://arxiv.org/abs/2301.07276). 

The vignettes and data associated with this package are stored in the associated ``countsplit.tutorials" package. To see the tutorials, please visit the updated tutorial website: [https://anna-neufeld.github.io/countsplit.tutorials/](https://anna-neufeld.github.io/countsplit.tutorials/). This change helps with overall package size and build time. 

How can I get countsplit?
-----

Make sure that ``remotes`` is installed by running ``install.packages("remotes")``, then type

```{r}
remotes::install_github("anna-neufeld/countsplit")
```

To also download the data needed to reproduce the package vignettes, be sure to also install the ``countsplit.tutorials" package.

```{r}
remotes::install_github("anna-neufeld/countsplit.tutorials"). 
```

Where can I learn more? 
-----

Please visit our tutorial website [https://anna-neufeld.github.io/countsplit.tutorials/](https://anna-neufeld.github.io/countsplit.tutorials/) to see an introduction to our framework on simple simulated data, as well as tutorials for integrating the count splitting package with common scRNA-seq analysis pipelines (Seurat, scran, and Monocle3). 

Please visit [https://github.com/anna-neufeld/countsplit_paper](https://github.com/anna-neufeld/countsplit_paper) for code to reproduce the figures and tables from our paper. 

References 
----

Neufeld, A.,Gao, L., Popp, J., Battle, A. & Witten, D. (2022), ‘Inference after latent variable estimation for single-cell RNA sequencing data’, arXiv.2207.00554

Neufeld, A.,Dharamshi, A., Gao, L., & Witten, D. (2023), ‘Data thinning for convolution-closed distributions’, https://arxiv.org/abs/2301.07276/ . 



