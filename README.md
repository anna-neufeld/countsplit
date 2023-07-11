What is countsplit?
-----

The ``countsplit`` R package splits an integer-valued matrix into multiple folds of data with the same dimensions. These folds will be independent under certain modeling assumptions, and can thus be used for cross validation. 

The motivation for this method in the setting where the data are Poisson distributed is described in Neufeld et al., 2022 [(link to paper)](http://arxiv.org/abs/2207.00554) in the context of inference after latent variable estimation for single cell RNA sequencing data. Briefly, count splitting allows users to perform differential expression analysis to see which genes vary across estimated cell types (such as those obtained via clustering) or along an estimated cellular trajectory (pseudotime). Neufeld et al., 2022 [(link to preprint)](???) extends the method to the setting where the data follow a negative binomial distributed, and provides additional settings where count splitting is useful. For example, count splitting is useful broadly for evaluating low-rank representations of the data.  

Recent package updates
-----

We have consolidated the functions in this package such that both Poisson and negative binomial thinning can be performed using the same function; `countsplit`. This function can also be used to create an arbitrary number of folds of data, rather than just a single train/test split. 

The vignettes and data associated with this package are stored in the associated ``countsplit.tutorials" package. To see the tutorials, please visit the updated tutorial website: [https://anna-neufeld.github.io/countsplit.tutorials/](https://anna-neufeld.github.io/countsplit.tutorials/). This change helps with overall package size and build time. *The negative binomial tutorials are still being updated, so most examples on the tutorial website use Poisson thinning only.* 

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

Please visit [https://github.com/anna-neufeld/countsplit_paper](https://github.com/anna-neufeld/countsplit_paper) for code to reproduce the figures and tables from our Poisson paper. 

Please visit [https://github.com/anna-neufeld/nbcs_paper_simulations](https://github.com/anna-neufeld/nbcs_paper_simualations) for code to reproduce the figures and tables from our negative binomial paper. 

References 
----

Neufeld, A.,Gao, L., Popp, J., Battle, A. & Witten, D. (2022), ‘Inference after latent variable estimation for single-cell RNA sequencing data’, arXiv.2207.00554

Neufeld, A.,Dharamshi, A., Gao, L., & Witten, D. (2023), ‘Data thinning for convolution-closed distributions’, https://arxiv.org/abs/2301.07276/ . 

Neufeld, A., Popp, J., Gao, L., Battle, A. & Witten, D. (2023), ‘Negative binomial count splitting for single-cell RNA sequencing data'. 



