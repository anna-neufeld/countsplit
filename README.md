What is countsplit?
-----

The ``countsplit`` R package splits an integer-valued matrix into a training matrix and a test matrix using binomial thinning. Under a Poisson assumption, the training and test matrices are independent. 

The motivation for this method is described in Neufeld et al., 2022 [(link to preprint)](http://arxiv.org/abs/2207.00554) in the context of inference after latent variable estimation for single cell RNA sequencing data. Briefly, count splitting allows users to perform differential expression analysis to see which genes vary across estimated cell types (such as those obtained via clustering) or along an estimated cellular trajectory (pseudotime). 


How can I get countsplit?
-----

Make sure that ``remotes`` is installed by running ``install.packages("remotes")``, then type

```{r}
remotes::install_github("anna-neufeld/countsplit")
```

Where can I learn more? 
-----

See the [introductory tutorial](https://anna-neufeld.github.io/countsplit/articles/countsplit_tutorial.html) tab for an introduction to our framework on simple simulated data. See the [seurat](https://anna-neufeld.github.io/countsplit/articles/seurat_tutorial.html),
[scran](https://anna-neufeld.github.io/countsplit/articles/scran_tutorial.html), and [monocle3](https://anna-neufeld.github.io/countsplit/articles/monocle3_tutorial.html) tutorials for examples of how the count splitting package can be integrated with common scRNA-seq analysis pipelines. 

Please visit [https://github.com/anna-neufeld/countsplit_paper](https://github.com/anna-neufeld/countsplit_paper) for code to reproduce the figures and tables from our paper. 

References 
----

Neufeld, A.,Gao, L., Popp, J., Battle, A. & Witten, D. (2022), ‘Inference after latent variable estimation for single-cell RNA sequencing data’, arXiv.2207.00554



