What is countsplit?
-----

The ``countsplit`` R package splits an integer-valued matrix into a training matrix and a test matrix using binomial thinning. Under a Poisson assumption, the training and test matrices are independent. 

The motivation for this method is described in [our preprint](XXXXXXX) in the context of inference after latent variable estimation for single cell RNA sequencing data. 

The tutorials on this website are a work in progress, and will continue to be updated over the next few weeks. 

How can I get countsplit?
-----

Make sure that ``remotes`` is installed by running ``install.packages("remotes")``, then type

```{r}
remotes::install_github("anna-neufeld/countsplit")
```

Where can I learn more? 
-----

See the [tutorial](https://anna-neufeld.github.io/countsplit/articles/countsplit_tutorial.html) tab for an introduction to our framework on simple simulated data. See the [scran](https://anna-neufeld.github.io/countsplit/articles/scran_tutorial.html),  [seurat](https://anna-neufeld.github.io/countsplit/articles/seurat_tutorial.html), and [monocle3](https://anna-neufeld.github.io/countsplit/articles/monocle3_tutorial.html) tutorials for examples of how the count splitting package can be integrated with common scRNA-seq analysis pipelines. 

References
----




