What is treevalues?
-----

The ``treevalues`` R package computes confidence intervals and p-values for the mean response within a region or the difference in mean response between two regions in a CART regression tree (built using the package ``rpart``). 

Because the regions in a regression tree are selected using the data, we cannot naively "double dip" in the same data to do inference on the means within these regions. 

The ``treevalues`` package implements a selective inference approach to conduct inference without double dipping in the data. 


How can I get treevalues?
-----

Make sure that ``remotes`` is installed by running ``install.packages("remotes")``, then type

```R
remotes::install_github("anna-neufeld/treevalues")
```

Where can I learn more? 
-----

See the [overview](https://anna-neufeld.github.io/treevalues/articles/overview.html) tab for a more detailed motivation for our framework, and the [tutorial](https://anna-neufeld.github.io/treevalues/articles/inference_tutorial.html) tab for instructions on how to use this package on real data. 

See [https://arxiv.org/abs/2106.07816](https://arxiv.org/abs/2106.07816) for the preprint that describes the selective inference methodology. 

See [https://github.com/anna-neufeld/treevalues-simulations](https://github.com/anna-neufeld/treevalues-simulations) for code to reproduce the experiments and figures in the preprint. 




