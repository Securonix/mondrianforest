# mondrianforest

This is the home page of my R implementation of the Mondrian forest algorithms introduced in [Lakshminarayanan at al., 2014](http://www.gatsby.ucl.ac.uk/~balaji/mondrian_forests_nips14.pdf). Mondrian forest is an online alternative to random forest. So far, it only supports classification. Please see the vignette for examples and explanations of its behavior.

Install with
```{r}
# install.packages('devtools')
devtools::install_github('millerjoey/mondrianforest')
```

For a demo, call
```{r}
vignette("mondrianforest")
```
