# transform.hazards

```R``` package for transforming cumulateve hazard estimates from Aalens additive hazard model<sup>[1](#transforming)</sup>.

## Getting Started


### Installing

Make sure to have the ``devtools`` package installed. Then run

```
devtools::install_github("palryalen/transform.hazards", build_vignettes=TRUE).
```
You can also download the source from this page and build the package manually.

### Vignette
Please find a detailed vignette with worked examples on how to utilize the package using the command ``` devtools::browseVignettes("transform.hazards") ```

## Examples

```
devtools::browseVignettes("transform.hazards")

library(transform.hazards)

example(pluginEstimate)
```

## Versioning

We use [SemVer](http://semver.org/) for versioning.

## Authors

* **Pål Ryalen**

## License

This project is licensed under the Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0) licence. You are therefore [free to](https://creativecommons.org/licenses/by-nc/4.0/)

* Share — copy and redistribute the material in any medium or format.

* Adapt — remix, transform, and build upon the material.


## References




<a name="transforming">1</a>: Pål Ryalen, Mats Stensrud, Kjetil Røysland, [*Transforming cumulative hazard estimates*](https://arxiv.org/abs/1710.07422v3)

<!---
<a name="additive_consistent">2</a>: Pål Ryalen, Mats Stensrud, Kjetil Røysland, [*The additive hazard estimator is consistent for continuous time marginal structural models*](https://arxiv.org/abs/1802.01946)
-->

