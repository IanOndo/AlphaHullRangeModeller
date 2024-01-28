
# AlphaHullRangeModeller

`AlphaHullRangeModeller`is an extension of the R package
[*rangeBuilder*](https://github.com/cran/rangeBuilder/tree/master) that
help to generate species range polygons using alpha-hulls, a
generalization of convex hulls that are particularly useful for
estimating species ranges whose habitat is irregularly shaped or where
populations are spatially structured. The package was used to generate
8,702 plant species ranges for the paper [*Areas of global importance
for conserving terrestrial biodiversity, carbon and
water*](https://www.nature.com/articles/s41559-021-01528-7) (Jung et al.
2021) and to help identify the model training area for 33,699 species
distribution models in [*The global distribution of plants used by
humans*](https://www.science.org/doi/10.1126/science.adg8028) (Pironon
and Ondo et al. 2024).   
The parameters for alpha-hulls are adaptively selected, starting with
initial alpha values - a parameter constraining the hull triangulation -
of 2 (or 3) as recommended by the IUCN Red List categories and criteria,
but adjusted for the distribution of records so that at least 95% of the
records were included within the estimated range.

## Enhancements

`AlphaHullRangeModeller` incorporates added features to help improving
and facilitate range estimates. 

- 1)  Variations in the parameter alpha can also affect subpopulation
      structure (i.e. number of subpopulations), so alpha-hulls are
      combined with the “1/10th max” circular buffer method (i.e. the
      buffer size is set to the tenth of the maximum interpoint
      distance) to better capture subpopulation structure (Rivers et al.
      2010).

- 2)  Polygon ranges crossing or intersecting the $180^{th}$ meridian
      (or international date line) are impossible to generate, so
      alpha-hulls are generated using geographically split occurrence
      records, then recombined.

- 3)  The main function, `Make_alpha_hulls`, can be run in parallel to
      help streamlining the generation of hundreds of thousands of
      species range polygons. It is very useful for developing modelling
      pipelines and thus processing large quantities of data
      simultaneously.

## *Installation*

Make sure to have [*R*](https://cloud.r-project.org/ "R") or
[*Rstudio*](https://rstudio.com/products/rstudio/download/ "Rstudio")
installed on your machine. Some R packages need to be compiled from
source, so if you are on Windows, you need to install
[*Rtools*](http://cran.r-project.org/bin/windows/Rtools/) too.  

Install *AlphaHullRangeModeller* with the following instructions:

``` r
devtools::install_github("IanOndo/AlphaHullRangeModeller")
library(AlphaHullRangeModeller)
```

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Jung2021" class="csl-entry">

Jung, Martin, Andy Arnell, Xavier de Lamo, Shaenandhoa García-Rangel,
Matthew Lewis, Jennifer Mark, Cory Merow, et al. 2021. “Areas of Global
Importance for Conserving Terrestrial Biodiversity, Carbon and Water.”
Journal Article. *Nature Ecology & Evolution* 5 (11): 1499–509.
<https://doi.org/10.1038/s41559-021-01528-7>.

</div>

<div id="ref-UsefulPlants" class="csl-entry">

Pironon and Ondo, M. Diazgranados, R. Allkin, A. C. Baquero, R.
Cámara-Leret, C. Canteiro, Z. Dennehy-Carr, et al. 2024. “The Global
Distribution of Plants Used by Humans.” *Science* 383 (6680): 293–97.
<https://doi.org/10.1126/science.adg8028>.

</div>

<div id="ref-Rivers2010SubpopulationsLA" class="csl-entry">

Rivers, Malin C., Steven P. Bachman, Thomas R. Meagher, Eimear M Nic
Lughadha, and Neil Brummitt. 2010. “Subpopulations, Locations and
Fragmentation: Applying IUCN Red List Criteria to Herbarium Specimen
Data.” *Biodiversity and Conservation* 19: 2071–85.
<https://api.semanticscholar.org/CorpusID:26800806>.

</div>

</div>
