---
editor_options: 
  markdown: 
    wrap: 72
---

# Install this forked version

```         
if (!require(devtools)) {
    install.packages('devtools')
}

devtools::install_github('eugloh/drawProteins')
```

[![Build
Status](https://travis-ci.org/brennanpincardiff/drawProteins.svg?branch=master)](https://travis-ci.org/brennanpincardiff/drawProteins)

[![Coverage
Status](https://coveralls.io/repos/github/brennanpincardiff/drawProteins/badge.svg?branch=master)](https://coveralls.io/github/brennanpincardiff/drawProteins?branch=master)

# Overview of drawProteins

This package has been created to allow the visualisation of protein
schematics based on the data obtained from the [Uniprot Protein
Database](http://www.uniprot.org/).

<img src="https://raw.githubusercontent.com/Bioconductor/BiocStickers/master/drawProteins/drawProteins.png" height="200" align="right"/>

## The basic workflow is:

-   to provide one or more Uniprot IDs

-   get a list of features from the Uniprot API

-   draw the basic chains of these proteins

-   add features as desired

drawProteins uses the package `httr` to interact with the Uniprot API
and extract a JSON object into R. The JSON object is used to create a
data.table.

The graphing package `ggplot2` is then used to create the protein
schematic.

The Vignette gives a good overview of the package.

Sample script on [R for Biochemists
blog](http://rforbiochemists.blogspot.co.uk/2017/11/using-drawproteins-for-draw-nf-kappab.html)

The `master` version of this package is available through
[Bioconductor](http://bioconductor.org/packages/release/bioc/html/drawProteins.html).
The `development` version is currently available here through GitHub.

# To install and use development branch:

This package is available through
[Bioconductor](https://bioconductor.org/packages/devel/bioc/html/drawProteins.html).
Installation instructions and documentation are
[available](https://bioconductor.org/packages/devel/bioc/html/drawProteins.html).

To install from Github:

```         
if (!require(devtools)) {
    install.packages('devtools')
}
dev_mode(on=TRUE)
devtools::install_github('brennanpincardiff/drawProteins')
```

------------------------------------------------------------------------

Please cite the following article when using `drawProteins`:

**Brennan P**. drawProteins: a Bioconductor/R package for reproducible
and programmatic generation of protein schematics [version 1; referees:
2 approved]. ***F1000Research*** 2018, 7:1105

------------------------------------------------------------------------

Feedback is <b>very</b> welcome. Please raise [Github
issues](https://github.com/brennanpincardiff/drawProteins/issues) to
provide bug reports, give feedback or request features.
