# LDtools
[![Build Status](https://travis-ci.org/DominikMueller64/LDtools.svg?branch=master)](https://travis-ci.org/DominikMueller64/LDtools)
[![codecov](https://codecov.io/gh/DominikMueller64/LDtools/branch/master/graph/badge.svg)](https://codecov.io/gh/DominikMueller64/LDtools)
[![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![Maintenance](https://img.shields.io/maintenance/yes/2016.svg)]()

---

LDtools provides functions of analysis of linkage disequilibrium and linkage
phases for haplotype (phased) and genotype (unphased) data.

---

### Installation

You can install LDtools from its [GitHub repository](https://github.com/DominikMueller64/LDtools). You first need to install the [devtools](https://github.com/hadley/devtools) package.

```r
install.packages('devtools')
```
Then you can easily install LDtools using the `install_github` function in the [devtools](https://github.com/hadley/devtools) package. (With
`build_vignettes=TRUE`, the vignettes will be built and installed.) 

```R
devtools::install_github("DominikMueller64/LDtools", build_vignettes = TRUE)
```

---

### Vignette

A vignette describing the functionality of the package is available from within R. Load the package and then use the vignette function.

```r
library(LDtools)
vignette('LDtools', package = 'LDtools')
```

---

### Author
Dominik Mueller 

---

### License

This package is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License, version 3, as
published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but
without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the GNU
General Public License for more details.

A copy of the GNU General Public License, version 3, is available at
<http://www.r-project.org/Licenses/GPL-3>




