### CATALYST - Cytometry dATa anALYsis Tools

To install the package from github, run the following code:

```r
library(devtools)
devtools::install_github(repo = "HelenaLC/CATALYST", 
                         auth_token = "65a3354f8323e33574e9659dfc9f639a47149e47")
```

To launch the Shiny GUI, run:

```r
CATALYST::launchGUI()
```

**02/27/17 update**
- *** added Shiny GUI ***
- fixed argument "method" in `computeSpillmat()`

**02/18/17 update**

- added validity check for "n_events" in `plotEvents()` and  
  warning about populations with insufficient event assignments

**02/18/17 update**

- counts and yields are now computed in `assignPrelim()`:  
  `estCutoffs()` may be skipped if `sep_cutoffs` are supplied manually
- added options "annotate" and "palette" in `plotSpillmat()`
- added validity check for "which" in `plotEvents()` and `plotYields()`

**02/17/17 update**

- *** added `plotMahal()` ***
- fixed bug when calculating yields
- fixed random sampling of colors when calling `plotEvents()`

**02/12/17 update**

- fixed bugs issued by Vito
- fixed titles in `plotYields()` and `plotEvents()`
- decreased running time by avoiding assigning IDs again after normalization
