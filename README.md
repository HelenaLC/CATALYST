### CATALYST - Cytometry dATa anALYsis Tools

To install the package from github, run the following code:

```r
library(devtools)
devtools::install_github(repo = "HelenaLC/CATALYST", 
                         auth_token = "65a3354f8323e33574e9659dfc9f639a47149e47")
```

**up next**

- compute counts and yields in `assignPrelim()` so that 
  `estCutoffs()` can be skipped if `sep_cutoffs` are supplied manually
- make `plotYields()` work if `sep_cutoffs` slot is empty
- check necessary slots are filled when running functions
- inform user if yield is below 10% in `estCutoffs()` and `applyCutoffs()`
  
**02/17/17 update**

- *** added `plotMahal()` ***
- fixed bug when calculating yields
- fixed random sampling of colors when calling `plotEvents()`

**02/12/17 update**

- fixed bugs issued by Vito
- fixed titles in `plotYields()` and `plotEvents()`
- decreased running time by avoiding assigning IDs again after normalization
