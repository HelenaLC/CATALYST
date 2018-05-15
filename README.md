### Welcome to CATALYST <img src="docs/logo.png" width="200" align="right" />

Mass cytometry (CyTOF) uses heavy metal isotopes rather than fluorescent tags as reporters to label antibodies, thereby substantially decreasing spectral overlap and allowing for examination of over 50 parameters at the single cell level. While spectral overlap is significantly less pronounced in CyTOF than flow cytometry, spillover due to detection sensitivity, isotopic impurities, and oxide formation can impede data interpretability. 

We designed CATALYST (Cytometry dATa anALYSis Tools) to provide:

- easy-to-use visualizations for differential analysis
- a pipeline for preprocessing of cytometry data, including
  - normalization using bead standards
  - single-cell deconvolution
  - bead-based compensation

[![platforms](http://bioconductor.org/shields/availability/3.7/CATALYST.svg)](http://bioconductor.org/packages/release/bioc/html/CATALYST.html#archives)&nbsp;
[![downloads](http://bioconductor.org/shields/downloads/CATALYST.svg)](http://bioconductor.org/packages/stats/bioc/CATALYST/)&nbsp;
[![posts](http://bioconductor.org/shields/posts/CATALYST.svg)](https://support.bioconductor.org/t/catalyst)&nbsp;
[![in Bioc](http://bioconductor.org/shields/years-in-bioc/CATALYST.svg)](http://bioconductor.org/packages/release/bioc/html/CATALYST.html#since)&nbsp;
[![build](http://bioconductor.org/shields/build/release/bioc/CATALYST.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/CATALYST)

### Quick links

* The preprint to this project is accesible on [BioRxiv preprint](https://doi.org/10.1101/185744).
* Our manuscript is available at [Cell Systems](https://doi.org/10.1016/j.cels.2018.02.010).
* Scripts to reproduce a majority of the figures from our paper can be found [here](https://github.com/BodenmillerGroup/cyTOFcompensation).
* To use our web application, go to [CATALYSTLite](http://imlspenticton.uzh.ch:3838/CATALYSTLite).

---

### Getting started

There are various entry points to using the CATALYST tools. If you are an R user, we recommend using the functions from the R package itself, from the command line or via a script. If you want to use the graphical environment, we recommend to run it locally on your own computer, as this will be faster and in terms of dataset sizes, is only limited by the resources on your computer. To do this, you will need to install recent versions of R and the necessary R packages (details below). If installing all the necessary R packages is a hurdle, you can use the hosted version.

- **Using the R package**
  - A stable release version is available at [Bioconductor](http://bioconductor.org/packages/CATALYST).
  - For detailed examples and usage instructions, see the package vignettes under **Articels**.
  
- **Using the GUI**
  - To use the Shiny app locally, run `CATALYST::launchGUI()` inside an R session.
  - Go to [CATALYSTLite](http://imlspenticton.uzh.ch:3838/CATALYSTLite) to use the app online.
  
---

### Submitting an issue or feature request

`CATALYST` is still under active development. We greatly welcome (and highly encourage!) all feedback, bug reports and suggestions for improvement. **Please make sure to raise issues with a reproducible example and, if you're running R, the output of your `sessionInfo()`.**

* R package and Shiny related issues should be raised [here](https://github.com/HelenaLC/CATALYST/issues).
* For general questions and feedback, please contact us directly via email.

---

### Contact

For software related issues:
* Helena L Crowell (crowellh@student.ethz.ch)
* Mark D Robinson (mark.robinson@imls.uzh.ch)

Regarding methodological questions:
* Stéphane Chevrier (stephane.chevrier@uzh.ch)
* Vito Zanotelli (vito.zanotelli@uzh.ch)