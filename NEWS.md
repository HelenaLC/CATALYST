# CATALYST 1.4.0

CATALYST now provides a series of visualization tools for differential analysis, most of which have been adapted from Nowicka et al. (2017), and which may be used in conjunction with the `diffcyt` R package (Lukas M. Weber, 2018) to provide a complete framework for differential analysis of CyTOF data.

<blockquote><p><h5>
Nowicka M, Krieg C, Weber LM et al. CyTOF workflow: differential discovery in high-throughput high-dimensional cytometry datasets [version 2; referees: 2 approved]. F1000Research 2017, 6:748 (doi: <a href="10.12688/f1000research.11622.2">10.12688/f1000research.11622.2</a>).
</h5></p><p><h5>
Weber L (2018). diffcyt: Differential discovery in high-dimensional cytometry via high-resolution clustering. R package version 1.0.0, (<a href="https://github.com/lmweber/diffcyt">https://github.com/lmweber/diffcyt</a>).
</h5></p></blockquote>

### The following functions are newly available:

- `plotCounts()`, `plotMDS()`, `plotExprs()`, `plotMedExprs()`: diagnostic plots
- `plotNRS()`: non-redundancy scores
- `tSNE()`, `plotSNE()`: run & visualize t-SNE dimentionality reduction
- `cluster()`, `mergeClusters()`: `FlowSOM` clustering and `ConsensusClusterPlus` metaclustering
- `plotCodes()`: t-SNE and PCA of SOM codes
- `plotAbundances()`: relative population abundances
- `plotClusterHeatmap()`, `plotDiffHeatmap()`: summary of (meta)clustering and differential testing results