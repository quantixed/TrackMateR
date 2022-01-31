# TrackMateR

Analysis of TrackMate XML outputs in R

[TrackMate](https://imagej.net/plugins/trackmate/) is a single-particle tracking plugin for ImageJ/Fiji [[1]](#1)[[2]](#2). The standard output from a tracking session is in TrackMate XML format.

The goal of this R package is to import all of the data associated with the final filtered tracks for further analaysis and visualization in R. This package builds on initial work by Julien Godet on [trackR](https://github.com/jgodet/trackR).

## References
<a id="1">[1]</a> 
Tinevez, J.-Y. et al. (2017)
TrackMate: An open and extensible platform for single-particle tracking.
Methods, 115, 80–90. doi:10.1016/j.ymeth.2016.09.016

<a id="2">[2]</a>
Ershov, D., et al. (2021)
Bringing TrackMate into the era of machine-learning and deep-learning.
bioRxiv. doi:10.1101/2021.09.03.458852