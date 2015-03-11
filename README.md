# QuickShiftClustering

[![Build Status](https://travis-ci.org/rened/QuickShiftClustering.jl.svg?branch=master)](https://travis-ci.org/rened/QuickShiftClustering.jl)

QuickShift [1] is a fast method for hierarchical clustering, which first constructs the clustering tree, and subsequently allows to quickly cut links in the tree which exceed a specified length. This second step can be performed for different link-lengths without having to re-run the clustering itself.

[1] [Quick Shift and Kernel Methods
for Mode Seeking](http://cronos.rutgers.edu/~meer/TEACH/ADD/vedaldiS08quick.pdf)

#### Example

```jl
using FunctionalData
data = @p map unstack(1:10) (x->10*randn(2,1).+randn(2,100)) | flatten
# data is nDim x nSamples, i.e. 2 x 1000

using QuickShiftClustering
a = quickshift(data)           # optional second parameter is sigma (refer to paper)
labels = quickshiftlabels(a)   # optional second param is max link length

quickshiftplot(a, data, labels)
```

![](example.png)
