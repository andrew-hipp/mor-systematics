---
title: Changing tip labels using morTipLabels
author: Andrew Hipp
date: 2021-01-28
output: pdf_document
---

The function `morTipLabels` is a pretty rudimentary method for relabeling tips for our own lab data, mainly the trees generated from RAD-seq data... but it's generalizable. To run it, first make sure you have installed `ape`, `magrittr`, and `openxlsx`. Then, you'll need to grab the function off of GitHub:


```r
source('https://raw.githubusercontent.com/andrew-hipp/mor-systematics/master/AAA.LAB/morTipLabels.R')
```

This reads the function into `R`. It also assumes you might want to use some of the specimen metadata files from the T: drive, so it _tries_ to read those in; it should work fine even if it doesn't succeed. To see whether it has succeeded, look at the contents of your workspace:


```r
ls()
```

```
## [1] "a"            "dat.acer"     "dat.malus"    "dat.tilia"    "morTipLabels" "tidyName"     "tr"
```

You should see `dat.acer`, `dat.malus`, and `dat.tilia` in your workspace. If you do, then hooray! If not, then the T: drive was probably not mounted at the time that you sourced in the file.
