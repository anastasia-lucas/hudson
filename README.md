# hudson
An R package for creating mirrored Manhattan plots

## Overview
hudson is an R package for creating mirrored Manhattan plots with a shared x-axis, similar to Figure 4 from Verma et al. shown [here](https://www.cell.com/ajhg/fulltext/S0002-9297(18)30062-4). The package includes functions for  You can simply specify a dataset for the top and bottom tracks to generate a basic plot, or provide meta information to annotate a more complex plot.

## Installation
As of now, there is only a development version of the package which can be installed using devtools as follows.

```devtoolls::install_github('anastasia-lucas/hudson')```

## Usage

### Create a mirrored Manhattan plot using GWAS data
```
library(hudson)
#Generate some data

#Create a basic plot
gmirror(top=t, bottom=b)

#Add meta information and text annotations
#Add shape information to our datasets
gmirror(top=t, bottom=b, line=0.0001, annotate_p=0.000001, chrcolor1="", chrcolor2="")
```

### Create a mirrored Manhattan plot using PheWAS data

```
library(hudson)
```

### Create a mirrored Manhattan plot using EWAS data
```
library(hudson)
```
