# hudson
An R package for creating mirrored Manhattan plots

## Overview
hudson is an R package for creating mirrored Manhattan plots with a shared x-axis, similar to Figure 4 from Verma et al. shown [here](https://www.cell.com/ajhg/fulltext/S0002-9297(18)30062-4) for position by position comparison of results. The package includes functions to visualize data from genome-wide, phenome-wide, and environment-wide association analyses (GWAS, PheWAS, EWAS, respectively) directly, though they may adaptable for other types of data such as beta or SNP intensity value, or other types of analyses. You can simply specify a dataset for the top and bottom tracks to generate a basic plot, or provide meta information to annotate a more complex plot.

## Installation
As of now, there is only a development version of the package which can be installed using devtools.

```devtools::install_github('anastasia-lucas/hudson')```

This package uses ggplot2, gridExtra, and RColorBrewer. ggrepel is suggested for improved text annotation, but not required.

## Usage

### Create a mirrored Manhattan plot using GWAS data
```
library(hudson)
#Generate some data
t <- data.frame(SNP=paste("rs", seq(1:5000), sep=""),
                CHR=rep(c(1:22, "X", "Y"), length.out=5000, each=200),
                POS=rep(seq(1, 10000, by = 200), length.out=5000),
                pvalue=runif(n=5000))
b <- data.frame(SNP=paste("rs", seq(1:5000), sep=""),
                CHR=rep(c(1:22, "X", "Y"), length.out=5000, each=200),
                POS=rep(seq(1, 10000, by = 200), length.out=5000),
                pvalue=runif(n=5000))
head(t)
  SNP CHR  POS    pvalue
1 rs1   1    1 0.7489581
2 rs2   1  201 0.6775381
3 rs3   1  401 0.4838361
4 rs4   1  601 0.1743456
5 rs5   1  801 0.1880735
6 rs6   1 1001 0.6829087
head(b)
  SNP CHR  POS    pvalue
1 rs1   1    1 0.3163273
2 rs2   1  201 0.5892635
3 rs3   1  401 0.7978299
4 rs4   1  601 0.7216598
5 rs5   1  801 0.7388363
6 rs6   1 1001 0.8230513

#Create a basic plot
gmirror(top=t, bottom=b)

```

![Imgur](https://i.imgur.com/5wYjMzJ.png)

```
#Add meta information and text annotations
#Add shape information to our datasets
t$Shape <- rep(paste("S", seq(1:3), sep="") , length.out=5000, each=1)
b$Shape <- rep(paste("S", seq(1:3), sep="") , length.out=5000, each=1)

#Create more informative plot
gmirror(top=t, bottom=b, tline=0.0001, bline=0.0001, annotate_p=0.0001, chrcolor1="#3FA7D6", chrcolor2="#FAC05E", toptitle="Top Plot", bottomtitle="Bottom Plot")
```
![Imgur](https://i.imgur.com/1yrlwsk.png)

### Create a mirrored Manhattan plot using PheWAS data

```
library(hudson)
#Add some phenotype information
t$PHE <- rep(paste("Pheno", seq(1:5), sep="") , length.out=5000, each=1)
b$PHE <- rep(paste("Pheno", seq(1:5), sep="") , length.out=5000, each=1)
t <- t[, c(6,1:5)]
b <- b[, c(6,1:5)]

#Create a plot with highlight and annotation by p-value threshold
phemirror(top=t, bottom=b, highlight_p = 0.0001, highlighter="green", annotate_p=0.0001)
```
![Imgur](https://i.imgur.com/XM9sJ4z.jpg)

### Create a mirrored Manhattan plot using EWAS data
```
library(hudson)
#Generate some data
t <- data.frame(Variable=paste("Var", seq(1:5000), sep=""), 
                pvalue=runif(n=5000), 
                Group=rep(paste("G", seq(1:6), sep=""), length.out=5000, each=1),
                Shape=rep(paste("S", seq(1:5), sep="") , length.out=5000, each=1))
b <- data.frame(Variable=paste("Var", seq(1:5000), sep=""),
                pvalue=runif(n=5000), 
                Group=rep(paste("G", seq(1:6), sep=""), length.out=5000, each=1),
                Shape=rep(paste("S", seq(1:5), sep="") , length.out=5000, each=1))
head(t)
  Variable    pvalue Group Shape
1     Var1 0.1430880    G1    S1
2     Var2 0.6055253    G2    S2
3     Var3 0.8830645    G3    S3
4     Var4 0.1560467    G4    S4
5     Var5 0.7784929    G5    S5
6     Var6 0.9709326    G6    S1
head(b)
  Variable    pvalue Group Shape
1     Var1 0.5381795    G1    S1
2     Var2 0.9858681    G2    S2
3     Var3 0.7146113    G3    S3
4     Var4 0.3732309    G4    S4
5     Var5 0.3375573    G5    S5
6     Var6 0.4545325    G6    S1

#Generate a plot and highlight by variable name
emirror(top=t, bottom=b, highlight_var = c("Var1266", "Var3682"), color1="#587291", color2="#2F97C1")

```

![Imgur](https://i.imgur.com/NUkQ9jb.png)

Note that although you can rotate the axis labels by changing the ```rotatelabel``` and ```labelangle``` parameters, you'll probaby want to keep your "Group" names pretty short if some of your categories don't have a lot of variables in them.

