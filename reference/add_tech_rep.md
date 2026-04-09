# Add technical replicates

Add counts across columns of a count matrix corresponding to technical
replicates.

## Usage

``` r
add_tech_rep(se, group)
```

## Arguments

- se:

  SummarizedExperiment object

- group:

  Character string vector, or factor, that indicates the grouping of
  technical replicates.

## Value

A SummarizedExperiment object, where technical replicates have been
added in the count matrix.

## Author

Robert Castelo

## Examples

``` r
library(SummarizedExperiment)

cnt <- matrix(1, nrow=3, ncol=9,
                dimnames=list(paste0("g", 1:3), paste0("s", 1:9)))
cdata <- data.frame(group=rep(LETTERS[1:3], each=3), row.names=colnames(cnt))
se <- SummarizedExperiment(assays=list(counts=cnt), colData=cdata)
dim(se)
#> [1] 3 9
assays(se)$counts
#>    s1 s2 s3 s4 s5 s6 s7 s8 s9
#> g1  1  1  1  1  1  1  1  1  1
#> g2  1  1  1  1  1  1  1  1  1
#> g3  1  1  1  1  1  1  1  1  1

se <- add_tech_rep(se, se$group)
dim(se)
#> [1] 3 3
assays(se)$counts
#>    A B C
#> g1 3 3 3
#> g2 3 3 3
#> g3 3 3 3
```
