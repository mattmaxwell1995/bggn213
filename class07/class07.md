Class 7 Functions and packages
================
Matt Maxwell
January 30, 2019

Functions revisited
-------------------

Load(i.e. **source**) our rescale() function from last class.

``` r
source("http://tinyurl.com/rescale-R")
```

we want to make this function more robust to these types of errors

``` r
#rescale2( c(1:5, "string"))
```

``` r
is.numeric( c(1:5, "string"))
```

    ## [1] FALSE

``` r
!is.numeric( c(1:5, "string"))
```

    ## [1] TRUE

``` r
!is.numeric(1:5)
```

    ## [1] FALSE

``` r
x <- c( 1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)
```

``` r
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE FALSE

``` r
is.na(x) & is.na(y)
```

    ## [1] FALSE FALSE  TRUE FALSE FALSE

``` r
sum( is.na(x) & is.na(y))
```

    ## [1] 1

``` r
both_na <- function(x,y) {
  ##check for NA elements in both input vectors
  sum( is.na(x) & is.na(y))
}
```

``` r
both_na(x, y)
```

    ## [1] 1

``` r
x <- c(NA, NA, NA)
y1 <- c( 1, NA, NA)
y2 <- c( 1, NA, NA, NA)
```

``` r
both_na(x, y1)
```

    ## [1] 2

``` r
both_na(x, y2)
```

    ## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
    ## shorter object length

    ## [1] 3

CRAN & Bioconductor
-------------------
