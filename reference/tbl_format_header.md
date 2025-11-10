# Format the header of a tibble

**\[experimental\]**

For easier customization, the formatting of a tibble is split into three
components: header, body, and footer. The `tbl_format_header()` method
is responsible for formatting the header of a tibble.

Override this method if you need to change the appearance of the entire
header. If you only need to change or extend the components shown in the
header, override or extend
[`tbl_sum()`](https://pillar.r-lib.org/reference/tbl_sum.html) for your
class which is called by the default method.

## Usage

``` r
# S3 method for class 'tidySeurat'
tbl_format_header(x, setup, ...)
```

## Arguments

- x:

  A tibble-like object.

- setup:

  A setup object returned from
  [`tbl_format_setup()`](https://pillar.r-lib.org/reference/tbl_format_setup.html).

- ...:

  These dots are for future extensions and must be empty.

## Value

A character vector.

## Examples

``` r
# TODO
```
