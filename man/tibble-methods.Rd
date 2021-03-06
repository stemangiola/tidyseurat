% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tibble_methods.R
\name{as_tibble}
\alias{as_tibble}
\alias{glimpse}
\title{Coerce lists, matrices, and more to data frames}
\arguments{
\item{rownames}{How to treat existing row names of a data frame or matrix:
 * `NULL`: remove row names. This is the default.
 * `NA`: keep row names.
 * A string: the name of a new column. Existing rownames are transferred
   into this column and the `row.names` attribute is deleted.
Read more in [rownames].}

\item{.name_repair}{see tidyr

  For compatibility only, do not use for new code.}

\item{x}{An object to glimpse at.}

\item{width}{Width of output: defaults to the setting of the option
`tibble.width` (if finite) or the width of the console.}

\item{...}{Unused, for extensibility.}
}
\value{
A tibble

x original x is (invisibly) returned, allowing `glimpse()` to be
  used within a data pipe line.
}
\description{
`r lifecycle::badge("maturing")`

`as_tibble()` turns an existing object, such as a data frame or
matrix, into a so-called tibble, a data frame with class [`tbl_df`]. This is
in contrast with [tibble()], which builds a tibble from individual columns.
`as_tibble()` is to [`tibble()`] as [base::as.data.frame()] is to
[base::data.frame()].

`as_tibble()` is an S3 generic, with methods for:
* [`data.frame`][base::data.frame()]: Thin wrapper around the `list` method
  that implements tibble's treatment of [rownames].
* [`matrix`][methods::matrix-class], [`poly`][stats::poly()],
  [`ts`][stats::ts()], [`table`][base::table()]
* Default: Other inputs are first coerced with [base::as.data.frame()].

`r lifecycle::badge("maturing")`

`glimpse()` is like a transposed version of `print()`:
columns run down the page, and data runs across.
This makes it possible to see every column in a data frame.
It's a little like [str()] applied to a data frame
but it tries to show you as much data as possible.
(And it always shows the underlying data, even when applied
to a remote data source.)

This generic will be moved to \pkg{pillar}, and reexported from there
as soon as it becomes available.
}
\section{Row names}{

The default behavior is to silently remove row names.

New code should explicitly convert row names to a new column using the
`rownames` argument.

For existing code that relies on the retention of row names, call
`pkgconfig::set_config("tibble::rownames" = NA)` in your script or in your
package's [.onLoad()]  function.
}

\section{Life cycle}{

Using `as_tibble()` for vectors is superseded as of version 3.0.0,
prefer the more expressive maturing `as_tibble_row()` and
`as_tibble_col()` variants for new code.
}

\section{S3 methods}{

`glimpse` is an S3 generic with a customised method for `tbl`s and
`data.frames`, and a default method that calls [str()].
}

\examples{
data("pbmc_small")
pbmc_small \%>\%  as_tibble()
data("pbmc_small")
pbmc_small \%>\% tidy \%>\% glimpse()


}
\seealso{
[tibble()] constructs a tibble from individual columns. [enframe()]
  converts a named vector to a tibble with a column of names and column of
  values. Name repair is implemented using [vctrs::vec_as_names()].
}
