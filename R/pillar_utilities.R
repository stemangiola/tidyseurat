NBSP <- "\U00A0"

pillar___format_comment = function (x, width)
{
  if (length(x) == 0L) {
    return(character())
  }
  map_chr(x, pillar___wrap, prefix = "# ", width = min(width, cli::console_width()))
}

#' @importFrom fansi strwrap_ctl
pillar___strwrap2 = function (x, width, indent)
{
  fansi::strwrap_ctl(x, width = max(width, 0), indent = indent,
                     exdent = indent + 2)
}


pillar___wrap = function (..., indent = 0, prefix = "", width)
{
  x <- paste0(..., collapse = "")
  wrapped <- pillar___strwrap2(x, width - get_extent(prefix), indent)
  wrapped <- paste0(prefix, wrapped)
  wrapped <- gsub(NBSP, " ", wrapped)
  paste0(wrapped, collapse = "\n")
}
