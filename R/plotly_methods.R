#' @name plotly
#' @rdname plotly
#' @inherit plotly::plot_ly
#' @return `plotly`
#' 
#' @examples
#' data(pbmc_small)
#' plot_ly(pbmc_small)
#' 
#' @importFrom plotly plot_ly
#' @export
plot_ly <- function(data=data.frame(), ..., type=NULL, name=NULL,
    color=NULL, colors=NULL, alpha=NULL,
    stroke=NULL, strokes=NULL, alpha_stroke=1,
    size=NULL, sizes=c(10, 100),
    span=NULL, spans=c(1, 20),
    symbol=NULL, symbols=NULL,
    linetype=NULL, linetypes=NULL,
    split=NULL, frame=NULL,
    width=NULL, height=NULL, source="A") {
    UseMethod("plot_ly")
}

#' @rdname plotly
#' @export
plot_ly.tbl_df <- function(data=data.frame(), ..., type=NULL, name=NULL,
    color=NULL, colors=NULL, alpha=NULL,
    stroke=NULL, strokes=NULL, alpha_stroke=1,
    size=NULL, sizes=c(10, 100),
    span=NULL, spans=c(1, 20),
    symbol=NULL, symbols=NULL,
    linetype=NULL, linetypes=NULL,
    split=NULL, frame=NULL,
    width=NULL, height=NULL, source="A") {
    data %>%
        # This is a trick to not loop the call
        drop_class("tbl_df") %>%
        plotly::plot_ly(...,
            type=type, name=name,
            color=color, colors=colors, alpha=alpha,
            stroke=stroke, strokes=strokes, alpha_stroke=alpha_stroke,
            size=size, sizes=sizes,
            span=span, spans=spans,
            symbol=symbol, symbols=symbols,
            linetype=linetype, linetypes=linetypes,
            split=split, frame=frame,
            width=width, height=height, source=source)
}

#' @rdname plotly
#' @export
plot_ly.Seurat <- function(data=data.frame(), ..., type=NULL, name=NULL,
    color=NULL, colors=NULL, alpha=NULL, 
    stroke=NULL, strokes=NULL, alpha_stroke=1,
    size=NULL, sizes=c(10, 100), 
    span=NULL, spans=c(1, 20),
    symbol=NULL, symbols=NULL,
    linetype=NULL, linetypes=NULL,
    split=NULL, frame=NULL, 
    width=NULL, height=NULL, source="A") {
  
    data %>%
        # This is a trick to not loop the call
        as_tibble() %>%
        plot_ly(..., type=type, name=name,
	        color=color, colors=colors, alpha=alpha, 
	        stroke=stroke, strokes=strokes, alpha_stroke=alpha_stroke,
	        size=size, sizes=sizes, 
	        span=span, spans=spans,
	        symbol=symbol, symbols=symbols, 
	        linetype=linetype, linetypes=linetypes,
	        split=split, frame=frame, 
	        width=width, height=height, source=source)    
}
