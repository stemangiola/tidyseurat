\name{NEWS}
\title{News for Package \pkg{tidyseurat}}

\section{Changes in version 0.5.1, Development}{
\itemize{
    \item Change default shape parameter in join_features() and join_transcripts() from "long" to "wide", resulting in a return type of Seurat by default
    \item Update documentation and tests accordingly
}}

\section{Changes in version 0.5.0, CRAN Release}{
\itemize{
    \item Rely of ttservice package for shared function with tidySingleCellExperiment to avoid clash
    \item Use .cell for cell column name to avoid errors when cell column is defined by the user
}}
