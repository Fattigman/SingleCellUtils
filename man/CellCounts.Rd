% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{CellCounts}
\alias{CellCounts}
\title{Cell Counts Analysis}
\usage{
CellCounts(seurat_object, meta_data_var = "orig.ident")
}
\arguments{
\item{seurat_object}{A Seurat object containing single-cell data.}

\item{meta_data_var}{A string specifying the metadata variable to analyze.}
}
\value{
A bar plot showing the cell count per metadata variable.
}
\description{
This function calculates the number of cells per metadata variable and creates a bar plot.
}
