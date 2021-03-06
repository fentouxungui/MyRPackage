% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scRNAseq_ScoreRegion.R
\name{scRNAseq_Score_Region_Check}
\alias{scRNAseq_Score_Region_Check}
\title{Trans row names of a data frame according to a meta file}
\usage{
scRNAseq_Score_Region_Check(
  BulkRNAseq.expr,
  Meta,
  from = "gene_id",
  to = "gene_name",
  output.lost.genes = FALSE
)
}
\arguments{
\item{BulkRNAseq.expr}{A data frame: Bulk RNA-seq gene expression values of each region, should have same gene names(row names) as in scRNA-seq}

\item{Meta}{A data frame with new name and old name}

\item{from}{old name}

\item{to}{new name}

\item{output.lost.genes}{output the genes exist in RNA-seq data while not exist in scRNA-seq data, Default FALSE.}
}
\value{
A data frame with row names changed (rows may not as long as before)
}
\description{
Trans row names of a data frame according to a meta file
}
\examples{
data(FlyGeneMeta)
data(RNAseq)
head(RNAseq$EC)
head(scRNAseq_Score_Region_Check(RNAseq$EC, FlyGeneMeta))
}
