% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scRNAseq_ScoreRegion.R
\name{scRNAseq_Score_Region_evaluate2}
\alias{scRNAseq_Score_Region_evaluate2}
\title{plot Euclidean Distance between each combination}
\usage{
scRNAseq_Score_Region_evaluate2(ScoreList)
}
\arguments{
\item{ScoreList}{A list: Output from scRNAseq_Score_Region function}
}
\value{
a heat-map plot of the distance of each combination(genes used and UMI cutoff)
}
\description{
plot Euclidean Distance between each combination
}
\examples{
data(scRNA)
data(RNAseq)
bulkRNAseq <- scRNAseq_Score_Region_Check(RNAseq$EC, FlyGeneMeta)
score.list <- scRNAseq_Score_Region(scRNA, bulkRNAseq)
scRNAseq_Score_Region_evaluate2(score.list)
}
