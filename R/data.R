#' RNA-seq RPKM value of four fly gut cell types from five gut regions
#'
#' This data set is a list with 4 data.frames, containing the RNA-seq RPKM value of four fly gut cell types(ISC/EB/EC/EE)
#' from five gut regions, which has a common column in all data.frames.
#'
#' @format list, with 4 data frames:
#' \enumerate{
#' \item \strong{ISC} A data frame with 15682 rows and 5 columns
#' \describe{
#'   \item{R1}{Gut Region 1: numeric}
#'   \item{R2}{Gut Region 2: numeric}
#'   \item{R3}{Gut Region 3: numeric}
#'   \item{R4}{Gut Region 4: numeric}
#'   \item{R5}{Gut Region 5: numeric}
#' }
#' \item \strong{EB} A data frame with 15682 rows and 5 columns
#' \describe{
#'   \item{R1}{Gut Region 1: numeric}
#'   \item{R2}{Gut Region 2: numeric}
#'   \item{R3}{Gut Region 3: numeric}
#'   \item{R4}{Gut Region 4: numeric}
#'   \item{R5}{Gut Region 5: numeric}
#' }
#' \item \strong{EC} A data frame with 15682 rows and 5 columns
#' \describe{
#'   \item{R1}{Gut Region 1: numeric}
#'   \item{R2}{Gut Region 2: numeric}
#'   \item{R3}{Gut Region 3: numeric}
#'   \item{R4}{Gut Region 4: numeric}
#'   \item{R5}{Gut Region 5: numeric}
#' }
#' \item \strong{EE} A data frame with 15682 rows and 5 columns
#' \describe{
#'   \item{R1}{Gut Region 1: numeric}
#'   \item{R2}{Gut Region 2: numeric}
#'   \item{R3}{Gut Region 3: numeric}
#'   \item{R4}{Gut Region 4: numeric}
#'   \item{R5}{Gut Region 5: numeric}
#' }
#' }
#' @source \url{http://flygutseq.buchonlab.com/resources}
"RNAseq"


#' Fly gut scRNA-seq by 10xGenomics
#'
#' Hung R J, Hu Y, Kirchner R, et al. A cell atlas of the adult Drosophila midgut. Proceedings of the National Academy
#' of Sciences, 2020, 117(3): 1514-1523.
#'
#' @format Seurat object, 2979 cells and 16960 features.
#' @source \url{https://www.pnas.org/doi/abs/10.1073/pnas.1916820117}
"scRNA"

#' Fly gene id and name correspondence
#'
#' Fly gene id - name transform
#'
#' @format data frame, 17558 rows and 2 variables:
#' \describe{
#'   \item{gene_id}{FBgn id}
#'   \item{gene_name}{Gene symbol}
#' }
#' @source \url{https://flybase.org/}
"FlyGeneMeta"
