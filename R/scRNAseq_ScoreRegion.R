#' Calculate Region preference for each cluster
#'
#' @param SeuratObj Seurat object
#' @param BulkRNAseq.expr A data frame: Bulk RNA-seq gene expression values of each region, should have same gene names(row names) as in scRNA-seq
#' @param UMI.gradient A numeric vector: UMI cut-off for selecting genes with enough summed expression
#' @param Genes.gradient A numeric vector: selected genes numbers for all regions.
#' @param minimal.RNAseq.value A numeric: remove genes with low expression value in bulk RNA-seq
#' @param scRNAseq.expression.cut  A numeric: cut off of the scaled value to define whether a gene is expressed or not.
#'
#' @return A list: cluster * Region preference for each UMI And Top N genes combination
#' @export
#'
#' @examples
#' data(scRNA)
#' data(bulkRNA)
#' scorelist <- scRNAseq_Score_Region(scRNA, bulkRNA)
scRNAseq_Score_Region <- function(SeuratObj,
                                  BulkRNAseq.expr,
                                  # filter genes by summed UMI count from all clusters, remove genes with low expression in scRNA-seq.
                                  UMI.gradient = c(10, 20, 30, 40, 50, 100, 200, 500, 1000, 1500, 2000),
                                  # select top genes used
                                  Genes.gradient = c(10, 20, 30, 40, 50, 100, 200, 300, 400, 500),
                                  # filter genes by expression value in Bulk RNA-seq
                                  minimal.RNAseq.value = 0.5,
                                  # cut off of the scaled value to define if a gene is expressed or not.
                                  scRNAseq.expression.cut = 0.5


){
  n.region <- ncol(BulkRNAseq.expr) # how many regions
  sorted.list <- list() # save the ranked fold change results for each region
  # for each region, choose genes with RPKM value >=3.5, and calculate the fold change compared with other regions.
  for (i in 1:n.region) {
    expr_tmp <- BulkRNAseq.expr[BulkRNAseq.expr[,i] >= 3.5,]  # remove genes bellow the minimal expression/RPKM value!
    expr_tmp <- dplyr::mutate(expr_tmp,fc = expr_tmp[,i]/apply(expr_tmp[,-i],1,function(x)mean(x))) # calculate fold change
    # results <- dplyr::arrange(expr_tmp,dplyr::desc(fc)) # rank the genes by fold change
    results <- expr_tmp[order(expr_tmp$fc, decreasing = TRUE),]
    sorted.list[[colnames(BulkRNAseq.expr)[i]]] <- results
  }
  # annotated with summed UMI counts from scRNA-seq data
  sorted.list <- lapply(sorted.list, function(x){
    x$UMIs <- apply(SeuratObj@assays$RNA@counts[rownames(x),],1,sum)
    return(x)
  })
  # list: cell * gene expression matrix of each cluster, using scaled data(why?)
  cluster.expr.list <- split(data.frame(t(SeuratObj@assays$RNA@scale.data)), as.character(SeuratObj@active.ident))
  # 计算每个UMI和Top x Genes数目组合下，每个cluster的Region打分
  main_fun <- function(cluster.expr.list = cluster.expr.list, # cell * gene
                       sorted.list = sorted.list,
                       UMI.gradient = UMI.gradient,
                       Genes.gradient = Genes.gradient){
    # for each combination of UMI and Top n genes, for each cluster, calculate the region preference.
    # return a score matrix: Cluster * Region
    region_cluster_mat <- function(cluster.expr.list = cluster.expr.list, gene.list = gene.list){
      # build results matrix
      res.sub <- matrix(0,
                        nrow = length(cluster.expr.list),
                        ncol = length(gene.list),
                        dimnames = list(names(cluster.expr.list),names(gene.list)))
      for (x in names(cluster.expr.list)) {
        for (y in names(gene.list)) {
          expr.sub <- cluster.expr.list[[x]]
          # remove genes not exist in scRNA-seq
          genes.region <- gene.list[[y]][gene.list[[y]] %in% colnames(expr.sub)]
          expr.sub <- expr.sub[,genes.region]
          # for each region, calculate the enrichment score of each scRNA-seq cluster.
          expr.sub <- expr.sub > scRNAseq.expression.cut # transfer to a expression binary matrix, 0:not expressed, 1:expressed.
          cell.percent <- apply(expr.sub, 2, mean) # mean expression percent of each gene
          res.sub[x,y] <- mean(cell.percent) # mean expression percent of all genes which is normalized by genes numbers, for different gene sets may have different gene numbers
        }
      }
      return(res.sub)
    }
    # save the final results for each combination
    res <- list()
    for (i in as.character(UMI.gradient)) {
      for (j in as.character(Genes.gradient)) {
        # for each combination, choose genes with enough UMIs, and only keep the top_n genes
        gene.list <- lapply(sorted.list, function(x){
          x <- x[x$UMIs > i,]
          return(rownames(x)[1:j])
        })
        # calculate the region preference of each cluster for each combination
        res[[i]][[j]] <- region_cluster_mat(cluster.expr.list = cluster.expr.list, gene.list = gene.list)
      }
    }
    return(res)
  }
  # run
  main_fun(cluster.expr.list = cluster.expr.list, # cell * gene
           sorted.list = sorted.list,
           UMI.gradient = UMI.gradient,
           Genes.gradient = Genes.gradient)
}

#' evaluate the UMI And top n genes combination by Gini index
#' hopefully find a combination with a relatively high Gini index
#'
#' @param ScoreList A list: Output from scRNAseq_Score_Region function
#'
#' @return a heat-map plot of each combination(x: genes used; y: UMI cutoff)
#' @export
#'
#' @examples
#' data(scRNA)
#' data(bulkRNA)
#' scorelist <- scRNAseq_Score_Region(scRNA, bulkRNA)
#' scRNAseq_Score_Region_evaluate(scorelist)
scRNAseq_Score_Region_evaluate <- function(ScoreList){
  gini.list <- lapply(ScoreList, function(x)(unlist(lapply(x, function(x){sum(apply(x, 2, ineq::ineq))}))))
  res.gini <- Reduce(rbind, gini.list)
  rownames(res.gini) <- names(gini.list)
  pheatmap::pheatmap(res.gini, cluster_rows = FALSE, cluster_cols = FALSE)
}


#' Heatmap of the region preference of each cluster
#' default use the UMI And Top n Genes with the bigest Gini index
#'
#' @param ScoreList A list: Output from scRNAseq_Score_Region function
#' @param UMI UMI cutoff defined in scRNAseq_Score_Region function - UMI.gradient
#' @param TopGene Top n genes defined in scRNAseq_Score_Region function - Genes.gradient
#'
#' @return A heat map of cluster preference in each region
#' @export
#'
#' @examples
#' data(scRNA)
#' data(bulkRNA)
#' scorelist <- scRNAseq_Score_Region(scRNA, bulkRNA)
#' scRNAseq_Score_Region_plot(scorelist)
#' scRNAseq_Score_Region_plot(scorelist, 100, 100)
scRNAseq_Score_Region_plot <- function(ScoreList,UMI = NULL, TopGene = NULL){
  if (is.null(UMI) | is.null(TopGene)) {
    gini.list <- lapply(ScoreList, function(x)(unlist(lapply(x, function(x){sum(apply(x, 2, ineq::ineq))}))))
    res.gini <- Reduce(rbind, gini.list)
    rownames(res.gini) <- names(gini.list)
    umi.gene.max <- which(res.gini == res.gini[which.max(res.gini)],arr.ind=T)
    message("Using UMI: ", umi.gene.max[1],"; Top Genes: ",umi.gene.max[2])
    pheatmap::pheatmap(ScoreList[[umi.gene.max[1]]][[umi.gene.max[2]]])
  }else{
    pheatmap::pheatmap(ScoreList[[as.character(UMI)]][[as.character(TopGene)]])
  }
}
