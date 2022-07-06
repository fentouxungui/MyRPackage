#' Trans row names of a data frame according to a meta file
#'
#' @param BulkRNAseq.expr A data frame: Bulk RNA-seq gene expression values of each region, should have same gene names(row names) as in scRNA-seq
#' @param Meta A data frame with new name and old name
#' @param from old name
#' @param to new name
#' @param output.lost.genes output the genes exist in RNA-seq data while not exist in scRNA-seq data, Default FALSE.
#'
#' @return A data frame with row names changed (rows may not as long as before)
#' @export
#'
#' @examples
#' data(FlyGeneMeta)
#' data(RNAseq)
#' head(RNAseq$EC)
#' head(scRNAseq_Score_Region_Check(RNAseq$EC, FlyGeneMeta))
scRNAseq_Score_Region_Check <- function(BulkRNAseq.expr,
                                        Meta,
                                        from = "gene_id",
                                        to = "gene_name",
                                        output.lost.genes = FALSE){
  # how many row names of BulkRNAseq.expr not exist in Meta.
  not.exists <- rownames(BulkRNAseq.expr)[!rownames(BulkRNAseq.expr) %in% Meta[,from]]
  message(paste0(length(not.exists), " features from data frame not exist in meta file!"))
  if (output.lost.genes) {
    cat("Genes exist in RNA-seq data while not exist in scRNAseq data:")
    cat(not.exists)
  }
  BulkRNAseq.expr <- BulkRNAseq.expr[rownames(BulkRNAseq.expr) %in%  Meta[,from],]
  mapping <- Meta[,to]
  names(mapping) <- Meta[,from]
  rownames(BulkRNAseq.expr) <- mapping[rownames(BulkRNAseq.expr)]
  return(BulkRNAseq.expr)
}

#' Check Bulk RNA-seq Region data file
#'
#' @param Adf A data frame: Bulk RNA-seq gene expression values of each region, should have same gene names(row names) as in scRNA-seq
#'
#' @return A data frame with all columns are numeric
#'
#' @examples
#' # Not Run
#' # data(RNAseq)
#' # bulkRNA <- check_RNAseq_df(RNAseq$EC)
check_RNAseq_df <- function(Adf){
  for (i in colnames(Adf)) {
    Adf[,i] <- as.numeric(Adf[,i])
  }
  return(Adf)
}
#' Calculate Region preference for each cluster
#'
#' @param SeuratObj Seurat object
#' @param BulkRNAseq.expr A data frame: Bulk RNA-seq gene expression values of each region, should have same gene names(row names) as in scRNA-seq
#' @param UMI.gradient A numeric vector: UMI cut-off for selecting genes with enough summed expression
#' @param Genes.gradient A numeric vector: selected genes numbers for all regions.
#' @param minimal.RNAseq.value A numeric: remove genes with low expression value in bulk RNA-seq
#' @param scRNAseq.expression.cut  A numeric: cut off of the scaled value to define whether a gene is expressed or not.
#' @param Output.lost.genes output the genes exist in RNA-seq data while not exist in scRNA-seq data, Default FALSE.
#' @param Customized.GeneList A list: customized regional markers
#'
#' @return A list: 1. selected genes for each region; 2. cluster * Region preference for each UMI And Top N genes combination
#' @export
#'
#' @examples
#' data(scRNA)
#' data(RNAseq)
#' bulkRNAseq <- scRNAseq_Score_Region_Check(RNAseq$EC, FlyGeneMeta)
#' score.list <- scRNAseq_Score_Region(scRNA, bulkRNAseq)
#' gene.list <- scRNAseq_Score_Region_ExtractFeatures(score.list)
#' GeneList <- list()
#' for(i in colnames(gene.list)){ GeneList[[i]] <- gene.list[,i]}
#' score.list.customized <- scRNAseq_Score_Region(SeuratObj = scRNA, Customized.GeneList = GeneList)
scRNAseq_Score_Region <- function(SeuratObj,
                                  BulkRNAseq.expr,
                                  # filter genes by summed UMI count from all clusters, remove genes with low expression in scRNA-seq.
                                  UMI.gradient = c(10, 20, 30, 40, 50, 100, 200, 500, 1000, 1500, 2000),
                                  # select top genes used
                                  Genes.gradient = c(10, 20, 30, 40, 50, 100, 200, 300, 400, 500),
                                  # filter genes by expression value in Bulk RNA-seq
                                  minimal.RNAseq.value = 3.5,
                                  # cut off of the scaled value to define if a gene is expressed or not.
                                  scRNAseq.expression.cut = 0.5,
                                  Output.lost.genes = FALSE,
                                  Customized.GeneList = NULL


){
  # list: cell * gene expression matrix of each cluster, using scaled data(why?)
  cluster.expr.list <- split(data.frame(t(SeuratObj@assays$RNA@scale.data)), as.character(SeuratObj@active.ident))
  # for each region in gene list and for each cluster, calculate the region preference.
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
    return(list(Genes = gene.list, Score = res.sub))
  }

  # 使用参数过滤的基因list
  if (is.null(Customized.GeneList)) {
    BulkRNAseq.expr <- check_RNAseq_df(BulkRNAseq.expr)
    # RNAseq里的基因是否都存在于scRNAseq中
    not.exists <- rownames(BulkRNAseq.expr)[!rownames(BulkRNAseq.expr) %in% rownames(SeuratObj)]
    message(paste0(length(not.exists), " features from RNA-seq not exist in scRNAseq!"))
    if (Output.lost.genes) {
      cat("Genes exist in RNA-seq data while not exist in scRNAseq data:")
      cat(not.exists)
    }
    BulkRNAseq.expr <- BulkRNAseq.expr[rownames(BulkRNAseq.expr) %in% rownames(SeuratObj),]
    n.region <- ncol(BulkRNAseq.expr) # how many regions
    sorted.list <- list() # save the ranked fold change results for each region
    # for each region, choose genes with RPKM value >=minimal.RNAseq.value, and calculate the fold change compared with other regions.
    for (i in 1:n.region) {
      expr_tmp <- BulkRNAseq.expr[BulkRNAseq.expr[,i] >= minimal.RNAseq.value,]  # remove genes bellow the minimal expression/RPKM value!
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
    # 计算每个UMI和Top x Genes数目组合下，每个cluster的Region打分
    main_fun <- function(cluster.expr.list = cluster.expr.list, # cell * gene
                         sorted.list = sorted.list,
                         UMI.gradient = UMI.gradient,
                         Genes.gradient = Genes.gradient){
      # for each combination of UMI and Top n genes, for each cluster, calculate the region preference.
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
    return(main_fun(cluster.expr.list = cluster.expr.list, # cell * gene
             sorted.list = sorted.list,
             UMI.gradient = UMI.gradient,
             Genes.gradient = Genes.gradient))
  }else{
    # Gene list里的基因是否都存在于scRNAseq中
    not.exists <- unname(unlist(Customized.GeneList))[!unname(unlist(Customized.GeneList)) %in% rownames(SeuratObj)]
    message(paste0(length(not.exists), " features from RNA-seq not exist in scRNAseq!"))
    if (Output.lost.genes) {
      cat("Genes exist in RNA-seq data while not exist in scRNAseq data:")
      cat(not.exists)
    }
    Customized.GeneList <- lapply(Customized.GeneList, function(x){x[x %in% rownames(SeuratObj)]})
    res <- list()
    res[["customized"]][["customized"]]<- region_cluster_mat(cluster.expr.list, Customized.GeneList)
    return(res)
  }
}

#' evaluate the UMI And top n genes combination by Gini index
#'
#' To find a combination with a relatively high Gini index
#'
#' @param ScoreList A list: Output from scRNAseq_Score_Region function
#' @param ... Other parameters passed to pheatmap() function
#'
#' @return a Gini index heatmap of all combinations (x: genes used; y: UMI cutoff)
#' @export
#'
#' @examples
#' data(scRNA)
#' data(RNAseq)
#' bulkRNAseq <- scRNAseq_Score_Region_Check(RNAseq$EC, FlyGeneMeta)
#' score.list <- scRNAseq_Score_Region(scRNA, bulkRNAseq)
#' scRNAseq_Score_Region_evaluate(score.list)
scRNAseq_Score_Region_evaluate <- function(ScoreList, ...){
  gini.list <- lapply(ScoreList, function(x)(unlist(lapply(x, function(x){sum(apply(x$Score, 2, ineq::ineq))}))))
  res.gini <- Reduce(rbind, gini.list)
  rownames(res.gini) <- names(gini.list)
  return(pheatmap::pheatmap(res.gini, ...))
}


#' Heatmap of the region preference of each cluster
#' default use the UMI And Top n Genes with the bigest Gini index
#'
#' @param ScoreList A list: Output from scRNAseq_Score_Region function
#' @param UMI UMI cutoff defined in scRNAseq_Score_Region function - UMI.gradient
#' @param TopGene Top n genes defined in scRNAseq_Score_Region function - Genes.gradient
#' @param ... Other parameters passed to pheatmap() function
#'
#' @return A heat map of cluster preference in each region
#' @export
#'
#' @examples
#' data(scRNA)
#' data(RNAseq)
#' bulkRNAseq <- scRNAseq_Score_Region_Check(RNAseq$EC, FlyGeneMeta)
#' score.list <- scRNAseq_Score_Region(scRNA, bulkRNAseq)
#' scRNAseq_Score_Region_plot(score.list)
#' scRNAseq_Score_Region_plot(score.list, 100, 100)
scRNAseq_Score_Region_plot <- function(ScoreList,UMI = NULL, TopGene = NULL, ...){
  if (is.null(UMI) | is.null(TopGene)) {
    # 按Region计算各群表达值的Gini系数，然后对5个Region求和。
    gini.list <- lapply(ScoreList, function(x)(unlist(lapply(x, function(x){sum(apply(x$Score, 2, ineq::ineq))}))))
    res.gini <- Reduce(rbind, gini.list)
    rownames(res.gini) <- names(gini.list)
    # 有时最大值可能会有多个，用第一个
    umi.gene.max <- which(res.gini == res.gini[which.max(res.gini)], arr.ind=T)[1,]
    message("Using UMI Cutoff: ", rownames(res.gini)[umi.gene.max[1]],"; Genes Used: ", colnames(res.gini)[umi.gene.max[2]])
    return(pheatmap::pheatmap(ScoreList[[umi.gene.max[1]]][[umi.gene.max[2]]]$Score, ...))
  }else{
    return(pheatmap::pheatmap(ScoreList[[as.character(UMI)]][[as.character(TopGene)]]$Score, ...))
  }
}

#' Extract features
#'
#' @param ScoreList A list: Output from scRNAseq_Score_Region function
#' @param UMI UMI cutoff defined in scRNAseq_Score_Region function - UMI.gradient
#' @param TopGene Top n genes defined in scRNAseq_Score_Region function - Genes.gradient
#'
#' @return 1. A venn diagram of the genes from different regions; 2.A data frame of the extracted genes for each region.
#' @export
#'
#' @examples
#' data(scRNA)
#' data(RNAseq)
#' bulkRNAseq <- scRNAseq_Score_Region_Check(RNAseq$EC, FlyGeneMeta)
#' score.list <- scRNAseq_Score_Region(scRNA, bulkRNAseq)
#' ## futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
#' scRNAseq_Score_Region_ExtractFeatures(score.list)
#' scRNAseq_Score_Region_ExtractFeatures(score.list, 100, 100)
scRNAseq_Score_Region_ExtractFeatures <- function(ScoreList,UMI = NULL, TopGene = NULL){
  if (is.null(UMI) | is.null(TopGene)) {
    # 按Region计算各群表达值的Gini系数，然后对5个Region求和。
    gini.list <- lapply(ScoreList, function(x)(unlist(lapply(x, function(x){sum(apply(x$Score, 2, ineq::ineq))}))))
    res.gini <- Reduce(rbind, gini.list)
    rownames(res.gini) <- names(gini.list)
    # 有时最大值可能会有多个，用第一个
    umi.gene.max <- which(res.gini == res.gini[which.max(res.gini)], arr.ind=T)[1,]
    UMI <- rownames(res.gini)[umi.gene.max[1]]
    TopGene <- colnames(res.gini)[umi.gene.max[2]]
  }
  message("Using UMI Cutoff: ", UMI,"; Genes Used: ", TopGene)
  Gene.list <- ScoreList[[as.character(UMI)]][[as.character(TopGene)]]$Genes
  # futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  print(grid::grid.draw( VennDiagram::venn.diagram(Gene.list, filename = NULL)))
  res <- Reduce(cbind, Gene.list)
  colnames(res) <- names(Gene.list)
  return(res)
}

#' plot Euclidean Distance between each combination
#'
#' @param ScoreList A list: Output from scRNAseq_Score_Region function
#'
#' @return a heat-map plot of the distance of each combination(genes used and UMI cutoff)
#' @export
#'
#' @examples
#' data(scRNA)
#' data(RNAseq)
#' bulkRNAseq <- scRNAseq_Score_Region_Check(RNAseq$EC, FlyGeneMeta)
#' score.list <- scRNAseq_Score_Region(scRNA, bulkRNAseq)
#' scRNAseq_Score_Region_evaluate2(score.list)
scRNAseq_Score_Region_evaluate2 <- function(ScoreList){
  res.list <- list()
  for (i in names(ScoreList)) {
    for (j in names(ScoreList[[1]])) {
      res.list[[paste0("UMI-",i,"-Genes-",j)]] <- ScoreList[[i]][[j]]$Score
    }
  }
  cal_dis <- function(Alist){
    dis.mat <- matrix(0, nrow = length(Alist), ncol = length(Alist), dimnames = list(names(Alist), names(Alist)))
    for (i in names(Alist)) {
      for (j in names(Alist)) {
        # Euclidean Distance
        dis <- sqrt(sum((Alist[[i]] - Alist[[j]])^2))
        # KLD 距离
        # dis <- KLD(as.matrix(cnn.res.list[[i]]),as.matrix(cnn.res.list[[j]]))$mean.sum.KLD
        dis.mat[i,j] <- dis
      }
    }
    return(dis.mat)
  }
  dis.mat <- cal_dis(res.list)
  return(pheatmap::pheatmap(dis.mat))
}

#' Predict Region preference for each cluster By gene expression correlation
#'
#' @param SeuratObj Seurat object
#' @param BulkRNAseq.expr A data frame: Bulk RNA-seq gene expression values of each region, should have same gene names(row names) as in scRNA-seq
#' @param Method method used in cor function, one of "spearman", "pearson", "kendall".
#' @param Genes.Selection use all genes or only selected genes
#' @param Top.foldchange if use selected genes, the cut off of fold change for each Region compared with other regions
#' @param Top.UMI.Cutoff if use selected genes, the cut off of UMI for each Region
#' @param Top.numbers if use selected genes, the cut off of top gene numbers for each Region
#' @param Output.lost.genes output the genes exist in RNA-seq data while not exist in scRNA-seq data, Default FALSE.
#'
#' @return A Region * Cluster matrix
#' @export
#'
#' @examples
#' data(scRNA)
#' data(RNAseq)
#' bulkRNAseq <- scRNAseq_Score_Region_Check(RNAseq$EC, FlyGeneMeta)
#' score.matrix <- scRNAseq_Score_Region2(scRNA, bulkRNAseq, Method = "spearman")
#' pheatmap::pheatmap(score.matrix)
#' score.matrix <- scRNAseq_Score_Region2(scRNA, bulkRNAseq, Method = "spearman",
#'                                        Genes.Selection = "Top")
#' pheatmap::pheatmap(score.matrix)
scRNAseq_Score_Region2 <- function(SeuratObj,
                                  BulkRNAseq.expr,
                                  Method = c("spearman", "pearson", "kendall"),
                                  Genes.Selection = c("Simple","Top"),
                                  Top.foldchange = 3.5,
                                  Top.UMI.Cutoff = 1,
                                  Top.numbers = 300,
                                  Output.lost.genes = FALSE
){
  BulkRNAseq.expr <- check_RNAseq_df(BulkRNAseq.expr)
  # RNAseq里的基因是否都存在于scRNAseq中
  not.exists <- rownames(BulkRNAseq.expr)[!rownames(BulkRNAseq.expr) %in% rownames(SeuratObj)]
  message(paste0(length(not.exists), " features from RNA-seq not exist in scRNAseq!"))
  if (Output.lost.genes) {
    cat("Genes exist in RNA-seq data while not exist in scRNAseq data:")
    cat(not.exists)
  }
  BulkRNAseq.expr <- BulkRNAseq.expr[rownames(BulkRNAseq.expr) %in% rownames(SeuratObj),]
  # average expression value for each cluster.
  scRNAseq.averageExpr <- Seurat::AverageExpression(SeuratObj)[[1]]
  scRNAseq.averageExpr <- scRNAseq.averageExpr[rownames(scRNAseq.averageExpr) %in% rownames(BulkRNAseq.expr),]
  # log transformation for bulk RNA-seq RPKM value
  BulkRNAseq.expr <- BulkRNAseq.expr[rownames(scRNAseq.averageExpr),]
  BulkRNAseq.expr <- log(BulkRNAseq.expr + 1)
  # Simple
  if (Genes.Selection[1] == "Simple") {
    # 仅仅保留表达值和均不为0的基因
    index <- !(apply(BulkRNAseq.expr,1,sum) == 0 | apply(scRNAseq.averageExpr,1,sum) == 0)
    matrix <- stats::cor(BulkRNAseq.expr[index,], scRNAseq.averageExpr[index,],method=Method[1],use="everything")
    return(matrix)
  # Top
  }else if (Genes.Selection[1] == "Top") {
    selected.list <- list()
    for (i in colnames(BulkRNAseq.expr)) {
      expr_tmp <- BulkRNAseq.expr[BulkRNAseq.expr[,i] >= Top.foldchange,]  # the minimal RPKM value!
      expr_tmp <- dplyr::mutate(expr_tmp,fc = expr_tmp[,i]/apply(expr_tmp,1,function(x)mean(x))) # calculate fold change
      # results <- dplyr::arrange(expr_tmp,dplyr::desc(fc)) # rank the genes by fold change
      results <- expr_tmp[order(expr_tmp$fc, decreasing = TRUE),]
      selected.list[[i]] <- results # save variant
    }
    # annotated with summed UMI counts from scRNA-seq
    selected.list <- lapply(selected.list, function(x){
      x$UMIs <- apply(scRNAseq.averageExpr[rownames(x),],1,sum)
      return(x)
    })
    selected.list <- lapply(selected.list,function(x){x[x$UMIs > Top.UMI.Cutoff,]})
    selected.genes <- unique(unlist(lapply(selected.list,function(x)rownames(x)[1:Top.numbers])))
    selected.genes <- selected.genes[!is.na(selected.genes)]
    matrix <- stats::cor(BulkRNAseq.expr[selected.genes,], scRNAseq.averageExpr[selected.genes,],method=Method[1],use="everything")
    return(matrix)
  }else{
    stop("Please check the Genes.Selection parameter!")
  }
}

#' Compare scRNAseq_Score_Region and scRNAseq_Score_Region2 outputs by correlation method
#'
#' @param ScoreList A list: Output from scRNAseq_Score_Region function
#' @param ScoreMatrix A matrix: Output from scRNAseq_Score_Region2 function
#' @param Method method used in cor function, one of "spearman", "pearson", "kendall".
#'
#' @return A named Vector
#' @export
#'
#' @examples
#' data(scRNA)
#' data(RNAseq)
#' bulkRNAseq <- scRNAseq_Score_Region_Check(RNAseq$EC, FlyGeneMeta)
#' score.list <- scRNAseq_Score_Region(scRNA, bulkRNAseq)
#' score.matrix <- scRNAseq_Score_Region2(scRNA, bulkRNAseq, Method = "spearman",
#'                                        Genes.Selection = "Top")
#' scRNAseq_Score_Compare(score.list,score.matrix)
scRNAseq_Score_Compare <- function(ScoreList,
                                     ScoreMatrix,
                                     Method = c("spearman", "pearson", "kendall")){
  res.list <- list()
  for (i in names(ScoreList)) {
    for (j in names(ScoreList[[1]])) {
      res.list[[paste0("UMI-",i,"-Genes-",j)]] <- ScoreList[[i]][[j]]$Score
    }
  }
  cor.res <- c()
  for (i in names(res.list)) {
    tmp.df <- res.list[[i]][colnames(ScoreMatrix),]
    cor.res <- append(cor.res, stats::cor(as.vector(tmp.df),as.vector(t(ScoreMatrix)), method = Method[1]))
  }
  names(cor.res) <- names(res.list)
  return(sort(cor.res, decreasing = TRUE))
}

