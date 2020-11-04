#' Plots the fraction of cells in each cluster compared between sample or group.
#' @param GEX.object Output of the automate_GEX function
#' @param per_sample Logical indicating if the grouping is done at the level of the samples or the groups
#' @export
#' @examples
#' \dontrun{
#' GEX_fraction_cells_per_cluster(automate_GEX.output, per_sample))
#' }

GEX_fraction_cells_per_cluster <- function(GEX.object, per_sample) {
  require(ggplot2)
  require(reshape2)
  require(stringr)
  if (per_sample) {
    temp_count <- list()
    for (i in unique(GEX.object[[1]]$sample_id)) {
      temp_count[[i]] <- data.frame(cbind(str_glue("sample_", i),
                                          table(GEX.object[[1]]$seurat_clusters[which(GEX.object[[1]]$sample_id==i)])))
    }
    temp_count <- do.call(rbind, temp_count)
    temp_count$cluster <- as.numeric(rownames(temp_count))
    colnames(temp_count) <-  c("sample_id", "count", "cluster")
    temp_count$count <- as.numeric(temp_count$count)
    print(ggplot(temp_count, aes(x = cluster, y = count, fill = sample_id)) + geom_bar(stat="identity", position=position_dodge()))
  }
  if (!per_sample) {
    temp_count <- list()
    for (i in unique(GEX.object[[1]]$group_id)) {
      temp_count[[i]] <- data.frame(cbind(str_glue("group_", i),
                                          table(GEX.object[[1]]$seurat_clusters[which(GEX.object[[1]]$group_id==i)])))
    }
    temp_count <- do.call(rbind, temp_count)
    temp_count$cluster <- as.numeric(rownames(temp_count))
    colnames(temp_count) <-  c("group_id", "count", "cluster")
    temp_count$count <- as.numeric(temp_count$count)
    print(ggplot(temp_count, aes(x = cluster, y = count, fill = group_id)) + geom_bar(stat="identity", position=position_dodge()))
  }
}
