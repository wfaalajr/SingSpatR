#' @title sc_preprocess
#'
#' @description
#' A wrapper that reads in a merged Seurat object and does cell-level and gene-level filtering.
#' Also runs normalization (sctransform) and cell cycle scoring. Provides some plots before and
#' after filtering if qc.plots=T.
#'
#' @param merged_seurat merged Seurat object.
#' @param grouping.var metadata column from the merged Seurat object to be used as a grouping variable.
#' @param qc.plots logical, set to TRUE to get before and after filtering plots.
#'
#' @return a filtered subset of the provided merged Seurat object.
#'
#' @examples
#' #install.packages("https://seurat.nygenome.org/src/contrib/ifnb.SeuratData_3.0.0.tar.gz",
#'                   repos = NULL, type = "source")
#' # library(ifnb.SeuratData)
#' #ifnb <- SeuratData::LoadData("ifnb")
#' #ifnb.list <- Seurat::SplitObject(ifnb, split.by = "orig.ident")
#' #ctrl <- ifnb.list$IMMUNE_CTRL
#' #stim <- ifnb.list$IMMUNE_STIM
#' #merged <- merge(ctrl, stim)
#' #ctrl_stim_filtered <- sc_preprocess(merged_seurat=merged, grouping.var="orig.ident", qc.plots=T)
#'
#' @export
#' @importFrom dplyr "%>%"
#'
sc_preprocess <- function(merged_seurat = NULL,
                          grouping.var = "orig.ident",
                          qc.plots = T) {
  if(missing(merged_seurat)) {
    stop("merged_seurat is missing with no defaults.")
  }
  if(!inherits(merged_seurat, "Seurat")) {
    stop("merged_seurat is not a Seurat object.")
  }
  if(missing(grouping.var)) {
    warning("grouping.var is missing. Using 'orig.ident' instead.")
  }
  if(!(grouping.var %in% colnames(merged_seurat@meta.data))) {
    stop(paste("No", grouping.var, "in", merged_seurat@grouping.var, "metadata", sep = " "))
  }
  if (requireNamespace("magrittr", quietly = TRUE)) {
    magrittr::"%>%"
  }
  #Create some more QC metrics
  merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)
  merged_seurat$mitoRatio <- merged_seurat$percent.mt / 100
  #Pre-filter quality checks and plots
  a <- merged_seurat@meta.data %>%
    ggplot(aes(x=!!as.name(grouping.var), fill=!!as.name(grouping.var))) +
    geom_bar() +
    geom_text(aes(label = ..count..), stat = "count", vjust = 1.5) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 0)) +
    theme(plot.title = element_text(hjust=0.5)) +
    ggtitle("NCells (before filtering)") +
    theme(legend.position = "none") +

    merged_seurat@meta.data %>%
    ggplot(aes(color=!!as.name(grouping.var), x=nCount_RNA, fill=!!as.name(grouping.var))) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    ylab("log10(cell density)") +
    geom_vline(xintercept = 500) +
    guides(color = "none") +
    theme(legend.position = "top", legend.text = element_text(size = 5)) +
    ggtitle("Transcripts per cell (before filtering)") +
    theme(plot.title = element_text(hjust=0.5)) +

    merged_seurat@meta.data %>%
    ggplot(aes(color=!!as.name(grouping.var), x=nFeature_RNA, fill=!!as.name(grouping.var))) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    scale_x_log10() +
    geom_vline(xintercept = 300) +
    theme(legend.position = "top", legend.text = element_text(size = 5)) +
    guides(color = "none") +
    ggtitle("Distribution of genes per cell (before filtering)") +
    theme(plot.title = element_text(hjust=0.5)) +

    merged_seurat@meta.data %>%
    ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) +
    geom_point() +
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method=lm) +
    scale_x_log10() +
    scale_y_log10() +
    theme_classic() +
    geom_vline(xintercept = 500) +
    geom_hline(yintercept = 250) +
    facet_wrap(formula(paste("~", grouping.var))) +
    ggtitle("nGenes vs nUMIs (before filtering)") +
    theme(plot.title = element_text(hjust=0.5)) +

    merged_seurat@meta.data %>%
    ggplot(aes(color=!!as.name(grouping.var), x=mitoRatio, fill=!!as.name(grouping.var))) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    geom_vline(xintercept = 0.2) +
    ggtitle("Mitochondrial counts ratio (before filtering)") +
    theme(plot.title = element_text(hjust=0.5)) +
    theme(legend.position = "top", legend.text = element_text(size = 5)) +
    guides(color = "none") +

    merged_seurat@meta.data %>%
    ggplot(aes(x=log10GenesPerUMI, color = !!as.name(grouping.var), fill=!!as.name(grouping.var))) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    geom_vline(xintercept = 0.8) +
    ggtitle("Complexity (before filtering)") +
    theme(plot.title = element_text(hjust=0.5)) +
    theme(legend.position = "top", legend.text = element_text(size = 5)) +
    guides(color = "none")
  if(qc.plots == T) {
    print(a)
  }

  #Do cell-level filtering:
  merged_seurat_f <- subset(merged_seurat,
                            subset = nFeature_RNA >= 200 &
                              nFeature_RNA <= 4000 &
                              nCount_RNA >= 500 &
                              log10GenesPerUMI > 0.8 &
                              mitoRatio < 0.2)
  #Do gene-level filtering:
  counts <- Seurat::GetAssayData(object = merged_seurat_f, slot = "counts")
  nonzero <- counts > 0
  keep_genes <- Matrix::rowSums(nonzero) >= 10
  filtered_counts <- counts[keep_genes, ]
  merged_seurat_f <- Seurat::CreateSeuratObject(filtered_counts, meta.data = merged_seurat_f@meta.data)

  #Re-plot QC metrices
  b <- merged_seurat_f@meta.data %>%
    ggplot(aes(x=!!as.name(grouping.var), fill=!!as.name(grouping.var))) +
    geom_bar() +
    geom_text(aes(label = ..count..), stat = "count", vjust = 1.5) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 0)) +
    theme(plot.title = element_text(hjust=0.5)) +
    ggtitle("NCells (after filtering)") +
    theme(legend.position = "none") +

    merged_seurat_f@meta.data %>%
    ggplot(aes(color=!!as.name(grouping.var), x=nCount_RNA, fill=!!as.name(grouping.var))) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    ylab("log10(cell density)") +
    geom_vline(xintercept = 500) +
    guides(color = "none") +
    theme(legend.position = "top", legend.text = element_text(size = 5)) +
    ggtitle("Transcripts per cell (after filtering)") +
    theme(plot.title = element_text(hjust=0.5)) +

    merged_seurat_f@meta.data %>%
    ggplot(aes(color=!!as.name(grouping.var), x=nFeature_RNA, fill=!!as.name(grouping.var))) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    scale_x_log10() +
    geom_vline(xintercept = 300) +
    theme(legend.position = "top", legend.text = element_text(size = 5)) +
    guides(color = "none") +
    ggtitle("Distribution of genes per cell (after filtering)") +
    theme(plot.title = element_text(hjust=0.5)) +

    merged_seurat_f@meta.data %>%
    ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) +
    geom_point() +
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method=lm) +
    scale_x_log10() +
    scale_y_log10() +
    theme_classic() +
    geom_vline(xintercept = 500) +
    geom_hline(yintercept = 250) +
    facet_wrap(formula(paste("~", grouping.var))) +
    ggtitle("nGenes vs nUMIs (after filtering)") +
    theme(plot.title = element_text(hjust=0.5)) +

    merged_seurat_f@meta.data %>%
    ggplot(aes(color=!!as.name(grouping.var), x=mitoRatio, fill=!!as.name(grouping.var))) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    geom_vline(xintercept = 0.2) +
    ggtitle("Mitochondrial counts ratio (after filtering)") +
    theme(plot.title = element_text(hjust=0.5)) +
    theme(legend.position = "top", legend.text = element_text(size = 5)) +
    guides(color = "none") +

    merged_seurat_f@meta.data %>%
    ggplot(aes(x=log10GenesPerUMI, color = !!as.name(grouping.var), fill=!!as.name(grouping.var))) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    geom_vline(xintercept = 0.8) +
    ggtitle("Complexity (after filtering)") +
    theme(plot.title = element_text(hjust=0.5)) +
    theme(legend.position = "top", legend.text = element_text(size = 5)) +
    guides(color = "none")
  if(qc.plots == T) {
    print(b)
  }
  return(merged_seurat_f)
}
