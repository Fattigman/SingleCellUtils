library(Seurat)
library(patchwork)
library(sctransform)
library(AUCell)
library(ggplot2)
library(glue)
library(dplyr)

#' Filter Differentially Expressed Genes
#'
#' This function filters differentially expressed genes based on average log2 fold change and cluster information. Very important!
#'
#' @param markers A data frame containing marker genes and their associated statistics.
#' @param seurat_object A Seurat object containing single-cell data.
#' @param n_genes An integer specifying the number of top genes to select per cluster.
#' @param p_filter A numeric specifying the p-value threshold for filtering.
#' @return A list of differentially expressed genes per cluster.
#' @export 
FilterDegs <- function(markers, seurat_object, n_genes = 6, p_filter = 0.05){
    # Initialize an empty list for storing DEGs (differentially expressed genes)
    deg_list <- list()
    if ("cluster" %in% names(markers)){
        clusters <- levels(seurat_object)
        # Iterate through clusters (0 to n)
        flag <- 0
        for (i in clusters) {
            # Filter markers based on avg_log2FC and cluster
            degs <- markers %>%
                group_by(cluster) %>%
                filter(abs(avg_log2FC) > 1, cluster == i) %>%
                arrange(p_val_adj) %>%
                filter(p_val_adj < p_filter) %>%
                slice_head(n = n_genes) %>%
                pull(gene)
            # Check if any DEGs were found
            degs = sort(degs, decreasing = TRUE)
            deg_list[[flag+1]] <- degs
            flag <- flag + 1
        }
    }else {
        
        degs <- markers %>%
            filter(abs(avg_log2FC) > 1) %>%
            arrange(p_val_adj) %>%
            filter(p_val_adj < p_filter) %>%
            slice_head(n = n_genes)
        
        degs <- rownames(degs)
        # Check if any DEGs were found
        deg_list[[1]] <- degs
    }
    
    return(deg_list)
}

#' Cluster Composition Analysis
#'
#' This function calculates the composition of each cluster in terms of original identities.
#'
#' @param seurat_object A Seurat object containing single-cell data.
#' @param meta_data A string specifying the metadata variable for original identities.
#' @param idents A vector specifying the cluster identities to include in the analysis.
#' @param cluster_name A string specifying the name of the cluster identity column in the Seurat object.
#' @param reduced A logical specifying whether each orig.ident should be reduced to the smallest cell count.
#' @return A plot showing the cluster composition in terms of original identities.
#' @export 
ClusterComp <-  function(seurat_object, meta_data= "orig.ident", idents = NULL, cluster_name = "seurat_clusters", reduced = FALSE){
    if (reduced) {
        least_cells_sample <- rownames(as.matrix(table(seurat_object$orig.ident)[table(seurat_object$orig.ident) == min(table(seurat_object$orig.ident))]))
        least_cells <- table(seurat_object$orig.ident)[table(seurat_object$orig.ident) == min(table(seurat_object$orig.ident))][1]
        cell_list <- c()
        for (sample in unique(seurat_object$orig.ident)) {
            # Create a subset of the data for the current sample
            if (sample == least_cells_sample){
                next
            }
            subset_sample <- seurat_object[,seurat_object$orig.ident == sample]
            
            # Sample 6413 column names from the subset
            sampled_columns <- sample(colnames(subset_sample), least_cells)
            
            # Append the sampled columns to the cell_list
            cell_list <- append(cell_list, sampled_columns)
        }
        gc()
        seurat_object <- seurat_object[,cell_list]
    }
    # Get the cluster identities and original identities
    cluster_ids <- seurat_object[[cluster_name]]
    orig_ids <- seurat_object@meta.data[[meta_data]]
    
    # Create a data frame
    df <- data.frame(Cluster = cluster_ids, OrigIdent = orig_ids)
    
    # Count the number of cells from each orig.ident in each cluster
    df_counts <- table(df)
    
    # Convert the table to a data frame
    df_counts <- as.data.frame(df_counts)
    
    # Rename the columns
    colnames(df_counts) <- c("Cluster", "OrigIdent", "Count")
    
    # Filter out cells if idents is specified
    if (!is.null(idents)) {
        df_counts <- df_counts %>%
            filter(Cluster %in% idents)
    }
    
    # Create the fraction column
    df_counts <- df_counts %>%
        group_by(Cluster) %>%
        mutate(TotalCount = sum(Count), Fraction = Count / TotalCount) %>%
        ungroup() # Ungroup to avoid issues in further data manipulations
    
    # Create the count plot
    count_plot <- ggplot(df_counts, aes(x = Cluster, y = Count, fill = OrigIdent)) +
        geom_bar(stat = "identity") +
        theme(legend.position = "none", axis.text = element_text(angle = 90)) +
        labs(x = "Cluster", y = "Cluster composition (Cell Count)", fill = "Sample", angle = 90,title = glue("{meta_data} Cluster composition (Cell Count)"))
    
    # Annotate the bars with counts
    count_plot <- count_plot
    
    # Create the fraction plot
    fraction_plot <- ggplot(df_counts, aes(x = Cluster, y = Fraction, fill = OrigIdent)) +
        geom_bar(stat = "identity") +
        theme(axis.text = element_text(angle = 90)) +
        labs(x = "Cluster", y = "Cluster Composition (%)", fill = "Sample", title = glue("{meta_data} cluster composition (%)"))
    
    count_plot + fraction_plot
}

#' Gene Ontology and KEGG Pathway Analysis
#'
#' This function performs gene ontology and KEGG pathway analysis on differentially expressed genes.
#'
#' @param seurat_object A Seurat object containing single-cell data.
#' @return A list of dot plots showing the enriched gene ontology and KEGG pathways.
#' @export 
GO_KEGG <- function(seurat_object){
    
    
    deg.ls <- split(rownames(seurat_object), f = seurat_object$seurat_clusters)
    
    geneid.ls <- deg.ls %>% map(~{
        
        # here for macaque
        gene.df <- select(org.Hs.eg.db,
                                            keys = .x,
                                            columns = c("ENTREZID", "SYMBOL"),
                                            keytype = "SYMBOL")
        
        gene <- gene.df$ENTREZID
        gene <- gene[which(!is.na(gene))]
        gene <- unique(gene)
        
    })
    
    gene.ls <- geneid.ls
    
    # her mcc for macaque
    compKEGG <- compareCluster(geneCluster   = gene.ls,
                                                         fun           = "enrichKEGG",
                                                         pvalueCutoff  = 0.05,
                                                         pAdjustMethod = "BH", 
                                                         organism = "hsa")
    
    compGO <- compareCluster(geneCluster   = gene.ls,
                                                     fun           = "enrichGO",
                                                     pvalueCutoff  = 0.05,
                                                     pAdjustMethod = "BH", 
                                                     OrgDb = org.Hs.eg.db, 
                                                     ont = 'BP')
    
    #  compPathway <- compareCluster(geneCluster   = gene.ls,
    #                                fun           = "enrichPathway",
    #                                pvalueCutoff  = 0.05,
    #                                pAdjustMethod = "BH")
    
    ## dot plot
    g1 <- dotplot(compGO, showCategory = 10, title = "GO Enrichment Analysis")
    #g2 <- dotplot(compPathway, showCategory = 10, title = "REACTOME Pathway Enrichment Analysis")
    g3 <- dotplot(compKEGG, showCategory = 10, title = "KEGG Pathway Enrichment Analysis")
    plots <- c()
    plots[[1]] <- g1
    plots[[2]] <- g3
    return (plots)
}

#' Create Merged Seurat Object
#'
#' This function creates a merged Seurat object from multiple samples.
#'
#' @param file_paths A character vector specifying the file paths of the single-cell data.
#' @return A merged Seurat object.
#' @export 
CreateSeuratList <- function(file_paths){
    # Initialize an empty list to store the Seurat objects for each sample
    seurat_list <- list()
    
    for (i in seq_along(file_paths)){
        seurat_data <- Read10X(data.dir = file_paths[[i]])
        sample_name <- strsplit(file_paths[[i]],"/")[[1]][[3]]
        seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                                                         min.features = 100, 
                                                                         min.cells = 3,
                                                                         project = sample_name)
        seurat_list[[i]] <- seurat_obj
    }
    
    merged_seurat <- merge(x = seurat_list[[1]], 
                                                 y = seurat_list[-1])
    
    # Remove objects to clear up memory
    rm(seurat_data)
    rm(seurat_list)
    rm(seurat_obj)
    
    # Compute percent mito ratio
    merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
    merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100
    
    
    # Add number of genes per UMI for each cell to metadata
    merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)
    
    return (merged_seurat)
}

#' Apply DoubletFinder
#'
#' This function applies DoubletFinder to identify potential doublets in the single-cell data.
#'
#' @param seurat_object A Seurat object containing single-cell data.
#' @return A modified Seurat object with doublet information.
#' @export 
ApplyDoubletFinder <- function (seurat_object){
    
    seurat_object$percent.mito <- PercentageFeatureSet(object = seurat_object, pattern = "^MT-")
    
    # run sctransform
    seurat_object <- SCTransform(seurat_object, vars.to.regress="percent.mito", assay="RNA",do.correct.umi=T,conserve.memory=T, verbose = T)
    
    seurat_object <- RunPCA(seurat_object)
    seurat_object <- RunUMAP(seurat_object, dims = 1:15)
    
    seurat_object <- FindNeighbors(seurat_object, dims = 1:15)
    
    seurat_object <- FindClusters(seurat_object, resolution = 0.3)
    
    sweep.res.list_seurat <- paramSweep(seurat_object, PCs = 1:15, sct = TRUE)
    sweep.stats_seurat <- summarizeSweep(sweep.res.list_seurat, GT = FALSE)
    bcmvn_seurat <- find.pK(sweep.stats_seurat)
    
    df <- data.frame(
        `cells.recovered` = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
        `multiplet.rate` = c(0.004, 0.008, 0.016, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076)
    )
    
    # 1. Fetch amount of cells from a seurat object
    num_cells <- length(seurat_object@meta.data$orig.ident)
    
    # 2. Check which row in the "# of Cells Recovered" column which this value is closest to
    closest_row_index <- which.min(abs(df$cells.recovered - num_cells))
    
    # 3. Define a new variable that matches the corresponding value on "Multiplet Rate (%)", that is on the same row.
    multiplet_rate <- df$multiplet.rate[closest_row_index]
    
    homotypic.prop <- modelHomotypic(seurat_object@meta.data$RNA_snn_res.0.3)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
    nExp_poi <- round(multiplet_rate*nrow(seurat_object@meta.data))
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    
    seurat_object <- doubletFinder(seurat_object, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
    
    return (seurat_object)
}

#' Cell Counts Analysis
#'
#' This function calculates the number of cells per metadata variable and creates a bar plot.
#'
#' @param seurat_object A Seurat object containing single-cell data.
#' @param meta_data_var A string specifying the metadata variable to analyze.
#' @return A bar plot showing the cell count per metadata variable.
#' @export 
CellCounts <- function(seurat_object, meta_data_var = "orig.ident") {
    # Extract the metadata from the Seurat object
    metadata <- seurat_object@meta.data
    
    # Count the number of cells per metadata var
    cell_counts <- metadata %>%
        group_by(.data[[meta_data_var]]) %>%
        summarise(n = n())
    
    # Create the bar plot
    p <- ggplot(cell_counts, aes(x = .data[[meta_data_var]], y = n)) +
        geom_bar(stat = "identity") +
        theme(legend.position = "none", axis.text = element_text(angle = 90)) +
        labs(x = meta_data_var, y = "Cell count", title = glue("Cell count per {meta_data_var}"))
    
    # Add text labels above each bar
    p + geom_text(aes(label = n), vjust = -0.5)
}

#' SCT Normalize Seurat Object
#'
#' This function performs SCT normalization on a Seurat object.
#'
#' @param seurat_object A Seurat object containing single-cell data.
#' @param vars.to.regress A character vector specifying the variables to regress out during normalization.
#' @param n.features An integer specifying the number of variable features to select.
#' @return A normalized Seurat object.
#' @export 
SCTNormalize <- function(seurat_object, vars.to.regress = NULL, n.features = 3000){
    # Get genes that aren't Mitochondrial or chrY
    gene.table <- read.table("workspace/data/raw_external/hg38_gencode_v27.txt")
    non.sex.genes <- gene.table %>% 
        filter(!V2 %in% c("chrM", "chrY")) %>%
        pull(V1)
    
    
    if (!is.null(vars.to.regress)){
        print(vars.to.regress)
        seurat_object <- SCTransform(seurat_object, do.correct.umi=T,conserve.memory=T, vars.to.regress = vars.to.regress, variable.features.n = n.features, )
    }
    else(
        seurat_object <- SCTransform(seurat_object, do.correct.umi=T,conserve.memory=T, variable.features.n = n.features)
    )
    
    print("SCTransform done")
    gc()
    # Set 1s to 0s 
    print("Joining layers")
    xr <- JoinLayers(seurat_object, assay = "RNA")[["RNA"]]$counts; gc()
    #xr <- seurat_object[["RNA"]]$counts; gc()
    wr<-xr==0
    rm(xr);gc()
    print("Creating wr matrix")
    wr<-wr[rownames(seurat_object[["SCT"]]@counts),];gc()
    
    print("Creating xs matrix")
    xs<-as.matrix(seurat_object[["SCT"]]@counts);gc()
    
    print("Reset true 0s in xs")
    for (i in seq(1, ncol(xs), by = 1000)){
        xs[, i:min(i+999, ncol(xs))][as.matrix(wr[, i:min(i+999, ncol(xs))])] <- 0
        gc()
    }
    rm(wr);gc()
    print("Reset true 0s in seurat object SCT counts")
    seurat_object[["SCT"]]@counts<-as.sparse(xs)
    rm(xs);gc()
    print("Overwrite SCT data slot in seurat object")
    #overwrite seurat_object slot by log2(seurat_object+1) instead of loge(seurat_object+1) from sctransform
    seurat_object[["SCT"]]@data<-as.sparse(log2(as.matrix(seurat_object[["SCT"]]@counts)+1))
    
    seurat_object <- ScaleData(seurat_object)
    
    genes <- c()
    for (i in seq_along(seurat_object[['SCT']]@var.features)){
        if (seurat_object[['SCT']]@var.features[i] %in% non.sex.genes){
            len <- length(genes)
            genes[len + 1] <- seurat_object[['SCT']]@var.features[i]
        }
    }
    seurat_object <- RunPCA(seurat_object, features = genes)
    seurat_object <- RunUMAP(seurat_object, dims = 1:30)
    return(seurat_object)
}

#' Calculate MetaGene Expression
#'
#' This function calculates the average expression of a gene list across cells.
#'
#' @param seurat_object A Seurat object containing single-cell data.
#' @param gene_list A character vector specifying the gene list to calculate the average expression.
#' @return A numeric vector representing the average expression of the gene list.
#' @export 
CalcMetaGene <- function(seurat_object, gene_list){
    return(colMeans(x = seurat_object@assays$SCT$data[gene_list, ], na.rm = TRUE))
}

#' Visualize Summarized Loadings
#'
#' This function visualizes the aggregated loadings across principal components.
#'
#' @param seurat_object A Seurat object containing single-cell data.
#' @param dims An integer specifying the number of principal components to consider.
#' @param genes An integer specifying the number of top genes to show in the plot.
#' @return A bar plot showing the aggregated loadings across principal components.
#' @export 
VizSummarizedLoadings <- function(seurat_object, dims = 15, genes = 20){
  loadings_matrix <- seurat_object@reductions$pca@feature.loadings[,1:dims]
  
  # Calculate absolute loadings across all PCs
  abs_loadings <- abs(loadings_matrix)
  
  # Sum up the absolute loadings for each feature
  aggregated_loadings <- rowSums(abs_loadings)
  
  # Rank features based on aggregated loadings
  sorted_features <- names(sort(aggregated_loadings, decreasing = TRUE))[1:genes]
  # Visualize aggregated loadings (example using ggplot2)
  
  df <- data.frame(Feature = sorted_features, Aggregated_Loadings = aggregated_loadings[sorted_features])
  df <- df[order(-df$Aggregated_Loadings), ]
  ggplot(df, aes(x = Feature, y = Aggregated_Loadings)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Features", y = "Aggregated Loadings", title = "Feature Importance Across PCs")
  
}
#' Calculate AUCell scores and add to Seurat object
#'
#' This function visualizes the aggregated loadings across principal components.
#'
#' @param seurat_object A Seurat object containing single-cell data.
#' @param gene_list A character vector specifying the gene list to calculate the module score.
#' @param col.name A string specifying the column name to store the AUCell scores.
#' @return A Seurat object with AUCell scores added as metadata.
#' @export 
CalculateAUCell <- function(seurat_object, gene_list, col.name){
  if (options()$Seurat.object.assay.version != 'v3'){
    seurat.option = options()$Seurat.object.assay.version
    options(Seurat.object.assay.version = 'v3')
  }else{
    seurat.option = 'v3'
  }
  gex_matrix <- seurat_object@assays$RNA@counts
  
  cells_rankings <- AUCell_buildRankings(gex_matrix, plotStats=FALSE, col.name)
  
  geneSets <- list(geneSet1=gene_list)
  
  cells_AUC <- AUCell_run(gex_matrix, geneSets, aucMaxRank=nrow(cells_rankings)*0.05)
  
  seurat_object <- AddMetaData(seurat_object, t(as.data.frame(cells_AUC@assays@data$AUC)), col.name = col.name)
  
  options(Seurat.object.assay.version = seurat.option)
  
  return(seurat_object)
}
#' Visualize Summarized Loadings
#'
#' This function visualizes the aggregated loadings across principal components.
#'
#' @param seurat_object A Seurat object containing single-cell data.
#' @param cell_names A character vector specifying the cell names to annotate.
#' @param dimplot_ident A string specifying the metadata variable to group by in the DimPlot.
#' @return A Seurat object with cell names annotated.
#' @export 
AnnotateCell <- function(seurat_object, cell_names, group_ident = "cell_types", dimplot_ident = "orig.ident"){
  plot <- DimPlot(seurat_object, group.by = dimplot_ident)
  selected.cells <- CellSelector(plot, seurat_object)
  cells <- WhichCells(object = selected.cells, idents = "SelectedCells")
  cell.df <- data.frame(cells, cell_names, row.names = cells)
  cell.df$cells <- NULL
  seurat_object <- AddMetaData(seurat_object, cell.df, col.name = group_ident)
}
