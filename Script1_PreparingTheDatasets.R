# Before running this script, you must create a working directory containing the following folders:
# 1) "Data", containing
    # the filtered_feature_bc_matrix folder from cell ranger for each library (Simon et al., in preparation)
    # the neural network classifier files (Apendix 1 from Ozel et al. 2021, https://pubmed.ncbi.nlm.nih.gov/33149298/)
# 2) "Results", which will contain the matrices and Seurat objects resulting from the analysis
# 3) "Plots"

#000000000000000000000#
#### Packages used ####
#000000000000000000000#

library(Seurat)
library(ggplot2)
library(cowplot)

#00000000000000000000000000000000000000000000000000#
#### Paths, to change before running the script ####
#00000000000000000000000000000000000000000000000000#

# Working directory
    path_wd = "..."
    
# Libraries
    path_Optix_Lib1 = "..."
    path_Optix_Lib2 = "..."
    path_Optix_Lib3 = "..."

    path_vOptix_Lib1 = "..."
    path_vOptix_Lib2 = "..."

    path_dOptix_Lib1 = "..."
    path_dOptix_Lib2 = "..."
    
    path_hh_Lib1 = "..."
    path_hh_Lib2 = "..."

    path_dpp_Lib1 = "..."
    path_dpp_Lib2 = "..."
  
    path_pxb = "..."

# The markers necessary to classify the cells.
# It is important to choose the right stage of development.
# These are in the neural network classifier files
    path_markers_Adult = "..."
    path_markers_P15 = "..."

#000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000#        
#### Parameters, to change before running the script, or when adding or removing libraries ####
#000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000#
    
# Dimensionality of the dataset (used to produce 2D visualizations)
    dims = 150
    
# Parameters used to create the libraries
    # "Include features detected in at least this many cells"
        min_cells = 3
    # "Include cells where at least this many features are detected"
        min_features = 200

# Parameters for quality control filtering of the different libraries
    percent.mt_Optix_Lib1 = 7
    nCount_RNA_Optix_Lib1 = 17000
    nFeature_RNA_Optix_Lib1 = 800
    
    percent.mt_Optix_Lib2 = 7
    nCount_RNA_Optix_Lib2 = 17000
    nFeature_RNA_Optix_Lib2 = 800
    
    percent.mt_Optix_Lib3 = 7
    nCount_RNA_Optix_Lib3 = 17000
    nFeature_RNA_Optix_Lib3 = 800
    
    percent.mt_vOptix_Lib1 = 10
    nCount_RNA_vOptix_Lib1 = 10000
    nFeature_RNA_vOptix_Lib1 = 500
    
    percent.mt_vOptix_Lib2 = 10
    nCount_RNA_vOptix_Lib2 = 10000
    nFeature_RNA_vOptix_Lib2 = 500
    
    percent.mt_dOptix_Lib1 = 5
    nCount_RNA_dOptix_Lib1 = 20000
    nFeature_RNA_dOptix_Lib1 = 1000
    
    percent.mt_dOptix_Lib2 = 5
    nCount_RNA_dOptix_Lib2 = 20000
    nFeature_RNA_dOptix_Lib2 = 1000
    
    percent.mt_hh_Lib1 = 10
    nCount_RNA_hh_Lib1 = 20000
    nFeature_RNA_hh_Lib1 = 700
    
    percent.mt_hh_Lib2 = 10
    nCount_RNA_hh_Lib2 = 20000
    nFeature_RNA_hh_Lib2 = 700
    
    percent.mt_dpp_Lib1 = 5
    nCount_RNA_dpp_Lib1 = 20000
    nFeature_RNA_dpp_Lib1 = 900
    
    percent.mt_dpp_Lib2 = 5
    nCount_RNA_dpp_Lib2 = 20000
    nFeature_RNA_dpp_Lib2 = 900
    
    percent.mt_pxb = 5
    nCount_RNA_pxb = 30000
    nFeature_RNA_pxb = 1300
    
# List of all filtering parameters, it is important to keep the same order in each vector
    percent.mt_list = c(percent.mt_Optix_Lib1, percent.mt_Optix_Lib2, percent.mt_Optix_Lib3,
                        percent.mt_vOptix_Lib1, percent.mt_vOptix_Lib2,
                        percent.mt_dOptix_Lib1, percent.mt_dOptix_Lib2,
                        percent.mt_hh_Lib1, percent.mt_hh_Lib2,
                        percent.mt_dpp_Lib1, percent.mt_dpp_Lib2,
                        percent.mt_pxb)            
    nCount_RNA_list = c(nCount_RNA_Optix_Lib1, nCount_RNA_Optix_Lib2, nCount_RNA_Optix_Lib3,
                        nCount_RNA_vOptix_Lib1, nCount_RNA_vOptix_Lib2,
                        nCount_RNA_dOptix_Lib1, nCount_RNA_dOptix_Lib2,
                        nCount_RNA_hh_Lib1, nCount_RNA_hh_Lib2,
                        nCount_RNA_dpp_Lib1, nCount_RNA_dpp_Lib2,
                        nCount_RNA_pxb)
    nFeature_RNA_list = c(nFeature_RNA_Optix_Lib1, nFeature_RNA_Optix_Lib2, nFeature_RNA_Optix_Lib3,
                          nFeature_RNA_vOptix_Lib1, nFeature_RNA_vOptix_Lib2,
                          nFeature_RNA_dOptix_Lib1, nFeature_RNA_dOptix_Lib2,
                          nFeature_RNA_hh_Lib1, nFeature_RNA_hh_Lib2,
                          nFeature_RNA_dpp_Lib1, nFeature_RNA_dpp_Lib2,
                          nFeature_RNA_pxb)

#0000000000000000000000000000000000000000000000000000000000000000000000000000000#
#### Creates the Seurat objects, to change when adding or removing libraries ####
#0000000000000000000000000000000000000000000000000000000000000000000000000000000#
    
setwd(dir = path_wd)
    
Optix_Lib1 = Read10X(data.dir = path_Optix_Lib1)
Optix_Lib1 = CreateSeuratObject(counts = Optix_Lib1, project = "Optix_Lib1", min.cells = min_cells, min.features = min_features)
Optix_Lib2 = Read10X(data.dir = path_Optix_Lib2)
Optix_Lib2 = CreateSeuratObject(counts = Optix_Lib2, project = "Optix_Lib2", min.cells = min_cells, min.features = min_features)
Optix_Lib3 = Read10X(data.dir = path_Optix_Lib3)
Optix_Lib3 = CreateSeuratObject(counts = Optix_Lib3, project = "Optix_Lib3", min.cells = min_cells, min.features = min_features)

vOptix_Lib1 = Read10X(data.dir = path_vOptix_Lib1)
vOptix_Lib1 = CreateSeuratObject(counts = vOptix_Lib1, project = "vOptix_Lib1", min.cells = min_cells, min.features = min_features)
vOptix_Lib2 = Read10X(data.dir = path_vOptix_Lib2)
vOptix_Lib2 = CreateSeuratObject(counts = vOptix_Lib2, project = "vOptix_Lib2", min.cells = min_cells, min.features = min_features)

dOptix_Lib1 = Read10X(data.dir = path_dOptix_Lib1)
dOptix_Lib1 = CreateSeuratObject(counts = dOptix_Lib1, project = "dOptix_Lib1", min.cells = min_cells, min.features = min_features)
dOptix_Lib2 = Read10X(data.dir = path_dOptix_Lib2)
dOptix_Lib2 = CreateSeuratObject(counts = dOptix_Lib2, project = "dOptix_Lib2", min.cells = min_cells, min.features = min_features)

hh_Lib1 = Read10X(data.dir = path_hh_Lib1)
hh_Lib1 = CreateSeuratObject(counts = hh_Lib1, project = "hh_Lib1", min.cells = min_cells, min.features = min_features)
hh_Lib2 = Read10X(data.dir = path_hh_Lib2)
hh_Lib2 = CreateSeuratObject(counts = hh_Lib2, project = "hh_Lib2", min.cells = min_cells, min.features = min_features)

dpp_Lib1 = Read10X(data.dir = path_dpp_Lib1)
dpp_Lib1 = CreateSeuratObject(counts = dpp_Lib1, project = "dpp_Lib1", min.cells = min_cells, min.features = min_features)
dpp_Lib2 = Read10X(data.dir = path_dpp_Lib2)
dpp_Lib2 = CreateSeuratObject(counts = dpp_Lib2, project = "dpp_Lib2", min.cells = min_cells, min.features = min_features)

pxb = Read10X(data.dir = path_pxb)
pxb = CreateSeuratObject(counts = pxb, project = "pxb", min.cells = min_cells, min.features = min_features)

#000000000000000000000000000000000000000000000000000000000000000000000000000000000#
#### Stores the objects in lists, to change when adding or removing libraries  #### 
#000000000000000000000000000000000000000000000000000000000000000000000000000000000#

# List of the Seurat Objects
# Must be in the same order as percent.mt_list, nCount_RNA_list, nFeature_RNA_list
    Seurat_objects_list = c(Optix_Lib1, Optix_Lib2, Optix_Lib3,
                            vOptix_Lib1, vOptix_Lib2,
                            dOptix_Lib1, dOptix_Lib2,
                            hh_Lib1, hh_Lib2,
                            dpp_Lib1, dpp_Lib2,
                            pxb)
    
# List of library names, must be in the same order as percent.mt_list, nCount_RNA_list, nFeature_RNA_list
    Library_names = c("Optix_Lib1", "Optix_Lib2", "Optix_Lib3",
                      "vOptix_Lib1", "vOptix_Lib2",
                      "dOptix_Lib1", "dOptix_Lib2",
                      "hh_Lib1", "hh_Lib2",
                      "dpp_Lib1", "dpp_Lib2",
                      "pxb")
    
# Same, but the libraries are separated by the stage at which they were acquired
    Library_names_Adult = c("Optix_Lib1", "Optix_Lib2", "Optix_Lib3",
                            "vOptix_Lib1", "vOptix_Lib2",
                            "hh_Lib1", "hh_Lib2")
    Library_names_P15 = c("dOptix_Lib1", "dOptix_Lib2",
                          "dpp_Lib1", "dpp_Lib2",
                          "pxb")

#000000000000000000000000000000000000000000000000000000000000000000000000000000#
#### Filters, normalizes, and runs TSNE, UMAP and PCA on the Seurat Objects ####
#000000000000000000000000000000000000000000000000000000000000000000000000000000#

# Creates a directory to store the Seurat Objects
    dir.create(path = "Results/Seurat_Objects/", showWarnings = TRUE, recursive = FALSE)
  
k = 0
for (Seurat_Object in Seurat_objects_list) {
    k = k + 1
    # Important that all these are in the same order, since they all use the same "k"
        # the library
        libraryname = Library_names[k]
        # the parameters to filter it
        nCount_RNA_Seurat_Object = nCount_RNA_list[k]
        nFeature_RNA_Seurat_Object = nFeature_RNA_list[k]
        percent.mt_Seurat_Object = percent.mt_list[k]
    # Creates a metadata field with mitochondrial genes
        Seurat_Object@meta.data[,"percent.mt"] <- PercentageFeatureSet(Seurat_Object, pattern = "^mt")
    # Visualizes the result of the quality control metrics chosen
        # Adds one column to the metadata, identifying the cells that will be discarded (value = 0) or kept (value = 1)
            Seurat_Object@meta.data$cells_kept = 1
            Seurat_Object@meta.data$cells_kept[Seurat_Object@meta.data$nCount_RNA > nCount_RNA_Seurat_Object | Seurat_Object@meta.data$nFeature_RNA < nFeature_RNA_Seurat_Object | Seurat_Object@meta.data$percent.mt > percent.mt_Seurat_Object] = 0
            n_cells_before = length(Seurat_Object@meta.data$cells_kept)
            n_cells_after = sum(Seurat_Object@meta.data$cells_kept)
        # Creates violin plots 
            a = VlnPlot(Seurat_Object, features = "nFeature_RNA", ncol = 1, pt.size = 0) + 
                geom_hline(yintercept=nFeature_RNA_Seurat_Object, linetype="dashed", color = "red", lwd=2) + 
                labs(title = "Gene nb") +
                theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5)) +
                NoLegend()
            b = VlnPlot(Seurat_Object, features = "nCount_RNA", pt.size = 0) + 
                geom_hline(yintercept=nCount_RNA_Seurat_Object, linetype="dashed", color = "red", lwd=2) + 
                labs(title = "RNA nb") +
                theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5)) +
                NoLegend()
            c = VlnPlot(Seurat_Object, features = "percent.mt", pt.size = 0) +
                geom_hline(yintercept=percent.mt_Seurat_Object, linetype="dashed", color = "red", lwd=2) + 
                labs(title = "Mito %") +
                theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5)) +
                NoLegend()
            plots = plot_grid(a, b, c, align= "h", nrow = 1)
            pdf(file = paste("Plots/", libraryname, "_QC_VlnPlot_", n_cells_before , "_to_", n_cells_after , "_cells.pdf", sep = ""))
            print(plots)
            dev.off()
    # Does the filtering, normalization, scaling
        Seurat_Object = subset(Seurat_Object, subset = percent.mt <= percent.mt_Seurat_Object & nCount_RNA <= nCount_RNA_Seurat_Object & nFeature_RNA >= nFeature_RNA_Seurat_Object)
        Seurat_Object = NormalizeData(Seurat_Object, normalization.method = "LogNormalize", scale.factor = 10000)
        Seurat_Object = FindVariableFeatures(Seurat_Object, selection.method = "vst", nfeatures = 2000)
        Seurat_Object = ScaleData(Seurat_Object)
    # Produces tsne/umap for visualization 
        Seurat_Object = RunPCA(Seurat_Object, npcs = dims)
        Seurat_Object = RunTSNE(Seurat_Object, reduction = "pca", dims = 1:dims)
        Seurat_Object = RunUMAP(Seurat_Object, reduction = "pca", dims = 1:dims)
    # Saves the filtered and normalized object
        saveRDS(object = Seurat_Object, file = paste("Results/Seurat_Objects/", libraryname, "_filtered_normalized.RDS", sep = ""))
}

#00000000000000000000000000000000000000000000000000000000000#
#### Creates the matrices for annotation of the datasets ####
#00000000000000000000000000000000000000000000000000000000000#

# Creates a directory to store the matrices for annotation by the neural network
    dir.create(path = "Results/Matrices_for_classification/", showWarnings = TRUE, recursive = FALSE)

# Loads the markers
    markers_Adult = readRDS(path_markers_Adult)
    markers_P15 = readRDS(path_markers_P15)

# Produces the matrices, for each stage independently
    for (libraryname in Library_names_Adult) {
        #The object to classify
            Object_filtered_normalized = readRDS(paste("Results/Seurat_Objects/", libraryname, "_filtered_normalized.RDS", sep = ""))
        # The data to classify    
            Matrix_for_classification = as.matrix(Object_filtered_normalized@assays$RNA@data[,])
        # Since all the markers necessary are not always found in the Seurat Object, this adds the missing markers with an expression value of 0
            missing_markers = setdiff(markers_Adult, rownames(Matrix_for_classification))
            To_add = matrix(data = 0, nrow = length(missing_markers), ncol = ncol(Matrix_for_classification))
            rownames(To_add) = missing_markers
            colnames(To_add) = colnames(Matrix_for_classification)
            Matrix_for_classification = rbind(Matrix_for_classification, To_add)
        # Keeps only the markers necessary for the classification, and puts them in the order required for the neural network
            if (length(setdiff(markers_Adult, rownames(Matrix_for_classification))) == 0) {
                Matrix_for_classification = Matrix_for_classification[markers_Adult, ]
                write.csv(x = Matrix_for_classification, file = paste("Results/Matrices_for_classification/", libraryname, "_NN_input.csv", sep = ""))
            }
    }
    
    for (libraryname in Library_names_P15) {
            Object_filtered_normalized = readRDS(paste("Results/Seurat_Objects/", libraryname, "_filtered_normalized.RDS", sep = ""))
            Matrix_for_classification = as.matrix(Object_filtered_normalized@assays$RNA@data[,])
            missing_markers = setdiff(markers_P15, rownames(Matrix_for_classification))
            To_add = matrix(data = 0, nrow = length(missing_markers), ncol = ncol(Matrix_for_classification))
            rownames(To_add) = missing_markers
            colnames(To_add) = colnames(Matrix_for_classification)
            Matrix_for_classification = rbind(Matrix_for_classification, To_add)
            if (length(setdiff(markers_P15, rownames(Matrix_for_classification))) == 0) {
                Matrix_for_classification = Matrix_for_classification[markers_P15, ]
                write.csv(x = Matrix_for_classification, file = paste("Results/Matrices_for_classification/", libraryname, "_NN_input.csv", sep = ""))
            }
    }
    
    