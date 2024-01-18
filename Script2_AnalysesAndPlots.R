# Before running this script, please classify the cells with the neural network classifier published in Ozel et al. 2021, https://pubmed.ncbi.nlm.nih.gov/33149298/

#000000000000000000000#
#### Packages used ####
#000000000000000000000#

library(Seurat)
library(readxl)
library(ggplot2)

#0000000000000000000000000000000000000000000#
#### To change before running the script ####
#0000000000000000000000000000000000000000000#

# Dataset that will be used as a reference to make the plots
    Reference_Stage = "Adult"

# Seurat object of the reference dataset, published in Ozel et al. 2021
    path_WholeOPC_Adult = "..."

# Working directory
    path_wd = "..."

# Dm clusters and their annotations
    # Ordered by abundance, including Dm15
    Dms_clusters =    c( "19", "128",  "35",  "225",  "226",   "15",  "124",   "12", "144",  "136",   "9",  "14")
    Dms_annotations = c("Dm8", "Dm2", "Dm3", "Dm3a", "Dm3b", "Dm10", "Dm15?", "Dm12", "Dm9", "Dm11", "Dm4", "Dm1")
    
# Clusters of immature cells:
    Imm_clusters =    c("240",  "241",  "242",  "243",  "244",  "245",  "246",  "247",  "248",  "249",  "250",   "251", "252",  "253",  "254",  "255",  "256", "257")
    Imm_annotations = c("Im01", "Im02", "Im03", "Im04", "Im05", "Im06", "Im07", "Im08", "Im09", "Im10", "Im11", "Im12", "Apop", "GMC1", "GMC2", "GMC3", "NB",  "LPC")
    
# Clusters used to normalize the data, and their annotation
    Norm_clusters =    c("1" , "142", "140", "118", "125", "62")
    Norm_annotations = c("T1", "Mi1", "Tm1", "Tm2", "Tm4", "Tm6")
    
# All files in the following lists should be in the same order
    # Paths to all files containing the predictions from the neural network
        paths_to_classification_prediction = c("...",
                                               "...",
                                               "...",
                                               "...",
                                               "...",
                                               "...",
                                               "...",
                                               "...",
                                               "...",
                                               "...",
                                               "...",
                                               "...")
    # Paths to all files containing the confidence in the predictions from the neural network
        paths_to_classification_confidence = c("...",
                                               "...",
                                               "...",
                                               "...",
                                               "...",
                                               "...",
                                               "...",
                                               "...",
                                               "...",
                                               "...",
                                               "...",
                                               "...")
    # Library names
        library_names = c("Optix_Lib1", "Optix_Lib2", "Optix_Lib3", "vOptix_Lib1", "vOptix_Lib2", "dOptix_Lib1", "dOptix_Lib2", "hh_Lib1", "hh_Lib2", "dpp_Lib1", "dpp_Lib2", "pxb")
    # Same sorted by stage
        library_names_Adult = c("Optix_Lib1", "Optix_Lib2", "Optix_Lib3", "vOptix_Lib1", "vOptix_Lib2", "hh_Lib1", "hh_Lib2")
        library_names_P15 = c("dOptix_Lib1", "dOptix_Lib2", "dpp_Lib1", "dpp_Lib2", "pxb")

# Parameters used to flag the classes (groups of cells with the same annotation) based on the confidence in their annotation
# If more than "confidence_proportion" of the cells in a class are classified with less than "confidence_threshold", the class will be flagged
    confidence_proportion = 0.8
    confidence_threshold = 0.5
    
# Each class containing this number of cells or less will also be flagged
    abundance_threshold = 3
    
# Order in which to plot the datasets on the graph showing the normalized abundances of each class in each dataset
    dataset_bargraph_order = c("Reference_dataset", "pxb_dataset", "Optix_Lib1_dataset", "Optix_Lib2_dataset", "Optix_Lib3_dataset", "vOptix_Lib1_dataset", "vOptix_Lib2_dataset", "dOptix_Lib1_dataset", "dOptix_Lib2_dataset", "dpp_Lib1_dataset", "dpp_Lib2_dataset", "hh_Lib1_dataset", "hh_Lib2_dataset")

#0000000000000000000000000000000000000000#
#### Adds predictions to the datasets ####
#0000000000000000000000000000000000000000#
    
setwd(path_wd)

k = 0
for (libraryname in library_names) {
    k = k + 1 
    Object_to_classify = readRDS(paste("Results/Seurat_Objects/", libraryname, "_filtered_normalized.RDS", sep = ""))
    # Creates a metadata entry with the neural network predictions
        Predictions = read.table(paths_to_classification_prediction[k])
        Object_to_classify@meta.data$NN_Cluster_Number = Predictions$V1
    # Creates a metadata entry with the confidence in the prediction for each cell
        Confidence = read.table(paths_to_classification_confidence[k])
        Object_to_classify@meta.data$Confidence_NN_Cluster_Number = Confidence$V1
    # Keeps only the Dm clusters, clusters of immature neurons, and clusters used to normalize
        Object_to_classify = subset(x = Object_to_classify, cells = colnames(Object_to_classify)[Object_to_classify@meta.data$NN_Cluster_Number %in% c(Dms_clusters, Norm_clusters, Imm_clusters)])
        # DimPlot(object = Object_to_classify, group.by = "NN_Cluster_Number")
    # Converts cluster numbers to annotations, and adds the annotations in a metadata field
    # Makes the annotations the active identities in the object
        Annotation = match(as.character(Object_to_classify@meta.data$NN_Cluster_Number), c(Dms_clusters, Norm_clusters, Imm_clusters))
        Annotation = c(Dms_annotations, Norm_annotations, Imm_annotations)[Annotation]
        Object_to_classify@meta.data$Annotation = Annotation
        Idents(Object_to_classify) = Annotation
        # DimPlot(object = Object_to_classify)
    # Saves the object
        saveRDS(object = Object_to_classify, file = paste("Results/Seurat_Objects/", libraryname, "_annotated.RDS", sep = ""))
}

#000000000000000000000000000000000000000000000000000000000000000#
#### Flags the classes annotated with a low confidence score ####
#000000000000000000000000000000000000000000000000000000000000000#

# This is too get non-overlapping labels on the 2D visualizations
    options(ggrepel.max.overlaps = Inf)

# Makes a table with a color for each class
    library(RColorBrewer)
    color_palette = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
    Class_colors = matrix(data = 0, nrow = length(c(Dms_annotations, Norm_annotations, Imm_annotations)), ncol = 2)
    colnames(Class_colors) = c("Colors", "Annotation")
    Class_colors[, "Annotation"] = c(Dms_annotations, Norm_annotations, Imm_annotations)
    Class_colors[, "Colors"] = sample(color_palette, length(c(Dms_annotations, Norm_annotations, Imm_annotations)), replace = F)
    
for (libraryname in library_names){
    Object_to_filter = readRDS(paste("Results/Seurat_Objects/", libraryname, "_annotated.RDS", sep="" ))
    # Makes a vector containing the colors to use for plotting, with each color named according to the class it corresponds to
        Colors = match(as.character(levels(Idents(Object_to_filter))), Class_colors[, "Annotation"])
        Colors = Class_colors[Colors, "Colors"]
        names(Colors) = as.character(levels(Idents(Object_to_filter)))
    # Flags the low confidence annotations
        for (cluster in levels(Idents(Object_to_filter))) {
            # Cells belonging to the cluster
                index = which(Idents(Object_to_filter) == cluster)
                cells = colnames(Object_to_filter)[index]
            # If more than "confidence_proportion" of the cells in a class are classified with less than "confidence_threshold" they will be flagged
            # When flagged, the name of the cluster is changed to LC_cluster, and the metadata field "cells_kept" is changed to "0"
                if ((sum(Object_to_filter$Confidence_NN_Cluster_Number[cells] < confidence_threshold) / length(cells)) >= confidence_proportion) {
                    Object_to_filter = SetIdent(object = Object_to_filter, cells = cells, value = paste("LC", cluster, sep = "_"))
                    Object_to_filter@meta.data$Annotation = Object_to_filter@active.ident
                    Object_to_filter@meta.data[cells, "cells_kept"] = 0
                }
        }
    # Violin plots of the confidence in the annotation of each cell, grouped by class
        pdf(file = paste("Plots/", libraryname ,"_classification_confidence_filtered.pdf", sep = ""), width = 40, height = 5)
        a = VlnPlot(object = Object_to_filter, features = "Confidence_NN_Cluster_Number", cols = Colors) + NoLegend()
        print(a)
        dev.off()
    # TSNE and UMAP of the classified and filtered objects
        #tsne
        pdf(file = paste("Plots/", libraryname, "_annotated_tsne_filtered.pdf", sep = ""))
        print(DimPlot(object = Object_to_filter, group.by = "Annotation", reduction = "tsne", label = T, label.size = 2, shuffle = T, raster = F, cols = Colors) + NoLegend())
        dev.off()
        # With repelling labels
        pdf(file = paste("Plots/", libraryname, "_annotated_tsne_filtered_repel.pdf", sep = ""))
        print(DimPlot(object = Object_to_filter, group.by = "Annotation", reduction = "tsne", label = T, label.size = 2, shuffle = T, raster = F, repel = T, cols = Colors) + NoLegend())
        dev.off()
        # Only the Dm and Imm cells
        index = which(Object_to_filter@active.ident %in% c(Dms_annotations, Imm_annotations, paste("LC_", Dms_annotations, sep = ""), paste("LC_", Imm_annotations, sep = "")))
        pdf(file = paste("Plots/", libraryname, "_annotated_tsne_filtered_repel_Dm_Imm.pdf", sep = ""))
        print(DimPlot(object = Object_to_filter, group.by = "Annotation", reduction = "tsne", label = T, label.size = 3, shuffle = T, raster = F, repel = T, cols = Colors, cells = colnames(Object_to_filter)[index]))
        dev.off()
        # umap
        pdf(file = paste("Plots/", libraryname, "_annotated_umap_filtered.pdf", sep = ""))
        print(DimPlot(object = Object_to_filter, group.by = "Annotation", reduction = "umap", label = T, label.size = 2, shuffle = T, raster = F, cols = Colors) + NoLegend())
        dev.off()
        # With repelling labels
        pdf(file = paste("Plots/", libraryname, "_annotated_umap_filtered_repel.pdf", sep = ""))
        print(DimPlot(object = Object_to_filter, group.by = "Annotation", reduction = "umap", label = T, label.size = 2, shuffle = T, raster = F, repel = T, cols = Colors) + NoLegend())
        dev.off()
        # Only the Dm and Imm cells
        index = which(Object_to_filter@active.ident %in% c(Dms_annotations, Imm_annotations, paste("LC_", Dms_annotations, sep = ""), paste("LC_", Imm_annotations, sep = "")))
        pdf(file = paste("Plots/", libraryname, "_annotated_umap_filtered_repel_Dm_Imm.pdf", sep = ""))
        print(DimPlot(object = Object_to_filter, group.by = "Annotation", reduction = "umap", label = T, label.size = 3, shuffle = T, raster = F, repel = T, cols = Colors, cells = colnames(Object_to_filter)[index]))
        dev.off()
    # TSNE and UMAP with the low confidence annotations highlighted
        pdf(file = paste("Plots/", libraryname, "_annotated_tsne_low_confidence_cells.pdf", sep = ""))
        print(DimPlot(object = Object_to_filter, group.by = "cells_kept", reduction = "tsne", label = F, label.size = 2, raster = F, cells.highlight = colnames(Object_to_filter)[which(Object_to_filter@meta.data$cells_kept == "0")]) + NoLegend())
        dev.off()
        pdf(file = paste("Plots/", libraryname, "_annotated_umap_low_confidence_cells.pdf", sep = ""))
        print(DimPlot(object = Object_to_filter, group.by = "cells_kept", reduction = "umap", label = F, label.size = 2, raster = F, cells.highlight = colnames(Object_to_filter)[which(Object_to_filter@meta.data$cells_kept == "0")]) + NoLegend())
        dev.off()
    # Saves the object
        saveRDS(object = Object_to_filter, file = paste("Results/Seurat_Objects/", libraryname, "_annotated.RDS", sep = ""))
}

#000000000000000000000000000000000000000000000000000000000000000000000000#
#### Calculates normalized class abundances for the reference dataset ####
#000000000000000000000000000000000000000000000000000000000000000000000000#

# Loads and annotates the reference dataset
    if (Reference_Stage == "Adult") {
        Reference_dataset = readRDS(path_WholeOPC_Adult)
    }
    if (Reference_Stage == "P15") {
        Reference_dataset = readRDS(path_WholeOPC_P15)
    }
    
# Filtering and annotation of the dataset
    clusters_to_keep = c(Dms_clusters, Norm_clusters)[c(Dms_clusters, Norm_clusters) %in% Idents(Reference_dataset)]
    Reference_dataset = subset(x = Reference_dataset, idents = clusters_to_keep)
    Annotation = match(as.character(Idents(Reference_dataset)), c(Dms_clusters, Norm_clusters))
    Annotation = c(Dms_annotations, Norm_annotations)[Annotation]
    Idents(Reference_dataset) = Annotation
    # DimPlot(object = Reference_dataset, label = T) + NoLegend()
     
# Creates a table with the normalized frequency of each class in the dataset (average for all libraries)
# as well as the min and max normalized frequency among all the libraries of the dataset  
    Reference_dataset_to_plot = matrix(data = 0, nrow = length(table(Idents(Reference_dataset))), ncol = 5)
    colnames(Reference_dataset_to_plot) = c("Object", "Cluster", "Frequency", "Max_frequency", "Min_frequency")
    Reference_dataset_to_plot[, "Object"] = "Reference_dataset"
    Reference_dataset_to_plot[, "Cluster"] = names(table(Idents(Reference_dataset)))
    for (i in 1:nrow(Reference_dataset_to_plot)) {
      cluster = Reference_dataset_to_plot[i, "Cluster"]
      # Nb of cells belonging to the cluster in each library
          index = which(Idents(Reference_dataset) == cluster)
          cluster_abundances  = table(Reference_dataset$orig.ident[index])
      # Normalization
          normalizing_abundances = 0
          for (Normalizing_cluster in Norm_annotations) {
            # Number of cells from this cluster in each library of the reference dataset
            temp = table(Reference_dataset$orig.ident[which(Idents(Reference_dataset) == Normalizing_cluster)])
            # Makes sure the vector is in the same order as cluster_abundances
            temp = temp[names(cluster_abundances)]
            normalizing_abundances = normalizing_abundances + temp
          }
          normalizing_abundances = normalizing_abundances / length(Norm_annotations)
          Normalized_freq = cluster_abundances / normalizing_abundances
      # Fills the table
        Reference_dataset_to_plot[i, "Frequency"] = mean(Normalized_freq)
        Reference_dataset_to_plot[i, "Max_frequency"] = max(Normalized_freq)
        Reference_dataset_to_plot[i, "Min_frequency"] = min(Normalized_freq)
    }
    
# Orders clusters for plotting
    Reference_dataset_to_plot = Reference_dataset_to_plot[match(c(Dms_annotations, Norm_annotations)[c(Dms_annotations, Norm_annotations) %in% Reference_dataset_to_plot[, "Cluster"]], Reference_dataset_to_plot[, "Cluster"]), ]
    
#0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000#
#### Calculates normalized class abundances for each library, and flags low abundance and low confidence classes ####
#0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000#
    
k= 0
tab = c()
for (libraryname in library_names){
    k=k+1
    Lib = readRDS(paste("Results/Seurat_Objects/", libraryname, "_annotated.RDS", sep="" ))
    # Creates the table that will contain the normalized frequency of each class in the dataset
    # Here, "Frequency", "Max_frequency" and "Min_frequency" are all the same, because there is only 1 library, but these columns are necessary to fuse with the reference dataset table
        Lib_to_plot = matrix(data = 0, nrow = length(table(Idents(Lib))), ncol = 5)
        colnames(Lib_to_plot) = c("Object", "Cluster", "Frequency", "Max_frequency", "Min_frequency")
        Lib_to_plot[, "Object"] = paste(libraryname,"_dataset", sep = "")
        Lib_to_plot[, "Cluster"] = names(table(Idents(Lib)))
        Lib_to_plot[, "Frequency"] = table(Idents(Lib))
    # Adds a column which flags all classes containing less than "abundance_threshold" cells
        Low_abundance = matrix(data = 0, nrow = nrow(Lib_to_plot), ncol = 1)
        Lib_to_plot = cbind(Lib_to_plot, Low_abundance)
        colnames(Lib_to_plot)[ncol(Lib_to_plot)] = "Low_abundance"
        for (cluster in Lib_to_plot[, "Cluster"]) {
            if (table(Idents(Lib))[cluster] <= abundance_threshold) {
              Lib_to_plot[Lib_to_plot[, "Cluster"] == cluster, "Low_abundance"] = 1
            }
        }
    # For the class flagged as low confidence in annotation
    # Reverts their name (removes the "LC_") but flags them in a new column
        index = which(grepl('LC_', Lib_to_plot[, "Cluster"]))
        Low_confidence = matrix(data = 0, nrow = nrow(Lib_to_plot), ncol = 1)
        Low_confidence[index] = 1
        Lib_to_plot = cbind(Lib_to_plot, Low_confidence)
        colnames(Lib_to_plot)[ncol(Lib_to_plot)] = "Low_confidence"
        Lib_to_plot[index, "Cluster"] = sub(".*LC_", "", Lib_to_plot[index, "Cluster"])
    # Adds clusters present in the reference dataset, but not in this library
        missing_clusters = setdiff(levels(Idents(Reference_dataset)), Lib_to_plot[, "Cluster"])
        Lib_to_add = matrix(data = 0, nrow = length(missing_clusters), ncol = 7)
        colnames(Lib_to_add) = c("Object", "Cluster", "Frequency", "Max_frequency", "Min_frequency", "Low_abundance", "Low_confidence")
        Lib_to_add[, "Object"] = paste(libraryname,"_dataset", sep = "")
        Lib_to_add[, "Cluster"] = missing_clusters
        Lib_to_plot = rbind(Lib_to_plot, Lib_to_add)
    # Normalization
        normalizing_abundance = 0
        for (Normalizing_cluster in Norm_annotations) {
          normalizing_abundance = normalizing_abundance + as.numeric(Lib_to_plot[which(Lib_to_plot[, "Cluster"] == Normalizing_cluster), "Frequency"])
        }
        normalizing_abundance = normalizing_abundance / length(Norm_annotations)
        Lib_to_plot[, "Frequency"] = (as.numeric(Lib_to_plot[, "Frequency"]) / normalizing_abundance)
    # Fills max and min columns, which are redundant with "Frequency" but necessary for downstream analysis
        Lib_to_plot[, "Max_frequency"] = Lib_to_plot[, "Frequency"]
        Lib_to_plot[, "Min_frequency"] = Lib_to_plot[, "Frequency"]
    # Orders the clusters for plotting
        Lib_to_plot = Lib_to_plot[match(c(Dms_annotations, Norm_annotations)[c(Dms_annotations, Norm_annotations) %in% Lib_to_plot[, "Cluster"]], Lib_to_plot[, "Cluster"]), ]
    # Removes clusters that are in this library but not in the reference dataset
        Lib_to_plot = Lib_to_plot[Lib_to_plot[, "Cluster"] %in% Idents(Reference_dataset), ]      
    # Stores the results
        tab[[k]] <- Lib_to_plot
}

#000000000000000000000000000#
##### Plots the results #####
#000000000000000000000000000#

    # Binds all tables, by first adding the columns "Low_abundance" and "Low_confidence" into the "Reference_dataset_to_plot" table
        To_plot = Reference_dataset_to_plot
        To_plot = cbind(To_plot, matrix(data = 0, nrow = nrow(To_plot), ncol = 1))
        colnames(To_plot)[ncol(To_plot)] = "Low_abundance"
        To_plot = cbind(To_plot, matrix(data = 0, nrow = nrow(To_plot), ncol = 1))
        colnames(To_plot)[ncol(To_plot)] = "Low_confidence"
        for (table_nb in 1:length(tab)) {
            To_plot = rbind(To_plot, tab[[table_nb]])
        }
      
    # adds if both low abundance and confidence
        To_plot = cbind(To_plot, matrix(data = 0, nrow = nrow(To_plot), ncol = 1))
        colnames(To_plot)[ncol(To_plot)] = "Low_abundance_confidence"
        index = which(To_plot[, "Low_abundance"] == 1 & To_plot[, "Low_confidence"] == 1)
        To_plot[index, "Low_abundance"] = 0
        To_plot[index, "Low_confidence"] = 0
        To_plot[index, "Low_abundance_confidence"] = 1
        
    # Prepares plotting
        To_plot = as.data.frame(To_plot)
        To_plot$Frequency = as.numeric(To_plot$Frequency)
        To_plot$Max_frequency = as.numeric(To_plot$Max_frequency)
        To_plot$Min_frequency = as.numeric(To_plot$Min_frequency)
        To_plot$Cluster = factor(To_plot$Cluster, levels = unique(To_plot$Cluster))
        To_plot$Object = factor(x = To_plot$Object, levels = dataset_bargraph_order)
        To_plot$Low_abundance = as.numeric(To_plot$Low_abundance)
        To_plot$Low_confidence = as.numeric(To_plot$Low_confidence)
        To_plot$Low_abundance_confidence = as.numeric(To_plot$Low_abundance_confidence)
    # Saves the table to plot
        dir.create(path = "Results/Class_frequencies_table/", showWarnings = TRUE, recursive = FALSE)
        write.csv(x = To_plot, file = paste("Results/Class_frequencies_table/Class_frequencies_table", Reference_Stage, "reference", paste(Norm_annotations, collapse = "_"), "normalization.csv", sep = "_"))
        
        
    # Plot with all the libraries shown
        pdf(file = paste("PLOTS/Class_frequencies", Reference_Stage, "reference", paste(Norm_annotations, collapse = "_"), "normalization.pdf", sep = "_"), width = 30, height = 7)
        ggplot(To_plot, aes(x=Cluster, y=Frequency, fill=Object)) +
          geom_bar(position = position_dodge(width = 0.95), width = 0.95, stat="identity") +
          geom_errorbar(aes(ymin=Min_frequency, ymax=Max_frequency), width=0, size = 0.2, position=position_dodge(0.95)) +
          scale_fill_manual(values=c("grey", "#E69F00", "#0072B2", "#0072B2", "#0072B2", "#009E73", "#009E73", "#56B4E9", "#56B4E9","#F0E442", "#F0E442", "#CC79A7", "#CC79A7")) +
          geom_text(aes(label = ifelse(Low_abundance, "A", "")), position = position_dodge(width = 0.95), vjust = -.5, size = 4) +
          geom_text(aes(label = ifelse(Low_confidence, "C", "")), position = position_dodge(width = 0.95), vjust = -.5, size = 4) +
          geom_text(aes(label = ifelse(Low_abundance_confidence, "AC", "")), position = position_dodge(width = 0.95), vjust = -.5, size = 3) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 15),
                panel.background = element_rect(fill = "white", colour = "grey", size = 1, linetype = "solid"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                #axis.ticks.y = element_blank(),
                legend.key = element_rect(fill = "white"))
        dev.off()
        
        
    # Simplified plot, with averaged libraries
        # First, summarizes the flags
            To_plot$Flag = 0
            index = which(To_plot$Low_abundance == 1 | To_plot$Low_confidence == 1 | To_plot$Low_abundance_confidence == 1)
            To_plot$Flag[index] = 1
            To_plot = To_plot[, c("Object", "Cluster", "Frequency", "Max_frequency", "Min_frequency", "Low_abundance", "Low_confidence", "Low_abundance_confidence", "Flag")]
        # Then averages for the libraries
            To_plot_simple = To_plot[To_plot$Object == "Reference_dataset", c("Object", "Cluster", "Frequency", "Max_frequency", "Min_frequency", "Flag")]
            for (cluster in c(Dms_annotations, Norm_annotations)[c(Dms_annotations, Norm_annotations) %in% To_plot$Cluster]) {
                # Subsets only the rows regarding this class and the libraries of the Optix region
                    index  = which(To_plot$Cluster == cluster & To_plot$Object %in% c("Optix_Lib1_dataset", "Optix_Lib2_dataset", "Optix_Lib3_dataset"))
                    # To_plot[index, ]
                    # To_plot[index, "Frequency"]
                    # To_plot[index, "Flag"]
                # Summarizes the values
                    to_add = matrix(data = 0, nrow = 1, ncol = ncol(To_plot_simple))
                    colnames(to_add) = c("Object", "Cluster", "Frequency", "Max_frequency", "Min_frequency", "Flag")
                    to_add[, "Object"] = "Optix"
                    to_add[, "Cluster"] = cluster
                    to_add[, "Frequency"] = mean(To_plot[index, "Frequency"])
                    to_add[, "Max_frequency"] = max(To_plot[index, "Frequency"])
                    to_add[, "Min_frequency"] = min(To_plot[index, "Frequency"])
                    if (sum(To_plot[index, "Flag"]) > 0) {to_add[, "Flag"] = 1}
                    To_plot_simple = rbind(To_plot_simple, to_add)
                    
                # Same for dOptix
                    index  = which(To_plot$Cluster == cluster & To_plot$Object %in% c("dOptix_Lib1_dataset", "dOptix_Lib2_dataset"))
                    to_add = matrix(data = 0, nrow = 1, ncol = ncol(To_plot_simple))
                    colnames(to_add) = c("Object", "Cluster", "Frequency", "Max_frequency", "Min_frequency", "Flag")
                    to_add[, "Object"] = "dOptix"
                    to_add[, "Cluster"] = cluster
                    to_add[, "Frequency"] = mean(To_plot[index, "Frequency"])
                    to_add[, "Max_frequency"] = max(To_plot[index, "Frequency"])
                    to_add[, "Min_frequency"] = min(To_plot[index, "Frequency"])
                    if (sum(To_plot[index, "Flag"]) > 0) {to_add[, "Flag"] = 1}
                    To_plot_simple = rbind(To_plot_simple, to_add)
                    
                # Same for vOptix
                    index  = which(To_plot$Cluster == cluster & To_plot$Object %in% c("vOptix_Lib1_dataset", "vOptix_Lib2_dataset"))
                    to_add = matrix(data = 0, nrow = 1, ncol = ncol(To_plot_simple))
                    colnames(to_add) = c("Object", "Cluster", "Frequency", "Max_frequency", "Min_frequency", "Flag")
                    to_add[, "Object"] = "vOptix"
                    to_add[, "Cluster"] = cluster
                    to_add[, "Frequency"] = mean(To_plot[index, "Frequency"])
                    to_add[, "Max_frequency"] = max(To_plot[index, "Frequency"])
                    to_add[, "Min_frequency"] = min(To_plot[index, "Frequency"])
                    if (sum(To_plot[index, "Flag"]) > 0) {to_add[, "Flag"] = 1}
                    To_plot_simple = rbind(To_plot_simple, to_add)

                # Same for pxb
                    index  = which(To_plot$Cluster == cluster & To_plot$Object %in% c("pxb_dataset"))
                    to_add = matrix(data = 0, nrow = 1, ncol = ncol(To_plot_simple))
                    colnames(to_add) = c("Object", "Cluster", "Frequency", "Max_frequency", "Min_frequency", "Flag")
                    to_add[, "Object"] = "pxb"
                    to_add[, "Cluster"] = cluster
                    to_add[, "Frequency"] = mean(To_plot[index, "Frequency"])
                    to_add[, "Max_frequency"] = max(To_plot[index, "Frequency"])
                    to_add[, "Min_frequency"] = min(To_plot[index, "Frequency"])
                    if (sum(To_plot[index, "Flag"]) > 0) {to_add[, "Flag"] = 1}
                    To_plot_simple = rbind(To_plot_simple, to_add)                    
                    
                # Same for hh
                    index  = which(To_plot$Cluster == cluster & To_plot$Object %in% c("hh_Lib1_dataset", "hh_Lib2_dataset"))
                    to_add = matrix(data = 0, nrow = 1, ncol = ncol(To_plot_simple))
                    colnames(to_add) = c("Object", "Cluster", "Frequency", "Max_frequency", "Min_frequency", "Flag")
                    to_add[, "Object"] = "hh"
                    to_add[, "Cluster"] = cluster
                    to_add[, "Frequency"] = mean(To_plot[index, "Frequency"])
                    to_add[, "Max_frequency"] = max(To_plot[index, "Frequency"])
                    to_add[, "Min_frequency"] = min(To_plot[index, "Frequency"])
                    if (sum(To_plot[index, "Flag"]) > 0) {to_add[, "Flag"] = 1}
                    To_plot_simple = rbind(To_plot_simple, to_add)
                    
                # Same for dpp
                    index  = which(To_plot$Cluster == cluster & To_plot$Object %in% c("dpp_Lib1_dataset", "dpp_Lib2_dataset"))
                    to_add = matrix(data = 0, nrow = 1, ncol = ncol(To_plot_simple))
                    colnames(to_add) = c("Object", "Cluster", "Frequency", "Max_frequency", "Min_frequency", "Flag")
                    to_add[, "Object"] = "dpp"
                    to_add[, "Cluster"] = cluster
                    to_add[, "Frequency"] = mean(To_plot[index, "Frequency"])
                    to_add[, "Max_frequency"] = max(To_plot[index, "Frequency"])
                    to_add[, "Min_frequency"] = min(To_plot[index, "Frequency"])
                    if (sum(To_plot[index, "Flag"]) > 0) {to_add[, "Flag"] = 1}
                    To_plot_simple = rbind(To_plot_simple, to_add)            
            }
            
    # Prepares plotting
        To_plot_simple = as.data.frame(To_plot_simple)
        To_plot_simple$Frequency = as.numeric(To_plot_simple$Frequency)
        To_plot_simple$Max_frequency = as.numeric(To_plot_simple$Max_frequency)
        To_plot_simple$Min_frequency = as.numeric(To_plot_simple$Min_frequency)
        To_plot_simple$Cluster = factor(To_plot_simple$Cluster, levels = unique(To_plot_simple$Cluster))
        To_plot_simple$Object = factor(x = To_plot_simple$Object, levels = c("Reference_dataset", "pxb", "Optix", "vOptix", "dOptix", "dpp", "hh"))
        To_plot_simple$Flag = as.numeric(To_plot_simple$Flag)
    # Saves the table to plot
        dir.create(path = "Results/Class_frequencies_table/", showWarnings = TRUE, recursive = FALSE)
        write.csv(x = To_plot_simple, file = paste("Results/Class_frequencies_table/Class_frequencies_table", Reference_Stage, "reference", paste(Norm_annotations, collapse = "_"), "normalization_simple.csv", sep = "_"))

    # Plots all libraries
        pdf(file = paste("PLOTS/Class_frequencies", Reference_Stage, "reference", paste(Norm_annotations, collapse = "_"), "normalization_simple.pdf", sep = "_"), width = 15, height = 6)
        ggplot(To_plot_simple, aes(x=Cluster, y=Frequency, fill=Object)) +
          geom_bar(position = position_dodge(width = 0.8), width = 0.8, stat="identity") +
          geom_errorbar(aes(ymin=Min_frequency, ymax=Max_frequency), width=0, size = 0.4, position=position_dodge(0.8)) +
          geom_text(aes(label = ifelse(Flag, "*", "")), position = position_dodge(width = 0.8), vjust = -.2, size = 7) +
          scale_fill_manual(values=c("grey", "#E69F00", "#0072B2", "#009E73", "#56B4E9", "#F0E442", "#CC79A7")) +
          #scale_fill_manual(values=c("grey", "#D55E00", "#E69F00", "#CC79A7", "#0072B2", "#56B4E9")) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
                panel.background = element_rect(fill = "white", colour = "grey", size = 1, linetype = "solid"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                #axis.ticks.y = element_blank(),
                legend.key = element_rect(fill = "white"))
        dev.off()
        
    # Plots all libraries except dpp
        To_plot_simple_no_dpp = To_plot_simple
        To_plot_simple_no_dpp = To_plot_simple_no_dpp[To_plot_simple_no_dpp$Object %in% c("Reference_dataset", "pxb", "Optix", "vOptix", "dOptix", "hh"), ]
        pdf(file = paste("PLOTS/Class_frequencies", Reference_Stage, "reference", paste(Norm_annotations, collapse = "_"), "normalization_simple_no_dpp.pdf", sep = "_"), width = 15, height = 6)
        ggplot(To_plot_simple_no_dpp, aes(x=Cluster, y=Frequency, fill=Object)) +
          geom_bar(position = position_dodge(width = 0.8), width = 0.8, stat="identity") +
          geom_errorbar(aes(ymin=Min_frequency, ymax=Max_frequency), width=0, size = 0.4, position=position_dodge(0.8)) +
          geom_text(aes(label = ifelse(Flag, "*", "")), position = position_dodge(width = 0.8), vjust = 0.2, size = 7) +
          scale_fill_manual(values=c("grey", "#0072B2", "#56B4E9", "#CC79A7", "dimgrey", "#D55E00")) +
          #scale_fill_manual(values=c("grey", "#D55E00", "#E69F00", "#CC79A7", "#0072B2", "#56B4E9")) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
                panel.background = element_rect(fill = "white", colour = "grey", size = 1, linetype = "solid"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                #axis.ticks.y = element_blank(),
                legend.key = element_rect(fill = "white"))
        dev.off()
        # same with a black background
        pdf(file = paste("PLOTS/Class_frequencies", Reference_Stage, "reference", paste(Norm_annotations, collapse = "_"), "normalization_simple_no_dpp_black.pdf", sep = "_"), width = 15, height = 6)
        ggplot(To_plot_simple_no_dpp, aes(x=Cluster, y=Frequency, fill=Object)) +
          geom_bar(position = position_dodge(width = 0.8), width = 0.8, stat="identity") +
          geom_errorbar(aes(ymin=Min_frequency, ymax=Max_frequency), width=0, size = 0.4, position=position_dodge(0.8), color = "white") +
          geom_text(aes(label = ifelse(Flag, "*", "")), position = position_dodge(width = 0.8), vjust = 0.2, size = 7, color = "white") +
          scale_fill_manual(values=c("grey", "#0072B2", "#56B4E9", "#CC79A7", "dimgrey", "#D55E00")) +
          #scale_fill_manual(values=c("grey", "#D55E00", "#E69F00", "#CC79A7", "#0072B2", "#56B4E9")) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
                panel.background = element_rect(fill = "black", colour = "grey", size = 1, linetype = "solid"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                #axis.ticks.y = element_blank(),
                legend.key = element_rect(fill = "white"))
        dev.off()
        
           
    # Plots dpp libraries
        To_plot_simple_dpp = To_plot_simple
        To_plot_simple_dpp = To_plot_simple_dpp[To_plot_simple_dpp$Object %in% c("Reference_dataset", "dpp"), ]
        pdf(file = paste("PLOTS/Class_frequencies", Reference_Stage, "reference", paste(Norm_annotations, collapse = "_"), "normalization_simple_dpp.pdf", sep = "_"), width = 10, height = 6)
        ggplot(To_plot_simple_dpp, aes(x=Cluster, y=Frequency, fill=Object)) +
          geom_bar(position = position_dodge(width = 0.8), width = 0.8, stat="identity") +
          geom_errorbar(aes(ymin=Min_frequency, ymax=Max_frequency), width=0, size = 0.4, position=position_dodge(0.8)) +
          geom_text(aes(label = ifelse(Flag, "*", "")), position = position_dodge(width = 0.8), vjust = -.2, size = 7) +
          scale_fill_manual(values=c("grey", "#009E73")) +
          #scale_fill_manual(values=c("grey", "#D55E00", "#E69F00", "#CC79A7", "#0072B2", "#56B4E9")) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
                panel.background = element_rect(fill = "white", colour = "grey", size = 1, linetype = "solid"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                #axis.ticks.y = element_blank(),
                legend.key = element_rect(fill = "white"))
        dev.off()
        # same with black background
        pdf(file = paste("PLOTS/Class_frequencies", Reference_Stage, "reference", paste(Norm_annotations, collapse = "_"), "normalization_simple_dpp_black.pdf", sep = "_"), width = 10, height = 6)
        ggplot(To_plot_simple_dpp, aes(x=Cluster, y=Frequency, fill=Object)) +
          geom_bar(position = position_dodge(width = 0.8), width = 0.8, stat="identity") +
          geom_errorbar(aes(ymin=Min_frequency, ymax=Max_frequency), width=0, size = 0.4, position=position_dodge(0.8), color = "white") +
          geom_text(aes(label = ifelse(Flag, "*", "")), position = position_dodge(width = 0.8), vjust = -.2, size = 7, color = "white") +
          scale_fill_manual(values=c("grey", "#009E73")) +
          #scale_fill_manual(values=c("grey", "#D55E00", "#E69F00", "#CC79A7", "#0072B2", "#56B4E9")) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
                panel.background = element_rect(fill = "black", colour = "grey", size = 1, linetype = "solid"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                #axis.ticks.y = element_blank(),
                legend.key = element_rect(fill = "white"))
        dev.off()
        