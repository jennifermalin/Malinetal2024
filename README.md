# What's in this repository?

## Domain-specific scRNA-seq processing

> Please refer to Simon et al. (In preparation) for the raw data

### Script1_PreparingTheDatasets.R

Prepare the dataset for quantification.

### Script2_AnalysesAndPlots.R

Quantification and visualization used in the manuscript.

## Neuroepithelium quantification

### The xls files starting with 20231115

Imaris quantification for segmented neuroepithelial cells (with Grh).
Channel 2 is staining for Vsx1 and Rx (both from the same species), while
Channel 3 is staining for Optix.


### process.R

The main script for loading, annotating, hypothesis testing, and visualization.

The script produces a raw quantification of the data used in Malin et al.
in quant.csv.

### quant.csv

Please see process.R.


