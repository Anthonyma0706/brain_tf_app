# brain_TF_app
## About  
This shiny app displays transcription factor activity inference data from a developmental timecourse of the mouse forebrain and pons in three main tabs including:

  - Tab 1 (tf and gene datatable and network graph): Given a list of TFs, display a table of those TFs with their inferred target genes, Normalized Enrichment Score (NES) representing the activity level, motifs, etc. Secondly, based the same gene regulatory network, a network graph visualization displaying the TF and targets, where overlaps in targets of different TFs can be visualized easily. Color and size of nodes designate different types of genes. Genes (nodes) which are part of a user-defined list can also be coloured to highlight members of a pathway.

  - Tab 2 (heatmap and UMAP clustering plot): For this tab, we visualize the activity levels of transcription at the cell and cluster level. Firstly, two heatmaps visualize the activity per cluster/cell, with cells annotated by cluster. Secondly, we show the TF activity in the cells in a two-dimensional embedding (UMAP), to highlight regions (clusters/cell types) with high activity.

  - Tab 3 (timeseries ribbon plot), we use the binarized TF activity to plot the proportion of cells where the TF is active at each timepoint during embryonic and postnatal development (E12.5, E15.5, P0, P3, P6). Moreover, we use plotly to display an interactive data visualization – hover over the plot to display the cluster label, and click on the cell type legend to restrict the visualization to selected clusters of interest.  


## Features
- The first tab generates a **cytoscape network graph**, displaying the input TFs, their target genes, the mutual target gene and gene pathway analysis. 
We color the input genes related to input tfs(purple nodes) in green. Further demonstration in terms of how to generate and manipulate the data required in the cytoscape network code is in the app_description.Rmd(html) file in this repository.

## Structure  
### data
We use the mouse brain data with cortex(forebrain) and pons region. Data should have the exactly same format among different regions in this app. 
Tab 1 and Tab 2 use more general data for cancer, so they can also apply to standard cancer samples, while tab3 uses more specific time couse data used in Kleinman lab that will be detached in a more general TF cancer app.
check out <https://github.com/fungenomics/TF_explorer_app> for a more general use and data transformation of tab1 and tab2 functionality using standard cancer sample data format.  
### User input and missing data  
User needs to select the region first, and then select the transcription factors. Same tf list will be used throught out the whole tab. 
However, sometimes you may see some tfs are on display in some tabs only (may not see the information of that transcription factor or the plot is not updated, etc.) This is because some transcription factors from your input may not have the corresponding data in the some tabs. Missing data of those tfs mean that they're not active in the timepoints during the collection of data.  
This app doesn't provide any error message if the tf data is empty.
### data structure
We generate all the data in data_prep.R which saves a list like this for each brain region:
*data_cortex <- list(
  "cell_metadata"  = forebrain_data,
  "TF_and_ext" = TF_and_ext,
  "TF_target_gene_info" = TF_target_gene,
  "unique_active_TFs_bare" = unique_TF,
  "active_TFs" = TF_active,
  "binary_active_TFs" = tf_df,
  "timeseries_input_meta" = timeseries_input_meta_cortex,
  "binary_activity" = binary_activity,
  "tfs_not_exist_timeseries" = l_nexist_cortex
)*
where all the names in the list will be the same for data_cortex and data_pons.
At the beginning of the app, we load those datasets as list, then based on the region selected by user, the corresponding dataset is assigned to the input_new() list in server.R,
which will be updated by the update button, see server.R for more details.















