This script would need data generated from lipidomics pipeline analysis and proteomics analysis.

1. Step 1: Put data under directory data/Pathway/input

  * For lipidomics data generated from lipidomics pipeline:
  
    * data/Quantification/aggregated_class.csv, 
    
    * data/Saturation/mean_median_filtered_lipidomics.csv,
    
    * data/QC/group_information.csv
    
  * For Proteomics data:
  
    *volcano_XXVSXX.txt, LipidMetPathwaysXX.xlsx (Please download the test version on github)
    
1. Step 2: Type control group name and contrast group name individually for pathway visualization in console.

1. Step 3: Type one method from mean or median for calculating relative fold change.

1. Step 4: Type the species from your experiment, e.g. Hs/Mm/Dm: Mm

1. Step 5: Type file name, e.g. flox_lko.csv for storing extacted information.

1. Step 6: Open the proteomics data which you put in data/Pathway/input, e.g. volcano_floxVSlko.txt.

1. Step 7: Check all the output under data/Pathway/output.

  * Please note that pdf or png files could be generated via gephi or cytoscape.
  
  * Gephi would need to plot two group pathway individually, load gephi_node_ctr.csv or and gephi_node_contrast with gephi_link.csv.
  
  * Cytoscape would need to import all node and edge information first via ct_table.csv including experiment groups information. Then load the edge table ct_link.csv. 
  For node size attributes, user need to load ct_node.csv. User need to add filter for filtering expriment group for analysis.
