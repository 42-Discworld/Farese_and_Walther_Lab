This script would need data generated from lipidomics pipeline analysis and proteomics analysis.

Step 1: Put data under directory data/Pathway/input

  * For lipidomics data generated from lipidomics pipeline:
  
    * data/Quantification/aggregated_class.csv, 
    
    * data/Saturation/mean_median_filtered_lipidomics.csv,
    
    * data/QC/group_information.csv
    
  * For Proteomics data:
  
    *volcano_XXVSXX.txt, LipidMetPathwaysXX.xlsx (Please download the test version on github)
    
Step 2: Type control group name and contrast group name individually for pathway visualization in console.

Step 3: Type one method from mean or median for calculating relative fold change.

Step 4: Type the species from your experiment, e.g. Hs/Mm/Dm: Mm

Step 5: Type file name, e.g. flox_lko.csv for storing extacted information.

Step 6: Open the proteomics data which you put in data/Pathway/input, e.g. volcano_floxVSlko.txt.

Step 7: Check all the output under data/Pathway/output.

  * Please note that pdf or png files could be generated via gephi or cytoscape.
  
  * Gephi would need to plot two group pathway individually, load gephi_node_ctr.csv or and gephi_node_contrast with gephi_link.csv.
  
  * Cytoscape_3.8.0 will need to load edges table (ct_link.csv) and node tables (ct_node.csv)
  
  
  
 Cytoscape version 
