# Script: FWL_lipidomics_3.2.2.R
# Author: Wenting, Niklas, Kenny
# Notes:  To start, typing command in the console-----> source("FWL_lipidomics_3.2.2.R") 
#         or press the source butthon. 
#         Please make sure Mac users installed XQuartz and gfortrans.
#         R version 3.6.2 , Rstudio 1.2.1335 for current work.
#         It is based on Niklas and Kenny's previous work (Their work files can be found in folders 
#         quality_control and statistics_quantification). Acknowledge to Weng's technical 
#         guidance, Laura's fatty acid saturation analysis project, Sebstian's shadow experiment.
#         This is the main script for running lipid data analysis.
#         
# Usage:  source("FWL_lipidomics_3.2.2.R") 
#      
# Versions:  LipidSearch 4.1.        
#####################################################################################
rm(list=ls())  ##### empty the environment first. If you need to keep former variables in workspace, just comment this part
setSessionTimeLimit(cpu = Inf, elapsed = Inf)
# source function from FWL_lipidomics.**.functions.R script
source("FWL_lipidomics_3.2.2.FUNCTIONS.R", echo = TRUE)

# options(warn=-1)

# Check directory existence, if not then make one
# all the plots will stored in plot directory and data in the data directory
dirs <- c("plot", "data", "plot/Quantification/classes/fc", "plot/Quantification/classes/svgs/fc", "plot/QC", "plot/Saturation", "plot/Ether", "plot/AcylLength",
          "plot/Volcano", "plot/Pathway", "data/Volcano","data/QC", "data/Quantification", "data/Saturation", "data/Ether", "data/AcylLength", "data/Pathway")
mkdirs(dirs)

message("Are you using PC or MAC?")
robot <- retype_choice("PC/MAC")
if(robot == "mac"){
  message("\n\nPlease make sure you installed Quartz and gfortran before running pipeline.
          \nDo you want to continue?\n")
  process <- retype_choice("Y/N")
}

message("What type of images you want to save? PNG or PDF?\n\n")
image_option <- retype_choice("PNG/PDF")



# read the file from csv directory
csv_files <- list.files(path = "converted", pattern = "\\.csv$")
csv_list <- c()
for(i in seq_along(csv_files)) csv_list[i] <- addquotes(!!as.character(i), " ", !!csv_files[i], "\n")
message("\nThe following files had been generated. \nSelect ONE for subsequent the list of file names:\n", csv_list)
file_option <- readline("Please input the index number of the file: ") %>% as.numeric()
# the file path and name
target_file <- paste("converted/", csv_files[file_option], sep = "")
print(target_file)
lipidomics <- read_csv(target_file, col_types = cols())


##########################################################################################
# QC part I
##########################################################################################
source("FWL_lipidomics_QC_3.2.2.R", echo = FALSE)


##########################################################################################
# Input group information
##########################################################################################
# Make groups
# input the group information and reterieve the data from the csv file.
samples <- Input(filtered_lipidomics2)
# sample info
sample_info <- samples[[1]]





##########################################################################################                                   
# PCA and correlation visualization
##########################################################################################
# pca and correlation plots
label <- "initial"
info_list <- PCA_pairs_Plot(sample_info, filtered_lipidomics2, label, image_option)


# making group repeats according to its position for making groups of later PCA
# sample_raw_list <- info_list[[1]]
# group_repeats <- info_list[[2]]

information <- retrieve_info(sample_info)
group_repeats <- information[[1]]
sample_raw_list <- information[[2]]
# make a data frame contains sample information and group information
group_info <- data.frame(samples=sample_raw_list, 
                         groups=group_repeats, 
                         stringsAsFactors = FALSE) %>% 
  group_by(groups) 
write_csv(group_info, "data/QC/group_information.csv")
group_names <- unique(group_repeats)
ngroups <- length(group_names)

###########################################################################################
# Background subtraction, filter potential invalid lipids
###########################################################################################
message("\nsubtract sample area from background/solvent run?" )
background_option <- retype_choice("Y/N")
filtered_lipidomics <- subtract_not(filtered_lipidomics2, sample_raw_list, background_option, group_info)


###########################################################################################
# Set color theme for plots
###########################################################################################
colors1 <- c("npg", "aaas", "nejm", "jama", "jco", "ucscgb", "d3 ", 
             "locuszoom", "igv", "uchicago", "startrek", "tron", 
             "futurama", "rickandmorty", "simpsons", "gsea", "lancet") 
colors1_no <- c(10, 10, 8, 7, 10, 15, 10, 7, 15, 15, 7, 7, 12, 12, 15, 12, 9) 
colors2 <-  c("BottleRocket1", "BottleRocket2", "Rushmore1", "Royal1", "Royal2", "Zissou1",
              "Darjeeling1", "Darjeeling2", "Chevalier1", "FantasticFox1", "Moonrise1", 
              "Moonrise2", "Moonrise3", "Cavalcanti1", "GrandBudapest1", "GrandBudapest2",
              "IsleofDogs1", "IsleofDogs2") 
colors2_no <- c(7, 5, 5, 4, 5, 5, 5, 5, 4, 5, 4, 4, 5, 5, 4, 4, 6, 5)
# dt1 <- data.table(t(data.frame(colors1, colors1_no)), row.names = c("ggsci_name", "color_number"))
# dt2 <- data.table(t(data.frame(colors2, colors2_no)), row.names = c("wesanderson", "color_number"))
# dt1 <- dt1 %>% select(length(colors1)+1, 1:length(colors1))
# dt2 <- dt2 %>% select(length(colors2)+1, 1:length(colors2))
# colnames(dt1) <- colnames(dt2) <- NULL
colors3 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
message("\nPlease pick a color from ggsci theme link: 
        \nhttps://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html
        \nOR
        \nhttps://github.com/karthik/wesanderson\n\n")
# message("\nPlease pick a color for you plots. 
# \nNotice that the color numbers can't exceed your experiment groups!!!!!
#         \nJust pick the color name from the two tables below.
#         \n OR you could pick 'colors3' theme which has 8 colors. \n\n")
# print(dt1, row.names = TRUE, col.names = "none")
# print(dt2, row.names = TRUE, col.names = "none", topn = 4)
# 

message("\nPlease pick a color for you plots. 
\nNotice that the color numbers can't exceed your experiment groups!!!!!
        \nJust pick the color name from the table in the Viewer panel\n\n")
x1 <- data.frame("color theme" = colors1, "color numbers" = colors1_no)
x2 <- data.frame("color theme" = colors2, "color numbers" = colors2_no)
x3 <- data.frame("color theme" = "colors3", "color numbers" = 8)
xx <- rbind(x1, x2, x3)
yy <- cbind(xx[1:18, ], xx[19:36,])
colnames(yy) <- c("color_theme[1]", "color_number[1]", "color_theme[2]", "color_number[2]")
yy %>% formattable(.)


color_choice <- set_color(colors1, colors2, colors3)
message("Please click plot panel for plot display.")

# visualize AUC of each sample in different lipid classes.
all_samples <- filtered_lipidomics %>%
  as.data.frame() %>%
  select(Class, contains("MainArea"), -"MainArea[c]") %>%
  group_by(Class)%>%
  summarise_at(vars(all_of(sample_raw_list)), list(~sum(., na.rm = TRUE))) %>%
  gather(SAMPLES, all_AUC, -Class) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(GROUPS = ifelse(SAMPLES %in% group_info$samples,
                         unlist(group_info[group_info$samples==SAMPLES, 2]), "NA")) %>%
  mutate(SAMPLES = str_remove_all(SAMPLES, "MainArea\\[") %>% str_remove_all(., "\\]")) %>%
  gather(type, value, -c("Class", "SAMPLES", "GROUPS"))

ordered_samples<- unique(all_samples$SAMPLES) %>% 
  str_remove_all(., "s") %>% 
  as.numeric() %>% 
  sort() %>% 
  paste0("s", .)
all_samples$SAMPLES <- factor(all_samples$SAMPLES, levels = rev(ordered_samples))

params <- c("SAMPLES", "value", "type")
p1 <- plot_all(data = all_samples, params) +
  geom_bar(stat = "identity") +
  facet_wrap(~Class, scales = "free") +
  theme(axis.text = element_text(angle = 30, size = 6, hjust = 1),
        axis.line = element_line(size = 0.2),
        axis.ticks.length = unit(0.6, "mm"),
        axis.ticks = element_line(size = 1)
        ) +
  scale_y_continuous(labels = scientific_format()) +
  coord_flip() + 
  labs(x = "experiment samples", y = "AUC", title = "aggregated AUC for each sample", fill = "") 
print(p1)
ggsave(paste0("plot/QC/raw_all_samples.", image_option), device = image_option, width=20, height = 20)

# # filter negative value
# filtered_samples <- all_samples %>% mutate(value=ifelse(value<=0, 0, value))
# p2 <- plot_all(data = filtered_samples, params) +
#   geom_bar(stat = "identity", position = "fill") +
#   facet_wrap(~Class, scales = "free") +
#   theme(axis.text.x = element_text(angle = 45, size = 8, hjust = 1),
#         axis.line = element_line(size = 0.2)) +
#   scale_y_continuous( expand = c(0, 0, 0.1, 0), labels = scales::percent_format()) +
#   labs(x = "experiment samples", y = "AUC", title = "AUC for each sample", fill = "") 
# print(p2)
# ggsave(paste0("plot/QC/raw_all_samples_percentage.", image_option), device = image_option, width = 20, height = 20)

# check if deleting samples needed and plot new PCA
message("\nPlots can visualized under 'plot' directory or r studio plots panel.\nDo you want to edit sample information for subsequent analyses?")
pca_check <- retype_choice("Y/N")
if(pca_check == "y"){
  info_list <- PCAcheck(pca_check, filtered_lipidomics2, image_option)
  # making group repeats according to its position for making groups of later PCA
  sample_raw_list <- info_list[[1]]
  group_repeats <- info_list[[2]]
  # make a data frame contains sample information and group information
  group_info <- data.frame(samples=sample_raw_list, 
                           groups=group_repeats, 
                           stringsAsFactors = FALSE) %>% 
    group_by(groups) 
  write_csv(group_info, "data/QC/group_information.csv")
  group_names <- unique(group_repeats)
  ngroups <- length(group_names)
  # background subtraction 
  filtered_lipidomics <- subtract_not(filtered_lipidomics2, sample_raw_list, background_option, group_info)
}

# delete unchoosed samples from sample list
total_samples <- filtered_lipidomics %>% select(contains("MainArea[s")) %>% colnames()
deleted_samples <- total_samples %>% subset(!total_samples %in% sample_raw_list) 
filtered_lipidomics <- filtered_lipidomics %>% select(-all_of(deleted_samples))
deleted <- deleted_samples %>% str_remove_all(., "MainArea\\[") %>% str_remove_all(., "\\]") %>% paste0(., sep = ",", collapse = "") %>% remove(.,  pos = -1)
if(length(deleted)>0){
  message("Samples ", deleted, " will be removed from analysis")
}
write_csv(filtered_lipidomics, "data/QC/filtered_lipidomics.csv")



###########################################################################################
# Quantification of lipid class and individual lipid molecules
###########################################################################################
source("FWL_lipidomics_QUANTIFICATION_3.2.2.R", echo = FALSE)


##########################################################################################
# saturation analysis
##########################################################################################
# calculate the saturation for different lipid class
message("The filtered data will be used for saturation analysis.\n")
source("FWL_FattyAcids_Saturation_3.3.3.R", echo = FALSE)
message("The SFA, MUFA, PUFA information will be stored in the count_lipid.csv and aggregated.csv\n")


##########################################################################################
# fatty acids length analysis
##########################################################################################
source("FWL_FattyAcids_Length_3.2.2.R", echo = FALSE)


##########################################################################################
# ether lipid analysis
##########################################################################################
message("\nThe filtered data will be used for ether lipid analysis.\n")
source("FWL_EtherLipids_3.2.2.R", echo = FALSE)


##########################################################################################
# Visualize random sample distribution
##########################################################################################
# visualize sample density plot
# log transformation and visualize approximately normal distribution of the transformed data
log2_lipids <- filtered_lipidomics %>% mutate_at(sample_raw_list, log2trans)
# log2_lipids %>% filter_at(sample_raw_list, any_vars(.<= 0)) 
# randomly choose a sample to see its distribution approximately normal distribution
i <- sample(length(unique(sample_raw_list)), 1)
message(paste("MainArea[s", i, "]", sep=""), " is chosen for plotting distribution")
# the data is approximately normal distribution
plot_all(log2_lipids, sample_raw_list[i]) +
  geom_density() +
  labs(x = "sample distribution (log transformed AUC value)", y = "lipid molecule count") 


##########################################################################################
# log transformation of subtracted/filtered data
##########################################################################################
write.csv(log2_lipids, "data/Volcano/log.molec.csv")


##########################################################################################
# volcano plot
##########################################################################################
# Create a design matrix 
samples <- factor(group_repeats)
design <- model.matrix(~0+samples)
colnames(design) <- levels(samples)
ncomparisons <- readline("How many volcano plots to generate: \n") %>% as.numeric()
impute_dt <- impute_not("y", filtered_lipidomics, sample_raw_list)
imputated_lipids <- impute_dt[[2]]
imputated_lipids <- imputated_lipids %>% ungroup() %>% as.data.frame(., stringsAsFactors=TRUE)
rownames(imputated_lipids) <- imputated_lipids %>% select(LipidMolec) %>% unlist()
processed_lipids <- imputated_lipids %>% select(-Class, -LipidMolec)

# Fit model and extract contrasts 
fit <- lmFit(processed_lipids, design)
# plot volcano graph
VolPlot(ncomparisons, fit, processed_lipids)

message("The pipeline for experiment reference is done here!")


