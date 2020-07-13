list(rm = ls())
source("FWL_lipidomics_3.2.2.FUNCTIONS.R")

# make directories in current folder
dirs <- c("data/Pathway/input", "data/Pathway/output")
mkdirs(dirs)
# put all needed data in the fold pathway/input
message("Please put the files listed in the data/Pathway/input folder")
# or your own pathway data 
message("For proteomics data, e.g. volcano_XXvsXX.txt, lipidpathway.xlsx")
message("\n\nFor lipidomics data generated from lipidomics pipeline:
data/Quantification/aggregated_class.csv, 
data/Saturation/mean_median_filtered_lipidomics.csv,
data/QC/group_information.csv\n")

enzyme_file <- list.files(path = "data/Pathway/input/", pattern = "LipidMetPathways")
proteomics <- read_xlsx(path = paste0("data/Pathway/input/", enzyme_file), sheet = 1) #%>% filter(`Enzyme Type` != "?")
# proteomics <- read_xlsx(file.choose(), sheet = 1) #%>% filter(`Enzyme Type` != "?")

hs_var <- proteomics$`Uniprot ID- Hs` %>% str_split(., ", ") %>% map_dbl(~length(.)) %>% max(.)
mm_var <- proteomics$`Uniprot ID- Mm` %>% str_split(., ", ") %>% map_dbl(~length(.)) %>% max(.)
dm_var <- proteomics$`Uniprot ID- Dm` %>% str_split(., ", ") %>% map_dbl(~length(.)) %>% max(.)
hs_uniprot_ids <- paste0("hs", 1:hs_var, "_uniprot_id")
hs_proteins <- paste0("hs", 1:hs_var,"_protein_name")
hs_genes <- paste0("hs", 1:hs_var, "_gene_name")

mm_uniprot_ids <- paste0("mm", 1:mm_var, "_uniprot_id")
mm_proteins <- paste0("mm", 1:mm_var, "_protein_name")
mm_genes <- paste0("mm", 1:mm_var, "_gene_name")

dm_uniprot_ids <- paste0("dm", 1:dm_var, "_uniprot_id")
dm_proteins <- paste0("dm", 1:dm_var, "_protein_name")
dm_genes <- paste0("dm", 1:dm_var, "_gene_name")

# get data for human sample
hs <- map2(str_split(proteomics$`Uniprot ID- Hs`, ",\\s*"), str_split(proteomics$`Protein name(s)- Hs`, ","), ~paste(.x, .y, sep = "_")) 
hss <- map2(hs, str_split(proteomics$`Gene name(s)- Hs`, ",\\s*"), ~paste(.x, .y, sep = "_"))
hs_list <- lapply(hss, function(x)as.list(x)) 
hss_list <- lapply(hs_list, function(x) paste(x, collapse = ", "))
hs_dt <-  data.frame(hs = matrix(unlist(hss_list), nrow=length(hss_list), byrow=T))
dt1 <- proteomics %>% select(`Arrow ID`, 'Substrate', 'Product', 'Enzyme Type') %>% bind_cols(., hs_dt) %>% gather(groups, type, -c(`Arrow ID`, 'Substrate', 'Product', 'Enzyme Type'))
options(warn = -1)
p1 <- dt1 %>% 
  separate(., type, into = hs_uniprot_ids, sep = ", ") %>% 
  gather(id, hs, -c("Arrow ID", 'Substrate', 'Product', 'Enzyme Type', "groups")) %>% 
  select(-id) %>% 
  filter(!str_detect(hs, "NA")) %>% 
  separate(., hs, into = c("Uniprot ID", "Protein name", "Gene name"), sep = "_\\s*")
options(warn = 0)

# get data for mouse
mm <- map2(str_split(proteomics$`Uniprot ID- Mm`, ",\\s*"), str_split(proteomics$`Protein name(s)- Mm`, ","), ~paste(.x, .y, sep = "_")) 
mms <- map2(mm, str_split(proteomics$`Gene name(s)- Mm`, ",\\s*"), ~paste(.x, .y, sep = "_"))
mm_list <- lapply(mms, function(x)as.list(x)) 
mms_list <- lapply(mm_list, function(x) paste(x, collapse = ", "))
mm_dt <-  data.frame(mm = matrix(unlist(mms_list), nrow=length(mms_list), byrow=T))

dt2 <- proteomics %>% select(`Arrow ID`, 'Substrate', 'Product', 'Enzyme Type') %>% bind_cols(., mm_dt) %>% gather(groups, type, -c(`Arrow ID`, 'Substrate', 'Product', 'Enzyme Type'))
options(warn = -1)
p2 <- dt2 %>% 
  separate(., type, into = mm_uniprot_ids, sep = ", ") %>% 
  gather(id, mm, -c("Arrow ID", 'Substrate', 'Product', 'Enzyme Type', "groups")) %>% 
  select(-id) %>% 
  filter(!str_detect(mm, "NA")) %>% 
  separate(., mm, into = c("Uniprot ID", "Protein name", "Gene name"), sep = "_\\s*")
options(warn = 0)

# get data for fly
dm <- map2(str_split(proteomics$`Uniprot ID- Dm`, ",\\s*"), str_split(proteomics$`Protein name(s)- Dm`, ","), ~paste(.x, .y, sep = "_")) 
dms <- map2(dm, str_split(proteomics$`Gene name(s)- Dm`, ",\\s*"), ~paste(.x, .y, sep = "_"))
dm_list <- lapply(dms, function(x)as.list(x)) 
dms_list <- lapply(dm_list, function(x) paste(x, collapse = ", "))
dm_dt <-  data.frame(dm = matrix(unlist(dms_list), nrow=length(dms_list), byrow=T))
dt3 <- proteomics %>% select(`Arrow ID`, 'Substrate', 'Product', 'Enzyme Type') %>% bind_cols(., dm_dt) %>% gather(groups, type, -c(`Arrow ID`, 'Substrate', 'Product', 'Enzyme Type'))
options(warn = -1)
p3<- dt3 %>% 
  separate(., type, into = dm_uniprot_ids, sep = ", ") %>% 
  gather(id, dm, -c("Arrow ID", 'Substrate', 'Product', 'Enzyme Type', "groups")) %>% 
  select(-id) %>% 
  filter(!str_detect(dm, "NA")) %>% 
  separate(., dm, into = c("Uniprot ID", "Protein name", "Gene name"), sep = "_\\s*")
options(warn = 0)
hs_mm <- rbind(p1, p2) %>% mutate(`Protein name` = str_to_upper(`Protein name`)) %>% arrange(`Protein name`)
write_csv(hs_mm, "data/Pathway/output/hs_mm.csv")

hs_mm_dm <- rbind(hs_mm, p3) %>% mutate(`Protein name` = str_to_upper(`Protein name`))  %>% arrange(`Protein name`)
write_csv(hs_mm_dm, "data/Pathway/output/hs_mm_dm.csv")
hs_mm_wd <- hs_mm %>% mutate(`Gene name` = str_to_upper(`Gene name`)) %>% spread(., groups, "Uniprot ID")
#hs_mm_dm_wd <- hs_mm_dm %>% mutate(`Gene name` = str_to_upper(`Gene name`)) %>% spread(., groups, "Uniprot ID")
write_csv(hs_mm_wd, "data/Pathway/output/hs_mm.wide.csv")
#write_csv(hs_mm_dm_wd, "data/Pathway/output/hs_mm_dm.wide.csv")

# read lipidomics data
data_class <- read_csv("data/Pathway/input/aggregated_class.csv", col_types = cols()) %>% arrange(Class)
# saturation data: data/Saturation/mean_median_filtered_lipidomcis.csv
saturation_data <- read_csv("data/Pathway/input/mean_median_filtered_lipidomics.csv", col_types = cols())
# load group information: data/QC/group_information.csv
group_information <- read_csv("data/Pathway/input/group_information.csv", col_types = cols())
group_nm <- group_information$groups %>% unique()
ctr <- check_group(group_nm, "control")
contrast <- check_group(group_nm, "contrast")
group_info <- group_information %>% filter(groups %in% c(ctr, contrast))
sample_raw_list <- group_info  %>% select(samples) %>% unlist()
names(sample_raw_list) <- NULL
group_names <- c(ctr, contrast)
ngroups <- 2

# choose control group value for normalization 
method <- retype_choice("MEAN/MEDIAN")
# select data for fold change calculation
control_name <- paste0(ctr, '_', method)

# calculate the relative fold change
dt <- cal_lipid_statistics(data_class, group_info, method, "Class")
dt_ctr <- dt[[1]] %>% select(Class, all_of(control_name)) %>% left_join(data_class, ., by = "Class")
normalized_dt <- norm_by_mean2(dt_ctr, "Class", sample_raw_list, control_name, "n")
dt_nm <- cal_lipid_statistics(normalized_dt, group_info, method, "Class")

# calculate saturation 
groups <- paste0(group_names, "_", method, "_SFA")
saturation_dt <- saturation_data %>% select(Class, contains(groups))
sdt <- FormatData2(saturation_dt, method) %>% rename(SFA = mean)
dd <- dt[[2]] %>% mutate_if(is.factor, as.character)
dt_st <- full_join(dd, sdt) %>% filter(!is.na(TYPE)) %>% mutate(saturation = SFA/mean)
samples <- paste0(sample_raw_list, "_SFA")
saturation_samples <- saturation_data %>% select(Class, all_of(samples)) %>% arrange(Class)
sample_SFA <- saturation_samples[, -1]/data_class[, -1] 
ds <- saturation_samples %>% select(Class) %>% bind_cols(., sample_SFA) 
colnames(ds) <- c("Class", sample_raw_list)
dts <- cal_lipid_statistics(ds, group_info, method, "Class")
dt_st <- dts[[2]] %>% rename(saturation = mean)

# merge foldchange and saturation information
dtn <- dt_nm[[2]] %>% mutate_if(is.factor, as.character)
dt_all <- dt_st %>% 
  mutate(Groups = as.character(Groups)) %>% 
  select(Class, Groups, saturation) %>% 
  left_join(dtn, .) %>% 
  mutate(# convert lipid class name to standards
    Class = replace(Class, Class == "DG", "DAG"),
    Class = replace(Class, Class == "SPH", "So"))

# store lipid data information under Pathway folder
lipid_name<- readline("Please type the file name for storage, e.g. flox_lko.csv:  ")
write_csv(dt_all, paste0("data/Pathway/output/", lipid_name))

# read proteomics data
#protein <- read_csv("data/Pathway/output/hs_mm.wide.csv", col_types = cols())
protein <- read_csv("data/Pathway/output/hs_mm_dm.csv", col_types = cols())
message("\nAre you using source from homo sapiens (Hs), Mus musculus (Nm) or fly (Dm) for experiment?")
species <- retype_choice("Hs/Mm/Dm")

# get proteomics data 
message("\nPlease open the proteomics data, e.g. volcano_floxVSlko.csv.\n\n")

sample_proteomics <- read.table(file.choose(), sep="\t", stringsAsFactors = FALSE)
sample <- sample_proteomics[-1, ] 
colnames(sample) <- sample_proteomics %>% slice(1)
# reformat
sample_id <- sample %>% 
  #select(Difference, "Prtein IDs", "Majority protein IDs") %>% 
  select(2:5) %>% 
  mutate(id = strsplit(`Protein IDs`, ";\\s*")) #%>% 
  #rename("-log(Pvalue)" = !!"-LOG(P-value)")
colnames(sample_id)[1] <- "-log(Pvalue)"

srows <- plyr::ldply(sample_id$id, rbind) %>% mutate_if(is.factor, as.character)

id <- sample_id %>% 
  select(-id) %>% 
  bind_cols(., srows) %>% 
  gather("type", "ID", -c("Protein IDs", "Difference", "Majority protein IDs", "-log(Pvalue)")) %>% 
  filter(!is.na(`ID`)) %>% 
  select(-type, -"Protein IDs")
write_csv(id, "data/Pathway/output/raw_detected.csv")


selected_info <- protein %>% filter(groups == species)
detected_proteins <- id %>% 
  rowwise() %>% 
  mutate(proteins = ifelse(ID %in% selected_info$`Uniprot ID`, unlist(selected_info[selected_info$`Uniprot ID`==`ID`, "Protein name"]), NA), 
            Arrows = ifelse(ID %in% selected_info$`Uniprot ID`, unlist(selected_info[selected_info$`Uniprot ID`==`ID`, "Arrow ID"]), NA),
            Enzyme = ifelse(ID %in% selected_info$`Uniprot ID`, unlist(selected_info[selected_info$`Uniprot ID`==`ID`, "Enzyme Type"]), NA)) %>% 
  filter(!is.na(proteins)) %>% 
  mutate(Difference = as.numeric(Difference)) %>% 
  ungroup()

detected_proteins
write_csv(detected_proteins, "data/Pathway/output/detected.csv")

# get enzyme list
protein_list <- protein$`Protein name` %>% unique() %>% na.omit()

all_nodes <- c(union(unique(protein$Substrate), unique(protein$Product)), "Acetyl-CoA", "Cholesterol", "CE")

link <- protein %>% select("Arrow ID", Substrate, Product, "Enzyme Type") %>% distinct() %>% rename(arrow_id = "Arrow ID", from = Substrate, to = Product, label = "Enzyme Type")
head <- data.frame(arrow_id = c(NA, NA), from = c("Acetyl-CoA","Acetyl-CoA"),  to = c('FA', 'Cholesterol'), label = c(NA, NA))
link <- bind_rows(head, link) %>% 
  # mutate(label = ifelse(label %in% detected_proteins$Enzyme, label, NA)) %>% 
  rowwise() %>% 
  mutate(text = detected_proteins[which(detected_proteins$Arrows == arrow_id), c("proteins", "Difference")] %>% 
           transmute(tx = paste0(proteins, ": ", round(Difference, 2), "\n\n")) %>% 
           select(tx) %>% 
           unlist() %>% 
           paste(., collapse = ""))

write_csv(link, "data/Pathway/output/edges.csv")

default_nodes <- setdiff(all_nodes, unique(dt_all$Class))
longley_nodes <- setdiff(unique(dt_all$Class), all_nodes)
miss_data <- tibble(Class = rep(default_nodes, each = ngroups),
                    Groups = rep(group_names, length(default_nodes)),
                    #mean = rep(c(1, rep(0, ngroups-1)), length(default_nodes)),
                    mean = rep(1, 2*length(default_nodes)),
                    saturation = rep(NA, ngroups*length(default_nodes)),
                    # color_saturation_code = rep(NA, ngroups*length(default_nodes)),
                    #fc = rep(c(.01, rep(0, ngroups-1)), length(default_nodes))
)
paths <- bind_rows(dt_all, miss_data)
path_data <- set_color_code(paths, "saturation", "color_saturation_code")
edge_data <- set_color_code(detected_proteins, "Difference", "color_edge_code")
path_data$fc <- path_data$mean*0.01
data <- path_data
write_csv(data, "data/Pathway/output/nodes_all.csv")



build_node <- function(info, default_node){
  node <- data.frame(
    id = info$Class, 
    label=info$Class, 
    size = info$mean, 
    scaled_size = info$fc, 
    shape = "dot" , 
    saturation = info$saturation, 
    color_saturation = info$color_saturation_code)  %>% 
    arrange(desc(color_saturation))
  node$color.background = c(viridis_pal()(nrow(node) - length(default_node)), 
                            rep(NA, length(default_node)))
  return(node)
}

info_ctr <- data %>% filter(Groups==ctr) 
node_ctr <- build_node(info_ctr, default_nodes)
info_contrast <- data %>% filter(Groups == contrast)
node_contrast <- build_node(info_contrast, default_nodes)
#node$color.background <- c(viridis_pal()(nrow(node) - length(default_nodes)), rep(NA, length(default_nodes))) # set node color

link$label <- link$text
link$label.cex <- seq(1, 2,length.out = nrow(link))

plot_pathway <- function(node, edge, group){
  node$shadow <- TRUE
  node$size <-  node$value # set fold change of group mean as node size
  node$color.border <- "black" 
  cs <- na.omit(node$color.background)
  lnode <- data.frame(label = c("Highest \nsaturation\n value", "Lowest\n saturation\n value"),
                      color = c(cs[1], cs[length(cs)]),
                      shape = "box",
                      font.color = "white",
                      font.size = 10)
  #view(node)
  network <- visNetwork(node, edge, main = addquotes(!!group, " pathway"), physics = TRUE) %>% 
    #  visOptions(manipulation = TRUE) %>%
    
    visEdges(arrows = "from") %>% 
    visLegend(addNodes = lnode) %>% 
    # visConfigure(enabled = TRUE) %>%
    visIgraphLayout(#layout = "layout_with_fr",
      smooth = list(enabled = TRUE,
                    type = "diagonalCross",
                    roundness = .7)) 
  
  return(network)
}
set.seed(1234)
path_ctr <- plot_pathway(node_ctr, link, ctr)
path_ctr 
visSave(path_ctr, file = paste0(ctr, ".html"))
path_contrast <- plot_pathway(node_contrast, link, contrast)
path_contrast 
visSave(path_contrast, file = paste0(contrast, ".html"))



# build node and edge for gephi
gephi_nodes_ctr <- node_ctr
gephi_nodes_contrast <- node_contrast
build_gephi <- function(gephi_nodes, link, group){
  colnames(gephi_nodes)[1] <- colnames(gephi_nodes)[1] <- "ID"
  gephi_nodes$color <- ifelse(is.na(gephi_nodes$saturation), -1, gephi_nodes$saturation)
  gephi_nodes$groups <- group
  gephi_links <- link
  colnames(gephi_links) <- case_when(colnames(gephi_links) == "from" ~"Source",
                                     colnames(gephi_links) == "to" ~ "Target", 
                                     TRUE~colnames(gephi_links))
  gephi_links <- gephi_links %>% mutate(text = str_replace_all(text, "\n", "|"))
  return(list(gephi_nodes, gephi_links))
}
gephi_info <- build_gephi(node_ctr, link, ctr) 
gephi_node_ctr <- gephi_info[[1]]
gephi_links <- gephi_info[[2]]
gephi_node_contrast <- build_gephi(node_contrast, link, contrast)[[1]] 
message("For Gephi, experiment control and contrast groups node data need be loaded separately with edges data.")
write_csv(gephi_node_ctr, "data/Pathway/output/gephi_node_ctr.csv")
write_csv(gephi_node_contrast, "data/Pathway/output/gephi_node_contrast.csv")
write_csv(gephi_links, "data/Pathway/output/gephi_link.csv")

# build node and edges for cytoscape
gephi_nodes <- bind_rows(gephi_node_ctr, gephi_node_contrast) 
cyto_node <- gephi_nodes %>% 
  select(ID, label,size, color.background, color, groups)  %>% 
  mutate(Source = label)

message("For cytoscape, you need to upload the ct_table.csv for all nodes and edges first.
Then you could upload ct_link.csv. If you want to modify the nodes attributes like size, please upload ct_node.csv")

cyto_link <- gephi_links %>% 
  select(Source, Target, text, label.cex) %>% 
  mutate(text = str_replace_all(text, "\n", " "),)
write_csv(cyto_node, "data/Pathway/output/ct_node.csv")
write_csv(cyto_link, "data/Pathway/output/ct_link.csv")
cyto_dt <- left_join(cyto_link, cyto_node, by = "Source") %>% ungroup()
write_csv(cyto_dt, "data/Pathway/output/ct_table.csv")

# # tiny shiny experiment
# ui <- fluidPage(
#   fluidRow(
#     column(6,
#            visNetworkOutput('vis1')),
#     column(6,
#            visNetworkOutput('vis2'))
#   )
# )
# 
# server <- function(input, output) {
#   output$vis1 <- renderVisNetwork(
#     {
#       plot_pathway(node_ctr, link, ctr) %>%
#         #visNodes(label = vis.node$id)
#         visOptions(manipulation = TRUE)
# 
#     })
#   output$vis2 <- renderVisNetwork(
#     {
#       plot_pathway(node_contrast, link, contrast) %>%
#         #visNodes(label = vis.node$id)
#         visOptions(manipulation = TRUE)
# 
#     })
# }
# 
# shinyApp(ui = ui, server = server)
