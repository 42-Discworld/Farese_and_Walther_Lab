#rm(list=ls())  ##### empty the environment first. If you need to keep former variables in workspace, just comment this part

library("FactoMineR")
library("tidyverse")
library("ggplot2")
library("magrittr")
library("ggrepel")
library("reshape")

# Read in file, and store the data into tible structure for the tidyverse to apply
lipidomics <- read_csv(file.choose())
# reformat the colnames for better dealing
colnames(lipidomics) <- make.names(colnames(lipidomics))  

###### if the file has extracted and unextracted samples
### either -contains(colnumber of these samples) or -contains(name of the samples) // contains(pattern)
# eg. 
extracted_P_No <- "APValue.s17"
unextracted_P_No <- "APValue.s18"

# select and count the sample grades you want, filter the p value which less than 0.001
lipidCount <- lipidomics %>% 
  select(contains("Grade"), -contains("Grade.c" ), -contains("extracted"), 
         contains("APValue"),-contains(extracted_P_No), -contains(unextracted_P_No))  %>% 
  transmute(A=rowSums(.=="A"), B=rowSums(.=="B"), 
            C=rowSums(.=="C"), D=rowSums(.=="D"), 
            No_grade=rowSums(.=="-"), APvalue_001=rowSums(. <= "0.001"))
  
# add this count into original strucutes
lipidSelect <- lipidomics %>% bind_cols(lipidCount)

# Filter the dataset based on your criteria  
# flexible parameter 4 which depends on the experiment for total number of A and B 
filtered_lipidomics <- lipidSelect %>% 
                        rowwise() %>% 
                        filter( Rej == 0 &
                                sum(A, B) >= 4 &
                                APvalue_001 >= 4)

# Check how many lipids per class that passed the filtering
summary(as.data.frame(filtered_lipidomics$Class))

#aa <- filtered_lipidomics %>% rowwise %>% filter(Class=="PS")

# QC PLOT 1 - Abundance vs. retention time [requires ggplot2]
ggplot(data=filtered_lipidomics, aes(x = log10(MainArea.s1._DMSO), y = BaseRt)) +
  geom_point() +
  theme_bw() +
  facet_grid(.~Class)


# QC plot for samples difference




# Preparation for pair-wise correlations
inf2NA <- function(x) { x[is.infinite(x)] <- NA; x }


panel.cor <- function(x, y, digits = 2, cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  # correlation coefficient
  r <- cor(inf2NA(x), inf2NA(y), use = "pairwise.complete.obs", method = "pearson")
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste("r= ", txt, sep = "")
  text(0.5, 0.6, txt)
  
  # p-value calculation
  p <- cor.test(x, y)$p.value
  txt2 <- format(c(p, 0.123456789), digits = digits)[1]
  txt2 <- paste("p= ", txt2, sep = "")
  if(is.na(p<0.01)) txt2 <- paste("p= ", "<0.01", sep = "")
  text(0.5, 0.4, txt2)
}

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

# QC PLOT 2 - Pair-wise correlation between replicates
pairs(log10(filtered_lipidomics[, 49:52]), lower.panel=panel.smooth,diag.panel=panel.hist, upper.panel = panel.cor, xlim=c(4,12), ylim=c(4,12))
pairs(log10(filtered_lipidomics[, 53:56]), lower.panel=panel.smooth,diag.panel=panel.hist, upper.panel = panel.cor, xlim=c(4,12), ylim=c(4,12))
pairs(log10(filtered_lipidomics[, 57:60]), lower.panel=panel.smooth,diag.panel=panel.hist, upper.panel = panel.cor, xlim=c(4,12), ylim=c(4,12))
pairs(log10(filtered_lipidomics[, 61:64]), lower.panel=panel.smooth,diag.panel=panel.hist, upper.panel = panel.cor, xlim=c(4,12), ylim=c(4,12))






# Formatting the table for PCA
filtered_lipidomics_PCA <- filtered_lipidomics %>% 
                           select(contains("MainArea.s"), -contains("extract")) %>% t %>% as.data.frame()

colnames(filtered_lipidomics_PCA) <- filtered_lipidomics$LipidMolec                           

log2_filtered_lipidomics_PCA <- log((filtered_lipidomics_PCA+1), 2)

filtered_lipidomics_PCA$Group <- c("DMSO", "DMSO", "DMSO", "DMSO", "DGAT2i", "DGAT2i", "DGAT2i", "DGAT2i", "GPAT3_4i", "GPAT3_4i", "GPAT3_4i", "GPAT3_4i", "DGAT2_GPAT3_4i", "DGAT2_GPAT3_4i", "DGAT2_GPAT3_4i", "DGAT2_GPAT3_4i")

log2_filtered_lipidomics_PCA$Group <- c("DMSO", "DMSO", "DMSO", "DMSO", "DGAT2i", "DGAT2i", "DGAT2i", "DGAT2i", "GPAT3_4i", "GPAT3_4i", "GPAT3_4i", "GPAT3_4i", "DGAT2_GPAT3_4i", "DGAT2_GPAT3_4i", "DGAT2_GPAT3_4i", "DGAT2_GPAT3_4i")
## Alternatively merge with a data frame containing the decode: filtered_lipidomics_PCA <- cbind(filtered_lipidomics_PCA, decode)


# Perform PCA [requires FactoMineR]
res.pca = PCA(filtered_lipidomics_PCA, scale.unit=TRUE, ncp=5, quali.sup=ncol(filtered_lipidomics_PCA), graph=T)

res.pca <- PCA(filtered_lipidomics_PCA[, -ncol(filtered_lipidomics_PCA)], scale.unit=TRUE, ncp=5, quali.sup=ncol(filtered_lipidomics_PCA), graph=T)

concat = cbind.data.frame(filtered_lipidomics_PCA[, ncol(filtered_lipidomics_PCA)], res.pca$ind$coord)

ellipse.coord = coord.ellipse(concat,bary=T)

plot.PCA(res.pca,habillage=ncol(filtered_lipidomics_PCA), ellipse=ellipse.coord, cex=0.8, label="all")

plot.PCA(res.pca,habillage=ncol(filtered_lipidomics_PCA), ellipse=ellipse.coord, cex=0.8, label="quali", legend = "false")


# Output the filtered lipidome
write.csv(filtered_lipidomics, "filtered_lipidomics.csv")


##### STATISTICS WITH LIMMA FOR INDIVIDUAL LIPIDS #####

# Create a design matrix 
samples <- factor(rep(c("wt_fasted", "ko_fasted", "wt_refed", "ko_refed"), rep(4,4)))

design <- model.matrix(~0+samples)

colnames(design) <- levels(samples)

# Reformat the filtered_lipidomics_PCA, consolidate duplicated lipids and log2+1 transform [requires magrittr, dplyr]
filtered_lipids <- filtered_lipidomics_PCA[, -ncol(filtered_lipidomics_PCA)] %>% t %>% as.data.frame()

row.names(filtered_lipids) <- colnames(subset(filtered_lipidomics_PCA, select = -Group))

log2_filtered_lipids <- log((filtered_lipids+1), 2)

write.csv(filtered_lipids, "filtered_lipids.csv")

# Fit model and extract contrasts [requires limma]
fit <- lmFit(log2_filtered_lipids, design)

cont.matrix <- makeContrasts(
  refed_vs_fasted_in_WT=wt_refed-wt_fasted,
  refed_vs_fasted_in_KO=ko_refed-ko_fasted,
  ko_vs_wt_in_fasted=ko_fasted-wt_fasted,
  ko_vs_wt_in_refed=ko_refed-wt_refed,
  Diff=(ko_refed-ko_fasted)-(wt_refed-wt_fasted),levels=design)

fit2 <- contrasts.fit(fit, cont.matrix)

fit2 <- eBayes(fit2)

# Use topTable to obtain significantly altered lipids
Diff <- topTable(fit2, coef="Diff", adjust.method = "fdr",
                 lfc=0, number=nrow(filtered_lipids))

refed_vs_fasted_in_WT <- topTable(fit2, coef="refed_vs_fasted_in_WT", adjust.method = "fdr",
                                  lfc=0, number=nrow(filtered_lipids))

refed_vs_fasted_in_KO <- topTable(fit2, coef="refed_vs_fasted_in_KO", adjust.method = "fdr",
                                  lfc=0, number=nrow(filtered_lipids))

ko_vs_wt_in_fasted <- topTable(fit2, coef="ko_vs_wt_in_fasted", adjust.method = "fdr",
                               lfc=0, number=nrow(filtered_lipids))

ko_vs_wt_in_refed <- topTable(fit2, coef="ko_vs_wt_in_refed", adjust.method = "fdr",
                              lfc=0, number=nrow(filtered_lipids))
# Write to .csv
write.csv(Diff, "Diff.csv")

write.csv(refed_vs_fasted_in_WT, "refed_vs_fasted_in_WT.csv")

write.csv(refed_vs_fasted_in_KO, "refed_vs_fasted_in_KO.csv")

write.csv(ko_vs_wt_in_fasted, "ko_vs_wt_in_fasted.csv")

write.csv(ko_vs_wt_in_refed, "ko_vs_wt_in_refed.csv")

# Vulcano plots with ggplot2 [requires ggplot2 and ggrepel]
refed_vs_fasted_in_WT$Significant <- ifelse(refed_vs_fasted_in_WT$adj.P.Val < 0.05, "Significant", "Not significant")

refed_vs_fasted_in_WT$LipidID <- row.names(refed_vs_fasted_in_WT)

refed_vs_fasted_in_WT$ClassID <- row.names(refed_vs_fasted_in_WT)

refed_vs_fasted_in_WT$ClassID <- sub("\\(.*", "", refed_vs_fasted_in_WT$ClassID)

ggplot(refed_vs_fasted_in_WT, aes(x = logFC , y = -log10(adj.P.Val))) +
  geom_point(aes(color = Significant, alpha = 0.1)) +
  facet_wrap(~ClassID) +
  theme_bw()

ggplot(refed_vs_fasted_in_WT, aes(x = logFC , y = -log10(adj.P.Val))) +
  geom_point(aes(color = Significant)) +
  facet_grid(~ClassID, scales = "free") +
  #xlim(-4, 4) +
  #ylim(-1,15) +
  #xlab("WT 25 mM / WT 0 mM (log2 fold change)") +
  ylab("adj. p-value (-log10)") +
  scale_color_manual(values = c("grey", "red")) +
  theme_bw(base_size = 12) + theme(legend.position = "false") +
  geom_text_repel(
    data = subset(refed_vs_fasted_in_WT, refed_vs_fasted_in_WT$Significant == "Significant" & abs(Diff$logFC) > 0.5),
    aes(label = LipidID),
    size = 3, force = 1) 


##### CLASS LEVEL SUMMARY, NO STATISTICS #####

# Remove standards and sum on a class level
standards <- c("TG(17:1/17:1/17:1)","PS(17:0/20:4)","PI(17:0/20:4)","PG(17:0/14:1)","PE(17:0/14:1)","LPS(17:1)","LPE(17:1)","LPE(17:1)","LPC(17:1)","DG(19:0/19:0)","CerG1(d18:1/12:0)","Cer(d18:1/17:0)","CL(14:1/14:1/15:1/14:1)")

filtered_lipids$Class <- row.names(filtered_lipids)

filtered_lipids <- filtered_lipids[!rownames(filtered_lipids) %in% standards, ]

row.names(filtered_lipids) <- filtered_lipids$Class

filtered_lipids$Class <- sub("\\(.*", "", filtered_lipids$Class)

class_sum <- filtered_lipids %>%
              group_by(Class) %>%
              summarise_all(funs(sum))

write.csv(class_sum, "class_sum.csv")

write.csv(filtered_lipids, "filtered_lipids_wo_standards.csv")

# Alternative 1 - Reshape and plot boxplots with ggplot [requires reshape2, ggplot2]
melt_class_sum <- melt(class_sum)

melt_class_sum$samplegroup <- melt_class_sum$variable

melt_class_sum$samplegroup <- sub(".*\\.", "", melt_class_sum$samplegroup)

melt_class_sum$samplegroup <- sub("_", "", melt_class_sum$samplegroup)

melt_class_sum$samplegroup <- factor(melt_class_sum$samplegroup, levels = c("WT_fasted", "WT_refed", "KO_fasted", "KO_refed"))

melt_class_sum$genotype <- melt_class_sum$samplegroup

melt_class_sum$genotype <- sub("\\_.*", "", melt_class_sum$genotype)

melt_class_sum$genotype <- factor(melt_class_sum$genotype, levels = c("WT", "KO"))

melt_class_sum$diet <- melt_class_sum$samplegroup

melt_class_sum$diet <- sub(".*_", "", melt_class_sum$diet)

ggplot(melt_class_sum) +
  geom_boxplot(aes(x = genotype, y = value, colour = diet)) + 
  facet_wrap(~Class, scales = "free") +
  expand_limits(y = 0) +
  scale_color_manual(values=c("black", "red")) +
  theme_bw()

# Alternative 2 - Create a heatmap with pheatmap [requires pheatmap, RcolorBrewer]
row.names(class_sum) <- class_sum$Class

class_sum$Class <- NULL

pheatmap(class_sum,
         scale = "row", 
         #clustering_method = "complete", 
         #clustering_distance_rows = "correlation",
         #fontsize_row = 4, fontsize_col = 4, 
         #cutree_rows = 7, cutree_cols = 4,
         #cluster_cols = F,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(15))



