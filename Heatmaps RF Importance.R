library(ggplot2)
library(pheatmap)
library(grid)

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(stringr)

setwd("C:/Users/Raphael/OneDrive/Área de Trabalho/(UFGRS) input_MicrobiomeAnalyst(non-contamination)-20240807T215735Z-001/Heatmap-RF/Heatmap 05")

rf_data <- read.csv2("rf_100_health_vs_caries-rankImportance.csv")

# Filtring only significant ASVs
rf_data_filtered <- rf_data %>% filter(Significant == TRUE)

asv_data <- read.delim("asv_RIGHT_NoNa-Genus.txt", header = TRUE, sep = "\t", fill = TRUE, check.names = FALSE, quote = "")

data_metadata <- read.delim("metadata-input_MicrobiomeAnalyst_COMESTUDO.txt", header = TRUE, sep = "\t", fill = TRUE, check.names = FALSE, quote = "")

# Joining two dataframes RandomForest + ASV_data (count_table)
colnames(asv_data)[1] <- "ASV" 
merged_data <- rf_data_filtered %>% inner_join(asv_data, by = "ASV")

# Merging data based in comuns ASVs
merged_data_summed <- merged_data %>%
  group_by(ASV) %>%
  summarise(
    across(where(is.numeric), sum, na.rm = TRUE),  # Soma colunas numéricas
    across(where(is.character), ~ paste(unique(.), collapse = ", ")),  # Concatenar valores únicos para colunas character
    .groups = "drop"  # Remover o agrupamento após a operação
  ) %>%
  select(ASV, freq, Caries, Health, MeanDecreaseAccuracy, MeanDecreaseGini, Genus, Phylum, Rank, RF_type, everything())  # Reorganizar as colunas

# Removing columns that aren't equal 0
merged_data_summed <- merged_data_summed %>%
  select(where(~ !is.numeric(.) || sum(., na.rm = TRUE) != 0))

merged_data_long <- merged_data_summed %>%
  pivot_longer(cols = starts_with("SRR"),  # Select columns that start with "SRR"
               names_to = "#NAME",          
               values_to = "Value")         

# Union data by id = "#NAME"
merged_final <- data_metadata %>%
  inner_join(merged_data_long, by = c("#NAME" = "#NAME"))

merged_final <- merged_final %>% select(-c(Rank, RF_type, MeanDecreaseGini, MeanDecreaseAccuracy, freq))
colnames(merged_final)[13] <- "Count"
merged_final$STUDY <- sub(",.*", "", merged_final$STUDY)


####################################
####### STARTING FUNCTIONS #########
####################################

save_pheatmap_png <- function(x, filename, width = 3200, height = 8000, res = 300) {
  # Adjust the width and height for landscape mode
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

graphs_heatmap <- function(heatmap_matrix, ann_col, ann_colors, study_name){
  
  ### Heatmap (whitout LEGEND and ANNOTATION) ###
  p <- pheatmap(heatmap_matrix, 
           annotation_col = ann_col,  # Use annotation_row instead of annotation_col
           cluster_rows = TRUE, 
           cluster_cols = TRUE, 
           scale = "none",
           main = "",
           color = colorRampPalette(c("white","yellow", "orange","red"))(200),
           annotation_colors = ann_colors,
           fontsize = 6, cellwidth = 7, cellheight = 7, width = 10, height = 15,
           annotation_legend = FALSE,
           legend = FALSE, 
           annotation_names_col = FALSE,
           gaps_col = 200)  # Adjust size for better visibility
  
  mh_name <- "heatmap"
  ext <- ".png"
  save_pheatmap_png(p, paste0("heatmap_", study_name, " - nonLegend_nonAnn", ext))
  
  ### Heatmap (whitout LEGEND = scale) ###
  p <- pheatmap(heatmap_matrix, 
                annotation_col = ann_col,  # Use annotation_row instead of annotation_col
                cluster_rows = TRUE, 
                cluster_cols = TRUE, 
                scale = "none",
                main = "",
                color = colorRampPalette(c("white","yellow", "orange","red"))(200),
                annotation_colors = ann_colors,
                fontsize = 6, cellwidth = 7, cellheight = 7, width = 10, height = 15,
                annotation_legend = TRUE,
                legend = FALSE, 
                annotation_names_col = FALSE,
                gaps_col = 200)  # Adjust size for better visibility
  
  mh_name <- "heatmap"
  ext <- ".png"
  save_pheatmap_png(p, paste0("heatmap_", study_name, " - nonLegend", ext))
  
  ### Heatmap (whitout ANNOTATION = Caries/Health legend condition) ###
  p <- pheatmap(heatmap_matrix, 
                annotation_col = ann_col,  # Use annotation_row instead of annotation_col
                cluster_rows = TRUE, 
                cluster_cols = TRUE, 
                scale = "none",
                main = "",
                color = colorRampPalette(c("white","yellow", "orange","red"))(200),
                annotation_colors = ann_colors,
                fontsize = 6, cellwidth = 7, cellheight = 7, width = 10, height = 15,
                annotation_legend = FALSE,
                legend = TRUE, 
                annotation_names_col = FALSE,
                gaps_col = 200)  # Adjust size for better visibility
  
  mh_name <- "heatmap"
  ext <- ".png"
  save_pheatmap_png(p, paste0("heatmap_", study_name, " - nonAnn", ext))
  
  ### Heatmap - Everything (LEGEND+ANN) ###
  p <- pheatmap(heatmap_matrix, 
                annotation_col = ann_col,  # Use annotation_row instead of annotation_col
                cluster_rows = TRUE, 
                cluster_cols = TRUE, 
                scale = "none",
                main = "",
                color = colorRampPalette(c("white","yellow", "orange","red"))(200),
                annotation_colors = ann_colors,
                fontsize = 6, cellwidth = 7, cellheight = 7, width = 10, height = 15,
                annotation_legend = TRUE,
                legend = TRUE, 
                annotation_names_col = FALSE,
                gaps_col = 200)  # Adjust size for better visibility
  
  mh_name <- "heatmap"
  ext <- ".png"
  save_pheatmap_png(p, paste0("heatmap_", study_name, ext))

}

heatmap_gen <- function(dataframe, study_name){
  workspace <- getwd()
  
  study_name <- dataframe$STUDY[1]
  study_name <- str_sub(study_name, end = -2)
  dir.create(study_name)
  setwd(study_name)
  
  heatmap_data <- dataframe %>%
    select(NAME, ASV, Count, STATUS) %>%
    pivot_wider(names_from = ASV, values_from = Count, values_fill = list(Count = 0))
  
  # Convert to data frame and set row names
  heatmap_data <- as.data.frame(heatmap_data)
  rownames(heatmap_data) <- heatmap_data$NAME
  heatmap_data <- heatmap_data %>% select(-NAME, -STATUS)  # Remove STATUS before conversion
  
  # Convert all values to numeric
  heatmap_matrix <- as.matrix(heatmap_data)  # No need for sapply(), already numeric
  #write.csv2(heatmap_matrix, file = "SumCounts_ASVxSAMPLE_NAME.csv") 
  
  heatmap_matrix <- log1p(heatmap_matrix)
  heatmap_matrix <- heatmap_matrix[, colSums(heatmap_matrix, na.rm = TRUE) != 0]
  write.csv(heatmap_matrix, file = paste0("log1pCounts_ASVxSAMPLE_NAME - ", study_name, ".csv"))
  
  # Extract STATUS for annotation (removing duplicates)
  status_annotation <- merged_final %>%
    distinct(NAME, STATUS) %>%
    column_to_rownames(var = "NAME")
  
  status_annotation <- status_annotation[rownames(heatmap_matrix), , drop = FALSE]
  # Order STATUS in graph (Didn't work!)
  status_annotation$STATUS <- factor(status_annotation$STATUS, levels = c("Health", "Caries"))
  ann_colors = list(STATUS = c(Health="#4ebecd", Caries="#f1908a"))
  heatmap_matrix_transposed <- t(heatmap_matrix)
  
  graphs_heatmap(heatmap_matrix_transposed, status_annotation, ann_colors, study_name)
  
  # Generate heatmap
  message(paste0("Generated heatmap for STUDY: ", study_name))
  setwd(workspace)
}

filter_study <- function(data, study_name) {
  filtered_data <- data %>%
    filter(STUDY == study_name)
  return(filtered_data)
}

##################################
####### ENDING FUNCTIONS #########
##################################


colnames(merged_final)[2] <- "NAME"

merged_final$Count <- as.numeric(merged_final$Count)
#write.csv(merged_final, file = "merged_final_counts_sum.csv" )

studys <- unique(merged_final$STUDY)

dir.create("Study Heatmaps")
setwd("Study Heatmaps")

for (study_name in studys) {
  df_study <- filter_study(merged_final, study_name)
  heatmap_gen(df_study, study_name)
}  

