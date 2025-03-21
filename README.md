# Heatmap_RandomForest

This code get Random Forest result table from Microbiota Analysis (16S), and generate a Heatmaps for each study (on Metadata_table) based on count_table (ASV_table).
In Heatmap was also included the information about Disease/Health from that sample.

To generate this heatmaps, the data from count table were transformed by log1p(x) function, to better visualization.
