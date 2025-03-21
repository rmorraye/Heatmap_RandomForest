# Heatmap_RandomForest

This code get the result table from Random Forest Analysis, of a Dental Microbiota Analysis (16S), and generate Heatmaps for each study (in Metadata_table) based on counts (ASV_table).
In Heatmap was also included the information about Disease/Health from that sample.

To generate this heatmaps, the data from count table were summed and transformed by log1p(x) function, to better visualization.
