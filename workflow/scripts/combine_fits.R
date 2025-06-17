library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
num_dfs <- length(args)-1
merged_df <- readRDS(args[1])
print(args[1])
print(head(merged_df))
merged_df$exp <- as.character(merged_df$exp)
merged_df$position <- as.numeric(merged_df$position)
for (i in 2:num_dfs) {
    df <- readRDS(args[i])
    df$exp <- as.character(df$exp)
    df$position <- as.numeric(df$position)
    print(head(df))
    intersect_cols <- intersect(colnames(merged_df), colnames(df))
    print(intersect_cols)
    for (j in 1:length(intersect_cols)) {
        col <- intersect_cols[j]
        if (any(sort(merged_df[[col]]) != sort(df[[col]]))) {
            print(paste("NOT MATCHING:", args[i], col))
        } else {
            print(paste("matching:", args[i], col))
        }
    }
    print(head(merged_df))
    print(nrow(merged_df))
    print(nrow(df))
    merged_df <- left_join(merged_df, df)
}

outfile <- args[length(args)]
saveRDS(merged_df, outfile)
