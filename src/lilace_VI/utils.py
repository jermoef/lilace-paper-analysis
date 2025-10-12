# lilace_VI/utils.py

def get_scores_df(lilace_VI_fit_df):
    columns_to_drop = [col for col in lilace_VI_fit_df.columns if col.startswith("c_")]
    columns_to_drop += ["rep", "n_counts"]
    score_df = lilace_VI_fit_df.drop(columns_to_drop, axis=1).drop_duplicates()
    return score_df

