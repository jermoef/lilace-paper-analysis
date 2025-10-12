# lilace_VI/data.py
import torch
import numpy as np
import pandas as pd

def preprocess_data(df):
    if "hgvs" in df.columns:
        var_col = "hgvs"
    else:
        var_col = "variant"
    count_cols = [col for col in df if col.startswith('c_')]
    K = len(count_cols)
    y = torch.tensor(df.loc[:,count_cols].values)
    y_syn = torch.tensor(df.loc[df["type"] == "synonymous",count_cols].values)
    n_reps = np.unique(df["rep"]).size
    n_pos = np.unique(df["position"]).size
    n_variants = np.unique(df[var_col]).size
    n_counts = torch.tensor(df["n_counts"].values)
    n_counts_syn = torch.tensor(df["n_counts"].values[df["type"] == "synonymous"])

    df["v_index"] = pd.factorize(df[var_col])[0]
    df["recoded_position"] = df["position"]
    df.loc[df["type"]=="synonymous", "recoded_position"] = 0
    vMAPp = torch.tensor(pd.Series(dict(zip(df["v_index"], df["recoded_position"]))).rank(method="dense").astype(int) - 1)
    n_syn = np.count_nonzero(vMAPp == 0)
    # syn_rep_nums = (df.loc[df["type"] == "synonymous", "rep"]).astype(str).str.replace("R", "").astype(int)
    # syn_rep_index = torch.tensor((syn_rep_nums - 1).values)
    syn_rep_index = torch.tensor(((df.loc[df["type"] == "synonymous", "rep"]).astype(str).str.replace("R", "").astype(int).rank(method="dense").astype(int)-1).values)
    # rep_nums = df["rep"].astype(str).str.replace("R", "").astype(int)
    # rep_index = torch.tensor((rep_nums - 1).values)
    rep_index = torch.tensor((df["rep"].astype(str).str.replace("R", "").astype(int).rank(method="dense").astype(int)-1).values)

    v_index = torch.tensor(df["v_index"].values)

    data_args = {
        "y": y, "y_syn": y_syn, "K": K, "vMAPp": vMAPp, "syn_rep_index": syn_rep_index, "rep_index": rep_index,
        "v_index": v_index, "n_reps": n_reps, "n_pos": n_pos, "n_variants": n_variants, "n_syn": n_syn,
        "n_counts": n_counts, "n_counts_syn": n_counts_syn
    }
    return data_args