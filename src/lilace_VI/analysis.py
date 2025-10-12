# lilace_VI/analysis.py
import torch
import pyro
import numpy as np

from .model import lilace_VI, lilace_VI_nopos

def get_posterior_summary(df, best_params_path, data_args, guide_type, n_batch_samples = 500, n_batches = 10, nopos = False):
    if nopos:
        model = lilace_VI_nopos
    else:
        model = lilace_VI
    
    pyro.clear_param_store() # Clear any lingering params
    pyro.get_param_store().load(best_params_path)

    n_samples = n_batch_samples * n_batches
    all_mu_samples = []
    if guide_type == "normal":
        final_guide = pyro.infer.autoguide.AutoNormal(model)
    elif guide_type == "multivariate_normal":
        final_guide = pyro.infer.autoguide.AutoMultivariateNormal(model)
    else:
        raise ValueError("guide_type not recognized")
    
    for i in range(n_batches):
        print(f"\nSampling posterior batch {i + 1}/{n_batches} ---")
        predictive = pyro.infer.Predictive(model, guide=final_guide, num_samples=n_batch_samples, return_sites=("mu",))
        posterior_samples = predictive(**data_args)
        all_mu_samples.append(posterior_samples['mu'])
    mu_samples = torch.cat(all_mu_samples, dim=0)
    mu_samples = mu_samples.numpy()

    mu_samples.mean(axis=0)
    mu_means = mu_samples.mean(axis=0)
    mu_stds = mu_samples.std(axis=0)

    p_less = np.mean(mu_samples < 0, axis=0)
    p_great = np.mean(mu_samples > 0, axis=0)
    mu_lfsr = np.minimum(p_less, p_great)

    is_syn = (data_args["vMAPp"] == 0)
    syn_indices = is_syn.nonzero().squeeze(-1)
    syn_mu = mu_samples[:, syn_indices]
    syn_mean = np.mean(syn_mu)


    n_syn_variants = syn_mu.shape[1]
    random_syn_indices = np.random.randint(0, n_syn_variants, size=mu_samples.shape)
    random_syn_samples = syn_mu[np.arange(n_samples)[:, np.newaxis], random_syn_indices]
    mu2 = mu_samples - random_syn_samples
    syn_recalibrated_mu_mean = np.mean(mu2, axis=0)
    syn_recalibrated_mu_sd = np.std(mu2, axis=0, ddof=1) # ddof=1 for sample std dev, like R's sd()
    p_less_recalibrated = np.mean(mu2 <= 0, axis=0)
    p_great_recalibrated = np.mean(mu2 >= 0, axis=0)
    syn_recalibrated_mu_lfsr = np.minimum(p_less_recalibrated, p_great_recalibrated)
    if nopos:
        df["lilace_nopos_VI_" + guide_type + "_unrecal_mu_mean"] = mu_means[df["v_index"].values]
        df["lilace_nopos_VI_" + guide_type + "_unrecal_mu_sd"] = mu_stds[df["v_index"].values]
        df["lilace_nopos_VI_" + guide_type + "_unrecal_mu_lfsr"] = mu_lfsr[df["v_index"].values]
        df["lilace_nopos_VI_" + guide_type + "_mu_mean"] = syn_recalibrated_mu_mean[df["v_index"].values]
        df["lilace_nopos_VI_" + guide_type + "_mu_sd"] = syn_recalibrated_mu_sd[df["v_index"].values]
        df["lilace_nopos_VI_" + guide_type + "_mu_lfsr"] = syn_recalibrated_mu_lfsr[df["v_index"].values]
    else:
        df["lilace_VI_" + guide_type + "_unrecal_mu_mean"] = mu_means[df["v_index"].values]
        df["lilace_VI_" + guide_type + "_unrecal_mu_sd"] = mu_stds[df["v_index"].values]
        df["lilace_VI_" + guide_type + "_unrecal_mu_lfsr"] = mu_lfsr[df["v_index"].values]
        df["lilace_VI_" + guide_type + "_mu_mean"] = syn_recalibrated_mu_mean[df["v_index"].values]
        df["lilace_VI_" + guide_type + "_mu_sd"] = syn_recalibrated_mu_sd[df["v_index"].values]
        df["lilace_VI_" + guide_type + "_mu_lfsr"] = syn_recalibrated_mu_lfsr[df["v_index"].values]

    return df