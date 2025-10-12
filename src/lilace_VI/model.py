# lilace_VI/model.py

import torch
import pyro

import pyro.distributions as dist
from torch.distributions import constraints


def lilace_VI(y, y_syn, K, vMAPp, syn_rep_index, rep_index, v_index,
              n_reps, n_pos, n_variants, n_syn, n_counts, n_counts_syn):
    std_norm = dist.Normal(0, 1)
    
    # priors
    with pyro.plate("replicate", n_reps):
        q = pyro.sample("q", dist.Dirichlet(torch.ones(K)))
        a = pyro.sample("a", dist.HalfNormal(0.1, 10)) # proxy for flat uniform prior TODO: consider param formulation
        b = pyro.sample("b", dist.HalfNormal(0.1, 10))

    with pyro.plate("positions", n_pos):
        theta = pyro.sample("theta", dist.Normal(0.,1.))
        sigma = pyro.sample("sigma", dist.InverseGamma(1.,1.))

    with pyro.plate("synonymous", n_syn):
        theta_syn = pyro.sample("theta_syn", dist.Normal(0.,1.))
    sigma_syn = pyro.sample("sigma_syn", dist.InverseGamma(1.,1.))

    with pyro.plate("variants", n_variants) as v:
        z = pyro.sample("z", dist.Normal(0.,1.))
        is_syn = (vMAPp == 0)
        syn_indices = is_syn.nonzero().squeeze(-1)
        nonsyn_indices = (~is_syn).nonzero().squeeze(-1)
        nonsyn_pos_indices = vMAPp[nonsyn_indices] - 1

        mu = torch.empty(n_variants)
        mu[syn_indices] = theta_syn + sigma_syn * z[syn_indices]
        mu[nonsyn_indices] = theta[nonsyn_pos_indices] + sigma[nonsyn_pos_indices] * z[nonsyn_indices]
        pyro.deterministic("mu", mu)

    # negative control module
    with pyro.plate("syn_obs_plate", len(y_syn)):
        phi_syn = a[syn_rep_index] * n_counts_syn + b[syn_rep_index]
        q_syn = q[syn_rep_index]
        pyro.sample("syn_obs", dist.DirichletMultinomial(phi_syn.unsqueeze(-1) * q_syn, validate_args=False), obs=y_syn)

    # full module
    t0 = std_norm.icdf(torch.cumsum(q[rep_index], dim=1)[:,:-1])
    with pyro.plate("obs_plate", len(y)):
        phi_obs = a[rep_index] * n_counts + b[rep_index]
        # get probabilities from mu vec
        obs_mu = mu[v_index]
        t2cdf = std_norm.cdf(t0 - obs_mu.unsqueeze(-1))
        # compute bin probabilities
        padded_t2cdf = torch.cat([torch.zeros(t2cdf.shape[0], 1), t2cdf, torch.ones(t2cdf.shape[0], 1)], dim=1)
        p0 = torch.diff(padded_t2cdf, dim=1)
        p0 = torch.clamp(p0, min=torch.finfo(p0.dtype).eps) # Add epsilon for numerical stability
        
        pyro.sample("obs", dist.DirichletMultinomial(phi_obs.unsqueeze(-1) * p0, validate_args=False), obs=y)


def lilace_VI_nopos(y, y_syn, K, vMAPp, syn_rep_index, rep_index, v_index,
              n_reps, n_pos, n_variants, n_syn, n_counts, n_counts_syn):
    std_norm = dist.Normal(0, 1)
    
    # priors
    with pyro.plate("replicate", n_reps):
        q = pyro.sample("q", dist.Dirichlet(torch.ones(K)))
        a = pyro.sample("a", dist.HalfNormal(0.1, 10)) # proxy for flat uniform prior TODO: consider param formulation
        b = pyro.sample("b", dist.HalfNormal(0.1, 10))

    with pyro.plate("variants", n_variants) as v:
        theta = pyro.sample("theta", dist.Normal(0.,1.))
        sigma = pyro.sample("sigma", dist.InverseGamma(1.,1.))
        mu = pyro.sample("mu", dist.Normal(theta, sigma))

    # negative control module
    with pyro.plate("syn_obs_plate", len(y_syn)):
        phi_syn = a[syn_rep_index] * n_counts_syn + b[syn_rep_index]
        q_syn = q[syn_rep_index]
        pyro.sample("syn_obs", dist.DirichletMultinomial(phi_syn.unsqueeze(-1) * q_syn, validate_args=False), obs=y_syn)

    # full module
    t0 = std_norm.icdf(torch.cumsum(q[rep_index], dim=1)[:,:-1])
    with pyro.plate("obs_plate", len(y)):
        phi_obs = a[rep_index] * n_counts + b[rep_index]
        # get probabilities from mu vec
        obs_mu = mu[v_index]
        t2cdf = std_norm.cdf(t0 - obs_mu.unsqueeze(-1))
        # compute bin probabilities
        padded_t2cdf = torch.cat([torch.zeros(t2cdf.shape[0], 1), t2cdf, torch.ones(t2cdf.shape[0], 1)], dim=1)
        p0 = torch.diff(padded_t2cdf, dim=1)
        p0 = torch.clamp(p0, min=torch.finfo(p0.dtype).eps) # Add epsilon for numerical stability
        
        pyro.sample("obs", dist.DirichletMultinomial(phi_obs.unsqueeze(-1) * p0, validate_args=False), obs=y)


def custom_hierarchical_guide(y, y_syn, K, vMAPp, syn_rep_index, rep_index, v_index,
                              n_reps, n_pos, n_variants, n_syn, n_counts, n_counts_syn):
    with pyro.plate("replicate", n_reps):
        q_concentration_param = pyro.param("q_concentration_param", torch.ones(n_reps, K),
                                           constraint=constraints.positive)
        pyro.sample("q", dist.Dirichlet(q_concentration_param))

        a_loc = pyro.param("a_loc", torch.zeros(n_reps))
        a_scale = pyro.param("a_scale", torch.ones(n_reps), constraint=constraints.positive)
        pyro.sample("a", dist.LogNormal(a_loc, a_scale))

        b_loc = pyro.param("b_loc", torch.zeros(n_reps))
        b_scale = pyro.param("b_scale", torch.ones(n_reps), constraint=constraints.positive)
        pyro.sample("b", dist.LogNormal(b_loc, b_scale))

    with pyro.plate("positions", n_pos):
        theta_loc = pyro.param("theta_loc", torch.zeros(n_pos))
        theta_scale = pyro.param("theta_scale", torch.ones(n_pos), constraint=constraints.positive)
        pyro.sample("theta", dist.Normal(theta_loc, theta_scale))

        # For the positive sigma, we again use a LogNormal guide.
        sigma_loc = pyro.param("sigma_loc", torch.zeros(n_pos))
        sigma_scale = pyro.param("sigma_scale", torch.ones(n_pos), constraint=constraints.positive)
        pyro.sample("sigma", dist.LogNormal(sigma_loc, sigma_scale))
    
    with pyro.plate("synonymous", n_syn):
        theta_syn_loc = pyro.param("theta_syn_loc", torch.zeros(n_syn))
        theta_syn_scale = pyro.param("theta_syn_scale", torch.ones(n_syn), constraint=constraints.positive)
        pyro.sample("theta_syn", dist.Normal(theta_syn_loc, theta_syn_scale))

    sigma_syn_loc = pyro.param("sigma_syn_loc", torch.tensor(0.0))
    sigma_syn_scale = pyro.param("sigma_syn_scale", torch.tensor(1.0), constraint=constraints.positive)
    pyro.sample("sigma_syn", dist.LogNormal(sigma_syn_loc, sigma_syn_scale))

    with pyro.plate("variants", n_variants):
        z_loc = pyro.param("z_loc", torch.zeros(n_variants))
        z_scale = pyro.param("z_scale", torch.ones(n_variants), constraint=constraints.positive)
        pyro.sample("z", dist.Normal(z_loc, z_scale))