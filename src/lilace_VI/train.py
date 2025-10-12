# lilace_VI/train.py

import torch
import pyro
import logging

from .model import lilace_VI, lilace_VI_nopos


def train_model(data_args, best_params_path, guide_type = "normal", lr = 0.001, n_steps=5000, n_restarts = 1, nopos=False):
    if nopos:
        model = lilace_VI_nopos
    else:
        model = lilace_VI

    best_loss = float('inf')
    STEP_INTERVAL = 100

    for i in range(n_restarts):
        print(f"\n--- Restart {i + 1}/{n_restarts} ---")
        pyro.clear_param_store()

        if guide_type == "normal":
            auto_guide = pyro.infer.autoguide.AutoNormal(model)
        elif guide_type == "multivariate_normal":
            auto_guide = pyro.infer.autoguide.AutoMultivariateNormal(model)
        else:
            raise ValueError("guide_type not recognized")

        adam = pyro.optim.Adam({"lr": lr})
        elbo = pyro.infer.Trace_ELBO()
        svi = pyro.infer.SVI(model, auto_guide, adam, elbo)


        losses = []
        for step in range(n_steps):  # Consider running for more steps
            loss = svi.step(**data_args)
            losses.append(loss)
            if step % 500 == 0:
                logging.info("Elbo loss: {}".format(loss))

            # if (step + 1) % STEP_INTERVAL == 0:
            #     scheduler.step()
        
        final_loss = losses[-1]
        print(f"Final loss for restart {i + 1}: {final_loss:.2f}")

        if final_loss < best_loss:
            best_loss = final_loss
            pyro.get_param_store().save(best_params_path)
    return best_params_path, losses



