import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

import logging
import time as tm

import argparse

from lilace_VI import data, train, analysis, utils

def main():
    parser = argparse.ArgumentParser(description="Run the complete Lilace VI workflow.")
    # --- Required Arguments ---
    parser.add_argument(
        "--input_file",
        type=str,
        required=True,
        help="Path to the input data CSV file formatted in Lilace format (see R packaage)."
    )
    parser.add_argument(
        "--output_file",
        type=str,
        required=True,
        help="Path to save the variant scores TSV."
    )

    parser.add_argument("--guide_type", type=str, default="normal", help="Guide type: 'normal' or 'multivariate_normal'. Default: 'normal'")
    parser.add_argument("--lr", type=float, default=0.001, help="Learning rate for the Adam optimizer. Default: 0.001")
    parser.add_argument("--n_steps", type=float, default=5000, help="Number of training steps. Default: 5000")
    parser.add_argument("--restarts", type=int, default=3, help="Number of random restarts for training. Default: 3")
    parser.add_argument("--n_batch_samples", type=int, default=500, help="Number of posterior samples to draw in each sampling batch. Default: 500")
    parser.add_argument("--n_batches", type=int, default=10, help="Number of batches of posterior samples to draw. Default: 10")
    parser.add_argument("--plot_dir", type=str, default="VI_plots", help="Path to save the loss curve plot. Default: VI_plots")
    parser.add_argument("--nopos", action='store_true', help="Path to save the loss curve plot. Default: VI_plots")

    args = parser.parse_args()
    lilace_df = pd.read_csv(args.input_file, header=0)

    data_args = data.preprocess_data(lilace_df)
    print(data_args)
    best_params_path = args.plot_dir + "/best_model_params.pt"
    best_params_path, losses = train.train_model(
        data_args,
        best_params_path,
        args.guide_type,
        args.lr,
        args.n_steps,
        args.restarts,
        args.nopos
        )
    
    plt.figure(figsize=(5, 2))
    plt.plot(losses)
    plt.xlabel("SVI step")
    plt.ylabel("ELBO loss")
    plt.savefig(args.plot_dir + "/loss_curve_" + args.guide_type + ".png")

    lilace_df_with_VI = analysis.get_posterior_summary(
        lilace_df,
        best_params_path,
        data_args,
        args.guide_type,
        args.n_batch_samples,
        args.n_batches,
        args.nopos
    )

    fitted_df_file = args.plot_dir + "/fitted_df_VI.csv"
    lilace_df_with_VI.to_csv(fitted_df_file, index=False)
    print(f"Fitted df saved to {fitted_df_file}")
    # save variant scores to variant_scores_VI.tsv
    print(lilace_df_with_VI.columns)
    scores_table = utils.get_scores_df(lilace_df_with_VI)
    scores_table.to_csv(args.output_file, sep='\t', index=False)
    print(f"Variant scores saved to {args.output_file}")

if __name__ == "__main__":
    main()

