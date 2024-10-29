import numpy as np 
import argparse
import pandas as pd 
import natsort 
import os

# Change working directory to reference folder
os.chdir("/ibex/user/iribarxm/spatial_transcriptomics/bidcell/BIDCell-main/breast_cancer_data/sc_references")

def main(config):
    #
    ref_df = pd.read_csv(config.fp_ref, index_col=0)
    n_genes = ref_df.shape[1] - 3
    print("Ref data shape", ref_df.shape)
    #
    cell_types = ref_df["cell_type"].tolist()
    cell_types = natsort.natsorted(list(set(cell_types)))
    print(cell_types)
    n_cell_types = len(cell_types)
    #
    ref_expr = ref_df.iloc[:, :n_genes].to_numpy()
    gene_names = ref_df.columns[:n_genes]
    #
    # Find genes with expressions in bottom 10% percentile for every ref cell type
    pct_10 = np.percentile(ref_expr, 10, axis=1, keepdims=True)
    pct_10 = np.tile(pct_10, (1, n_genes))
    low_expr_true = np.zeros(pct_10.shape)
    low_expr_true[ref_expr <= pct_10] = 1
    #
    # Find overlap for different ref samples of the same cell type
    ct_idx = ref_df["ct_idx"].to_numpy()
    low_expr_true_agg = np.zeros((n_cell_types, n_genes))
    for ct in range(n_cell_types):
        rows = np.where(ct_idx == ct)[0]
        low_expr_true_ct = low_expr_true[rows]
        low_expr_true_agg[ct, :] = np.prod(low_expr_true_ct, axis=0)
    #
    # Set overlaps to 0
    overlaps = np.sum(low_expr_true_agg, 0)
    too_many = np.where(overlaps > config.max_overlaps_neg)[0]
    low_expr_true_agg[:, too_many] = 0
    #
    df_neg = pd.DataFrame(low_expr_true_agg, index=cell_types, columns=gene_names)
    #
    # Find genes with expressions in top 90% percentile for every ref cell type
    pct_90 = np.percentile(ref_expr, 90, axis=1, keepdims=True)
    pct_90 = np.tile(pct_90, (1, n_genes))
    high_expr_true = np.zeros(pct_90.shape)
    high_expr_true[ref_expr >= pct_90] = 1
    #
    # Find overlap for different ref samples of the same cell type
    ct_idx = ref_df["ct_idx"].to_numpy()
    high_expr_true_agg = np.zeros((n_cell_types, n_genes))
    for ct in range(n_cell_types):
        rows = np.where(ct_idx == ct)[0]
        high_expr_true_ct = high_expr_true[rows]
        high_expr_true_agg[ct, :] = np.prod(high_expr_true_ct, axis=0)
    #
    # Set overlaps to 0
    overlaps = np.sum(high_expr_true_agg, 0)
    too_many = np.where(overlaps > config.max_overlaps_pos)[0]
    high_expr_true_agg[:, too_many] = 0
    #
    df_pos = pd.DataFrame(high_expr_true_agg, index=cell_types, columns=gene_names)
    #
    df_pos.to_csv(config.fp_pos)
    df_neg.to_csv(config.fp_neg)
    #
    print("Done")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    #
    parser.add_argument(
        "--fp_ref",
        default="/ibex/user/iribarxm/spatial_transcriptomics/bidcell/BIDCell-main/breast_cancer_data/sc_references/sc_breast.csv",
        type=str,
        help="ref data",
    )
    #
    parser.add_argument(
        "--max_overlaps_pos",
        default=4,
        type=int,
        help="no more than n cell types can have the same markers",
    )
    #
    parser.add_argument(
        "--max_overlaps_neg",
        default=15,
        type=int,
        help="no more than n cell types can have the same markers",
    )
    #
    parser.add_argument(
        "--fp_pos", default="markers_pos.csv", type=str, help="positive markers"
    )
    parser.add_argument(
        "--fp_neg", default="markers_neg.csv", type=str, help="negative markers"
    )
    #
    config = parser.parse_args()
    main(config)