from __future__ import annotations

from src_prediction.config import LIVER_POSITIVE_EDGES_OUT, PRED_OUTPUTS_DIR
from src_prediction.io_utils import read_csv_auto, write_csv_gz


def main() -> None:
    positive_edges = read_csv_auto(LIVER_POSITIVE_EDGES_OUT)

    kinase_counts = positive_edges.groupby("site_node_id")["kinase_node_id"].transform("count")
    multi_kinase_edges = positive_edges.loc[kinase_counts > 1].copy().reset_index(drop=True)
    multi_kinase_edges["trial_index"] = range(1, len(multi_kinase_edges) + 1)

    out_path = PRED_OUTPUTS_DIR / "frozen_trials_multi_kinase.csv.gz"
    write_csv_gz(multi_kinase_edges, out_path)

    n_sites = multi_kinase_edges["site_node_id"].nunique()
    print(f"Multi-kinase LOO trials: {len(multi_kinase_edges):,} edges across {n_sites:,} sites")
    print(f"Wrote: {out_path}")


if __name__ == "__main__":
    main()
