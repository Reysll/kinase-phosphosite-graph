from __future__ import annotations

from src_prediction.core.config import LIVER_POSITIVE_EDGES_OUT, PRED_OUTPUTS_DIR
from src_prediction.core.io_utils import read_csv_auto, write_csv_gz


def main() -> None:
    """Freeze ALL kinase-site positive edges as LOO trials, including single-kinase sites."""
    positive_edges = read_csv_auto(LIVER_POSITIVE_EDGES_OUT)
    all_edges = positive_edges.copy().reset_index(drop=True)
    all_edges["trial_index"] = range(1, len(all_edges) + 1)

    out_path = PRED_OUTPUTS_DIR / "frozen_trials_all.csv.gz"
    write_csv_gz(all_edges, out_path)

    n_sites = all_edges["site_node_id"].nunique()
    n_multi = (all_edges.groupby("site_node_id")["kinase_node_id"].transform("count") > 1).sum()
    n_single = len(all_edges) - n_multi
    print(f"All-sites LOO trials: {len(all_edges):,} edges across {n_sites:,} sites")
    print(f"  Multi-kinase edges: {n_multi:,}")
    print(f"  Single-kinase edges: {n_single:,}")
    print(f"Wrote: {out_path}")


if __name__ == "__main__":
    main()
