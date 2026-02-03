from __future__ import annotations
import pandas as pd

def check_site_collisions(nodes_df: pd.DataFrame) -> None:
    sites = nodes_df[nodes_df["node_type"] == "site"].copy()

    # Split "PROTEIN-S298" into protein and site_label
    split = sites["node_id"].str.rsplit(pat="-", n=1, expand=True)
    sites["protein"] = split[0]
    sites["site_label"] = split[1]

    shared_counts = sites.groupby("site_label")["protein"].nunique()
    shared_labels = shared_counts[shared_counts > 1].sort_values(ascending=False)

    print(f"Total site nodes: {len(sites)}")
    print(f"Site labels shared across proteins (expected): {len(shared_labels)}")

    if len(shared_labels) > 0:
        print("Examples (site_label -> #proteins):")
        for label, nprot in shared_labels.head(10).items():
            print(f"  {label} -> {nprot}")

    dup_within = sites.groupby(["protein", "site_label"]).size()
    dup_within = dup_within[dup_within > 1]

    if not dup_within.empty:
        print("WARNING: duplicate site nodes within same protein detected!")
        print(dup_within.head(10))
    else:
        print("No duplicate sites within proteins.")
