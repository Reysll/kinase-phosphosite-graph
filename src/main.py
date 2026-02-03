# src/main.py
from __future__ import annotations

import os

from src.builders.kinase_substrate import build_kinase_substrate_graph
from src.sanity_checks import check_site_collisions
from src.config import KINASE_SUBSTRATE_FILE, OUTPUT_DIR


def main() -> None:
    nodes_df, edges_df = build_kinase_substrate_graph(KINASE_SUBSTRATE_FILE)

    # Sanity checks, including your requirement:
    # "Same site label (ex: S298) can exist on multiple proteins"
    check_site_collisions(nodes_df)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    nodes_df.to_csv(os.path.join(OUTPUT_DIR, "nodes.csv.gz"), index=False, compression="gzip")
    edges_df.to_csv(os.path.join(OUTPUT_DIR, "edges.csv.gz"), index=False, compression="gzip")

    print("Wrote outputs:")
    print(f"  {os.path.join(OUTPUT_DIR, 'nodes.csv.gz')}")
    print(f"  {os.path.join(OUTPUT_DIR, 'edges.csv.gz')}")


if __name__ == "__main__":
    main()
