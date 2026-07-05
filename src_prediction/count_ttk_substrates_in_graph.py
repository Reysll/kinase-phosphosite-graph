# src_prediction/count_ttk_substrates_in_graph.py

from pathlib import Path
import pandas as pd

EDGES_FILE = Path("outputs/edges.csv.gz")
TTK = "PROTEIN:TTK"


def main():
    if not EDGES_FILE.exists():
        raise FileNotFoundError(
            f"Missing file: {EDGES_FILE}\n"
            f"Run the generic graph build first if needed."
        )

    edges = pd.read_csv(EDGES_FILE)

    required = {"source", "target", "relation"}
    missing = required - set(edges.columns)
    if missing:
        raise ValueError(f"edges file missing required columns: {sorted(missing)}")

    ttk_edges = edges[
        (edges["source"] == TTK) &
        (edges["relation"] == "phosphorylates")
    ].copy()

    n_edges = len(ttk_edges)
    n_sites = ttk_edges["target"].nunique()

    print("TTK substrate count in graph")
    print("=" * 40)
    print(f"TTK phosphorylates edges in graph: {n_edges}")
    print(f"Unique phosphosite substrates for TTK in graph: {n_sites}")
    print()

    if not ttk_edges.empty:
        print("Example TTK substrate sites:")
        print(ttk_edges["target"].drop_duplicates().head(20).to_string(index=False))

        out_path = Path("outputs_prediction/ttk_substrates_in_graph.csv")
        ttk_edges.to_csv(out_path, index=False)
        print()
        print(f"Wrote: {out_path}")


if __name__ == "__main__":
    main()