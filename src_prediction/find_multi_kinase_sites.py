# src_prediction/find_multi_kinase_sites.py

from pathlib import Path
import pandas as pd


INPUT_TSV = Path("data/raw/Kinase_Substrate_Dataset.tsv")
OUTPUT_DIR = Path("outputs_prediction")


def normalize_site_token(site: str) -> str:
    """
    Normalize site tokens like:
    S52 -> S52
    T308 -> T308
    Y204 -> Y204
    Ser52 -> S52
    Thr308 -> T308
    Tyr204 -> Y204
    """
    s = str(site).strip()

    replacements = {
        "Ser": "S",
        "SER": "S",
        "Thr": "T",
        "THR": "T",
        "Tyr": "Y",
        "TYR": "Y",
    }
    for old, new in replacements.items():
        s = s.replace(old, new)

    s = s.replace(" ", "")
    return s.upper()


def is_single_site(site: str) -> bool:
    """
    Keep only rows that appear to contain a single phosphosite.
    Exclude things like:
    S10;S12
    S10,S12
    S10/S12
    """
    s = str(site).strip()
    if s == "" or s.lower() == "nan":
        return False

    bad_separators = [";", ",", "/"]
    return not any(sep in s for sep in bad_separators)


def build_site_node_id(substrate_gene: str, site_token: str) -> str:
    return f"SITE:{str(substrate_gene).upper()}-{site_token}"


def build_kinase_node_id(kinase: str) -> str:
    return f"PROTEIN:{str(kinase).upper()}"


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if not INPUT_TSV.exists():
        raise FileNotFoundError(
            f"Could not find input file: {INPUT_TSV}\n"
            f"Make sure Kinase_Substrate_Dataset.tsv is located there."
        )

    print(f"Reading: {INPUT_TSV}")
    df = pd.read_csv(INPUT_TSV, sep="\t", dtype=str).fillna("")

    required_cols = [
        "KINASE",
        "KIN_ORGANISM",
        "SUB_GENE",
        "SUB_ORGANISM",
        "SUB_MOD_RSD",
    ]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"Missing required columns: {missing}\n"
            f"Found columns: {list(df.columns)}"
        )

    print(f"Total rows in raw file: {len(df):,}")

    # Keep only human kinase and human substrate rows
    human_df = df[
        df["KIN_ORGANISM"].str.lower().eq("human")
        & df["SUB_ORGANISM"].str.lower().eq("human")
    ].copy()

    print(f"Rows after human-human filter: {len(human_df):,}")

    # Keep only single-site rows
    human_df = human_df[human_df["SUB_MOD_RSD"].map(is_single_site)].copy()

    print(f"Rows after single-site filter: {len(human_df):,}")

    # Normalize identifiers
    human_df["site_token"] = human_df["SUB_MOD_RSD"].map(normalize_site_token)
    human_df["site_node_id"] = human_df.apply(
        lambda r: build_site_node_id(r["SUB_GENE"], r["site_token"]),
        axis=1,
    )
    human_df["kinase_node_id"] = human_df["KINASE"].map(build_kinase_node_id)

    # Keep only the columns we need and drop exact duplicates
    pairs = human_df[
        [
            "KINASE",
            "SUB_GENE",
            "SUB_MOD_RSD",
            "site_token",
            "site_node_id",
            "kinase_node_id",
            "CST_CAT#",
        ]
    ].drop_duplicates()

    print(f"Unique human kinase-site pairs: {len(pairs):,}")

    # Group by site to find all reported kinases per site
    grouped = (
        pairs.groupby("site_node_id")
        .agg(
            substrate_gene=("SUB_GENE", "first"),
            site_token=("site_token", "first"),
            original_site_value=("SUB_MOD_RSD", "first"),
            n_true_kinases=("kinase_node_id", lambda s: len(sorted(set(s)))),
            true_kinase_node_ids=("kinase_node_id", lambda s: ";".join(sorted(set(s)))),
            true_kinase_symbols=("KINASE", lambda s: ";".join(sorted(set(s)))),
            cst_cat_values=("CST_CAT#", lambda s: ";".join(sorted(set(x for x in s if str(x).strip() != "")))),
        )
        .reset_index()
    )

    grouped = grouped.sort_values(
        by=["n_true_kinases", "site_node_id"],
        ascending=[False, True]
    )

    # Save all sites
    all_sites_path = OUTPUT_DIR / "all_human_sites_kinase_counts.csv"
    grouped.to_csv(all_sites_path, index=False)

    # Save only multi-kinase sites
    multi = grouped[grouped["n_true_kinases"] > 1].copy()
    multi_path = OUTPUT_DIR / "multi_kinase_sites.csv"
    multi.to_csv(multi_path, index=False)

    print()
    print(f"Total unique human sites: {len(grouped):,}")
    print(f"Sites with >1 reported kinase: {len(multi):,}")
    print(f"Wrote: {all_sites_path}")
    print(f"Wrote: {multi_path}")

    if len(multi) > 0:
        print()
        print("Top multi-kinase sites:")
        print(
            multi[
                [
                    "site_node_id",
                    "n_true_kinases",
                    "true_kinase_symbols",
                ]
            ].head(20).to_string(index=False)
        )


if __name__ == "__main__":
    main()