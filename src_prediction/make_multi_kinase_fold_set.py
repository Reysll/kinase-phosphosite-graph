# src_prediction/make_multi_kinase_fold_set.py

from pathlib import Path
import pandas as pd


INPUT_TSV = Path("data/raw/Kinase_Substrate_Dataset.tsv")
MULTI_SITES_CSV = Path("outputs_prediction/multi_kinase_sites.csv")
OUTPUT_DIR = Path("outputs_prediction")


def normalize_site_token(site: str) -> str:
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
    s = str(site).strip()
    if s == "" or s.lower() == "nan":
        return False
    return not any(sep in s for sep in [";", ",", "/"])


def build_site_node_id(substrate_gene: str, site_token: str) -> str:
    return f"SITE:{str(substrate_gene).upper()}-{site_token}"


def build_kinase_node_id(kinase: str) -> str:
    return f"PROTEIN:{str(kinase).upper()}"


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if not INPUT_TSV.exists():
        raise FileNotFoundError(f"Missing input TSV: {INPUT_TSV}")

    if not MULTI_SITES_CSV.exists():
        raise FileNotFoundError(
            f"Missing multi-site file: {MULTI_SITES_CSV}\n"
            f"Run python -m src_prediction.find_multi_kinase_sites first."
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
        raise ValueError(f"Missing required columns: {missing}")

    df = df[
        df["KIN_ORGANISM"].str.lower().eq("human")
        & df["SUB_ORGANISM"].str.lower().eq("human")
    ].copy()

    df = df[df["SUB_MOD_RSD"].map(is_single_site)].copy()

    df["site_token"] = df["SUB_MOD_RSD"].map(normalize_site_token)
    df["site_node_id"] = df.apply(
        lambda r: build_site_node_id(r["SUB_GENE"], r["site_token"]),
        axis=1,
    )
    df["kinase_node_id"] = df["KINASE"].map(build_kinase_node_id)

    multi = pd.read_csv(MULTI_SITES_CSV, dtype=str).fillna("")
    target_sites = set(multi["site_node_id"])

    df = df[df["site_node_id"].isin(target_sites)].copy()

    df = df.drop_duplicates(subset=["kinase_node_id", "site_node_id"])

    folds = df[["kinase_node_id", "site_node_id"]].copy()
    folds.insert(0, "fold_index", range(1, len(folds) + 1))
    folds["held_out_relation"] = "phosphorylates"

    folds = folds[["fold_index", "kinase_node_id", "site_node_id", "held_out_relation"]]

    out_path = OUTPUT_DIR / "frozen_folds_multi_kinase_sites.csv.gz"
    folds.to_csv(out_path, index=False, compression="gzip")

    print()
    print(f"Multi-kinase kinase-site pairs written: {len(folds):,}")
    print(f"Unique multi-kinase sites included: {folds['site_node_id'].nunique():,}")
    print(f"Wrote: {out_path}")
    print()
    print("Preview:")
    print(folds.head(10).to_string(index=False))


if __name__ == "__main__":
    main()