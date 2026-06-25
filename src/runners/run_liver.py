from pathlib import Path
from src.graph.liver_network import run

# Builds the liver-specific graph with ALL THREE correlation edge types embedded
# in a single output file.  The LOO evaluation scripts select which type to use
# via allowed_relations — no need to rebuild the graph for each experiment.
#
# Edge types added:
#   site_corr_fc_pos / site_corr_fc_neg       fold-change correlation (original)
#   site_corr_ctrl_pos / site_corr_ctrl_neg    control/healthy raw abundance corr
#   site_corr_cancer_pos / site_corr_cancer_neg tumor raw abundance corr
#
# All three use 80th percentile threshold and require ≥6 samples per site.

run(
    liver_xlsx=Path("data/liver/LiverCancer_ProtExp_Phospho_casecntrl.xlsx"),
    protein_sheet="ProteinExpression",
    phospho_sheet="Phosphorylation",
    generic_nodes=Path("outputs/nodes.csv.gz"),
    generic_edges=Path("outputs/edges.csv.gz"),
    out_nodes=Path("outputs/Liver_network_nodes.csv.gz"),
    out_edges=Path("outputs/Liver_network_edges.csv.gz"),
    ptmsigdb_path="data/raw/PTMsigDB.txt",
    msigdb_path="data/raw/MsigDB.txt",
    ppi_path="data/raw/high_confidence_score.csv",
    gene_map_path="data/raw/gene_protein.csv",
    ptmcode_path="data/raw/PTMcode2_associations_between_proteins.txt.gz",
    site_corr_percentile=80.0,
    add_ctrl_corr=True,
    add_cancer_corr=True,
    add_protein_fc_corr=False,
)
