from pathlib import Path
from src.liver_network import run

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
)
