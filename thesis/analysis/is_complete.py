"""
Completeness probe used by run_overnight.ps1.

Exits 0 if the named artefact exists AND is complete, 1 otherwise, so the
driver can skip finished work without re-running it.

Usage
-----
    python -m thesis.analysis.is_complete loo   <results.csv.gz>  [min_trials]
    python -m thesis.analysis.is_complete ranks <ranks.csv.gz>
"""

from __future__ import annotations

import sys
from pathlib import Path


def main() -> int:
    if len(sys.argv) < 3:
        return 1
    kind, path = sys.argv[1], Path(sys.argv[2])
    if not path.exists():
        return 1
    try:
        import pandas as pd
        df = pd.read_csv(path)
    except Exception:
        return 1

    if kind == "loo":
        min_trials = int(sys.argv[3]) if len(sys.argv) > 3 else 6000
        ok = len(df) >= min_trials and "adjusted_held_out_rank" in df.columns
    elif kind == "ranks":
        ok = len(df) > 0 and {"site_node_id", "rank_in_site"} <= set(df.columns)
    else:
        ok = False
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
