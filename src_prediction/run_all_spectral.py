"""
Run all four networks with SpectralEmbeddingStrategy + category ranking.
Executes sequentially: generic -> control -> cancer -> liver FC.
Each network takes ~45-55 minutes (6,909 trials, 8 workers).
"""

from __future__ import annotations

import time

from src_prediction.run_multi_kinase_generic_spectral import main as run_generic
from src_prediction.run_multi_kinase_control_spectral import main as run_control
from src_prediction.run_multi_kinase_cancer_spectral import main as run_cancer
from src_prediction.run_multi_kinase_liver_spectral import main as run_liver


def main() -> None:
    t0 = time.time()

    print("\n" + "=" * 70)
    print("SPECTRAL RUN 1/4: GENERIC")
    print("=" * 70)
    run_generic()

    print("\n" + "=" * 70)
    print("SPECTRAL RUN 2/4: CONTROL")
    print("=" * 70)
    run_control()

    print("\n" + "=" * 70)
    print("SPECTRAL RUN 3/4: CANCER")
    print("=" * 70)
    run_cancer()

    print("\n" + "=" * 70)
    print("SPECTRAL RUN 4/4: LIVER FC")
    print("=" * 70)
    run_liver()

    elapsed = (time.time() - t0) / 60
    print(f"\nAll four spectral runs complete. Total wall time: {elapsed:.1f} min")


if __name__ == "__main__":
    main()
