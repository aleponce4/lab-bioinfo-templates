"""
Simulate haplotype composition data for template 14:
  - data/haplotype_frequencies.csv — per-sample haplotype frequencies
  - data/manifest.csv             — sample metadata (sample_id, dpi, replicate, route)

Run from the template folder:
    python data/simulate_data.py
"""

import csv
import os
import random

SEED = 42
random.seed(SEED)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TEMPLATE_DIR = os.path.dirname(SCRIPT_DIR)
DATA_DIR = os.path.join(TEMPLATE_DIR, "data")
os.makedirs(DATA_DIR, exist_ok=True)

DPI_VALUES = [3, 5]
ROUTES = ["Intranasal", "Subcutaneous"]
REPS_PER_DPI_ROUTE = 3
BASE_HAPLOTYPES = [
    {"name": "Wild-type",       "snps": (), "base_freq": 0.55},
    {"name": "Hap_Cluster_1",   "snps": (241, 3401, 9102), "base_freq": 0.18},
    {"name": "Hap_Cluster_2",   "snps": (5200, 6780), "base_freq": 0.12},
    {"name": "Hap_Cluster_3",   "snps": (2400, 3401, 10800), "base_freq": 0.08},
    {"name": "Hap_Cluster_4",   "snps": (150, 10100), "base_freq": 0.05},
    {"name": "Hap_Cluster_5",   "snps": (3401,), "base_freq": 0.02},
]


def generate():
    rows = []
    manifest_rows = []
    sample_idx = 1

    for route in ROUTES:
        route_factor = 0.7 if route == "Intranasal" else 1.0
        for dpi_val in DPI_VALUES:
            dpi_factor = 0.8 if dpi_val == 3 else 1.0
            for rep in range(1, REPS_PER_DPI_ROUTE + 1):
                sample_id = f"S{sample_idx}"
                sample_idx += 1

                manifest_rows.append({
                    "sample_id": sample_id,
                    "dpi": dpi_val,
                    "replicate": rep,
                    "route": route,
                })

                freqs = []
                for hap in BASE_HAPLOTYPES:
                    f = hap["base_freq"] * route_factor * dpi_factor * random.gauss(1.0, 0.15)
                    f = max(0.005, f)
                    freqs.append(f)

                # Normalise so each sample's haplotype frequencies sum to 1.0.
                # Round first, then absorb the rounding residual into the largest
                # haplotype so the written values sum to exactly 1.0000.
                total = sum(freqs)
                rounded = [round(f / total, 4) for f in freqs]
                biggest = max(range(len(rounded)), key=rounded.__getitem__)
                rounded[biggest] = round(rounded[biggest] + (1.0 - sum(rounded)), 4)

                for hap, freq in zip(BASE_HAPLOTYPES, rounded):
                    rows.append({
                        "sample_id": sample_id,
                        "haplotype": hap["name"],
                        "snp_positions": ",".join(str(s) for s in hap["snps"]) if hap["snps"] else "",
                        "frequency": freq,
                    })

    with open(os.path.join(DATA_DIR, "haplotype_frequencies.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    with open(os.path.join(DATA_DIR, "manifest.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(manifest_rows[0].keys()))
        w.writeheader()
        w.writerows(manifest_rows)
    print(f"  Wrote data/haplotype_frequencies.csv ({len(rows)} rows) "
          f"and data/manifest.csv ({len(manifest_rows)} samples)")


if __name__ == "__main__":
    generate()
