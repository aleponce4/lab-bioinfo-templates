"""
Simulate haplotype composition data for template 13:
  - data/haplotype_frequencies.csv — per-sample haplotype frequencies
  - data/gene_coords.csv          — genome annotation
  - data/manifest.csv             — sample metadata
"""

import csv, random, math, os, itertools

SEED = 42
random.seed(SEED)

GENES = {
    "nsp1": (1, 1800), "nsp2": (1801, 4500), "nsp3": (4501, 5700),
    "nsp4": (5701, 7500), "capsid": (7501, 8300), "E3": (8301, 9000),
    "E2": (9001, 10200), "6K": (10201, 10400), "E1": (10401, 11400)
}

DPI_VALUES = [3, 5]
ROUTES = ["Intranasal", "Subcutaneous"]
REPS_PER_DPI_ROUTE = 3
BASE_HAPLOTYPES = [
    {"name": "Ref (Wild-type)", "snps": (), "base_freq": 0.55},
    {"name": "Hap_Cluster_1",   "snps": (241, 3401, 9102), "base_freq": 0.18},
    {"name": "Hap_Cluster_2",   "snps": (5200, 6780), "base_freq": 0.12},
    {"name": "Hap_Cluster_3",   "snps": (2400, 3401, 10800), "base_freq": 0.08},
    {"name": "Hap_Cluster_4",   "snps": (150, 10100), "base_freq": 0.05},
    {"name": "Hap_Cluster_5",   "snps": (3401,), "base_freq": 0.02},
]


def generate():
    # Write gene coords
    with open(os.path.join(os.path.dirname(__file__), "gene_coords.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "start", "end"])
        for name, (s, e) in GENES.items():
            w.writerow([name, s, e])

    rows = []
    manifest_rows = []
    sample_idx = 1

    for route in ROUTES:
        study_code = "45" if route == "Intranasal" else "47"
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
                    "study_code": study_code
                })

                freqs = []
                for hap in BASE_HAPLOTYPES:
                    f = hap["base_freq"] * route_factor * dpi_factor * random.gauss(1.0, 0.15)
                    f = max(0.005, f)
                    freqs.append(f)

                total = sum(freqs)
                for hap, freq in zip(BASE_HAPLOTYPES, freqs):
                    rows.append({
                        "sample_id": sample_id,
                        "dpi": dpi_val,
                        "replicate": rep,
                        "route": route,
                        "study_code": study_code,
                        "haplotype": hap["name"],
                        "snp_positions": ",".join(str(s) for s in hap["snps"]) if hap["snps"] else "",
                        "frequency": round(freq / total, 4)
                    })

    with open(os.path.join(os.path.dirname(__file__), "haplotype_frequencies.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)

    with open(os.path.join(os.path.dirname(__file__), "manifest.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=manifest_rows[0].keys())
        w.writeheader()
        w.writerows(manifest_rows)
    print(f"  Wrote haplotype_frequencies.csv ({len(rows)} rows) and manifest.csv")


if __name__ == "__main__":
    generate()
