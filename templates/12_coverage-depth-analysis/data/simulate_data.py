"""
Simulate genome coverage depth data for template 11:
  - data/gene_coords.csv     — genome annotation
  - data/coverage.csv        — per-position depth per sample
  - data/manifest.csv        — sample metadata
"""

import csv, random, math, os

SEED = 42
random.seed(SEED)

GENOME_LEN = 11447
GENES = {
    "nsp1": (1, 1800), "nsp2": (1801, 4500), "nsp3": (4501, 5700),
    "nsp4": (5701, 7500), "capsid": (7501, 8300), "E3": (8301, 9000),
    "E2": (9001, 10200), "6K": (10201, 10400), "E1": (10401, 11400)
}

DPI_REPS = {1: 3, 3: 3, 5: 3}
BASE_DEPTH = 8000
AMPLICON_DROPS = [(2000, 2200), (4300, 4450), (7000, 7150)]


def generate():
    rows = []
    for dpi_val, n_reps in DPI_REPS.items():
        for rep in range(1, n_reps + 1):
            sample_id = f"DPI{dpi_val}_R{rep}_Lung"
            # Simulate genome coverage with amplicon dropout sites and noise
            for pos in range(1, GENOME_LEN + 1):
                depth = BASE_DEPTH + random.gauss(0, 800)
                # Amplicon dropout regions
                for ds, de in AMPLICON_DROPS:
                    if ds <= pos <= de:
                        depth *= random.uniform(0.05, 0.40)
                # 3' bias toward structural genes
                if pos > 7500:
                    depth *= random.uniform(0.85, 1.15)
                else:
                    depth *= random.uniform(0.70, 1.10)
                depth = max(20, int(depth))
                rows.append({
                    "CHROM": "VEEV",
                    "Position": pos,
                    "Depth": depth,
                    "sample_id": sample_id,
                    "dpi": dpi_val,
                    "replicate": rep
                })

    with open(os.path.join(os.path.dirname(__file__), "coverage.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    print(f"  Wrote coverage.csv ({len(rows)} rows)")

    with open(os.path.join(os.path.dirname(__file__), "gene_coords.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "start", "end"])
        for name, (s, e) in GENES.items():
            w.writerow([name, s, e])

    with open(os.path.join(os.path.dirname(__file__), "manifest.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample_id", "dpi", "replicate"])
        for dpi_val, n_reps in DPI_REPS.items():
            for rep in range(1, n_reps + 1):
                w.writerow([f"DPI{dpi_val}_R{rep}_Lung", dpi_val, rep])


if __name__ == "__main__":
    generate()
