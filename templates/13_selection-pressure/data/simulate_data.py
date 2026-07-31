"""
Simulate selection pressure (piN-piS) data for template 12:
  - data/delta_per_sample.csv   — piN, piS, delta per gene/sample/DPI
  - data/selection_gene_key.csv — limma & Kruskal FDR per gene
"""

import csv, random, math, os

SEED = 42
random.seed(SEED)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TEMPLATE_DIR = os.path.dirname(SCRIPT_DIR)
os.chdir(TEMPLATE_DIR)
DATA_DIR = os.path.join(TEMPLATE_DIR, "data")
os.makedirs(DATA_DIR, exist_ok=True)


GENES = ["nsp1", "nsp2", "nsp3", "nsp4", "capsid", "E3", "E2", "6K", "E1"]
DPI_VALUES = [1, 3, 5]
REPLICATES = 3

GENE_DELTA_EFFECTS = {
    "nsp1": 0.001, "nsp2": 0.003, "nsp3": 0.008, "nsp4": 0.001,
    "capsid": 0.002, "E3": 0.012, "E2": 0.010, "6K": 0.015, "E1": 0.005
}
GENE_EFFECT_DPI_BOOST = {"nsp3": (3, 2.0), "E2": (5, 2.5), "E3": (5, 1.8), "6K": (3, 1.5)}


def generate():
    rows = []
    for dpi_val in DPI_VALUES:
        for rep in range(1, REPLICATES + 1):
            sample_id = f"DPI{dpi_val}_R{rep}"
            for gene in GENES:
                base = GENE_DELTA_EFFECTS[gene]
                boost = 1.0
                for g, (boost_dpi, factor) in GENE_EFFECT_DPI_BOOST.items():
                    if g == gene and dpi_val == boost_dpi:
                        boost = factor
                delta = base * boost * random.gauss(1.0, 0.3)
                pi_n = random.gauss(0.012, 0.003)
                pi_s = pi_n - delta * random.gauss(1.0, 0.2)
                rows.append({
                    "threshold": "minfreq_0p01",
                    "sample": sample_id,
                    "dpi": f"dpi{dpi_val}",
                    "product": gene,
                    "piN": round(max(0.001, pi_n), 6),
                    "piS": round(max(0.001, pi_s), 6),
                    "delta_piN_minus_piS": round(delta, 6)
                })

    with open(os.path.join(os.path.dirname(__file__), "delta_per_sample.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)

    key_rows = []
    for gene in GENES:
        vals = [float(r["delta_piN_minus_piS"]) for r in rows if r["product"] == gene]
        if vals:
            mean_val = sum(vals) / len(vals)
            fdr = min(0.99, max(0.001, math.exp(-abs(mean_val) * 200)))
            key_rows.append({
                "threshold": "minfreq_0p01",
                "product": gene,
                "mean_delta": round(mean_val, 6),
                "limma_overall_bh_fdr": round(fdr, 6),
                "kruskal_bh_fdr": round(fdr * random.uniform(0.5, 1.5), 6)
            })

    with open(os.path.join(os.path.dirname(__file__), "selection_gene_key.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=key_rows[0].keys())
        w.writeheader()
        w.writerows(key_rows)
    print(f"  Wrote delta_per_sample.csv ({len(rows)} rows) and selection_gene_key.csv")


if __name__ == "__main__":
    generate()
