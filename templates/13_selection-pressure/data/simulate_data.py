"""
Simulate selection pressure (piN-piS) data for template 13:
  - data/delta_per_sample.csv   — piN, piS and delta = piN - piS per gene/sample/DPI
  - data/selection_gene_key.csv — per-gene mean delta (descriptive key, no statistics)

Statistics are deliberately NOT computed here: template.qmd runs the per-gene
Kruskal-Wallis tests and the Benjamini-Hochberg correction on these per-sample
values, so the p-values on the published figure come from a real test.

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


GENES = ["nsp1", "nsp2", "nsp3", "nsp4", "capsid", "E3", "E2", "6K", "E1"]
DPI_VALUES = [1, 3, 5]
REPLICATES = 3

# Baseline synonymous nucleotide diversity (piS) shared by every gene.
PIS_MEAN = 0.012
PIS_SD = 0.003

# Per-gene excess of non-synonymous over synonymous diversity, i.e. the
# expected delta = piN - piS. Larger values = stronger positive selection.
GENE_DELTA_EFFECTS = {
    "nsp1": 0.001, "nsp2": 0.003, "nsp3": 0.008, "nsp4": 0.001,
    "capsid": 0.002, "E3": 0.012, "E2": 0.010, "6K": 0.015, "E1": 0.005
}
# gene -> (dpi at which the effect peaks, multiplier applied at that dpi)
GENE_EFFECT_DPI_BOOST = {"nsp3": (3, 2.0), "E2": (5, 2.5), "E3": (5, 1.8), "6K": (3, 1.5)}


def generate():
    rows = []
    for dpi_val in DPI_VALUES:
        for rep in range(1, REPLICATES + 1):
            sample_id = f"DPI{dpi_val}_R{rep}"
            for gene in GENES:
                base = GENE_DELTA_EFFECTS[gene]
                boost_dpi, factor = GENE_EFFECT_DPI_BOOST.get(gene, (None, 1.0))
                boost = factor if dpi_val == boost_dpi else 1.0

                # Draw piS and piN independently, then let delta fall out of them
                # so that delta == piN - piS exactly, including after clamping.
                pi_s = max(0.0, random.gauss(PIS_MEAN, PIS_SD))
                pi_n = max(0.0, pi_s + base * boost * random.gauss(1.0, 0.3))

                pi_n = round(pi_n, 6)
                pi_s = round(pi_s, 6)
                rows.append({
                    "threshold": "minfreq_0p01",
                    "sample": sample_id,
                    "dpi": f"dpi{dpi_val}",
                    "product": gene,
                    "piN": pi_n,
                    "piS": pi_s,
                    "delta_piN_minus_piS": round(pi_n - pi_s, 6),
                })

    with open(os.path.join(DATA_DIR, "delta_per_sample.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    # Gene-level key: descriptive means only. No p-values, no FDR — the template
    # derives those from the per-sample table above.
    key_rows = []
    for gene in GENES:
        vals = [r["delta_piN_minus_piS"] for r in rows if r["product"] == gene]
        if not vals:
            continue
        key_rows.append({
            "threshold": "minfreq_0p01",
            "product": gene,
            "n_samples": len(vals),
            "mean_delta": round(sum(vals) / len(vals), 6),
        })

    with open(os.path.join(DATA_DIR, "selection_gene_key.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(key_rows[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(key_rows)

    print(f"  Wrote delta_per_sample.csv ({len(rows)} rows) "
          f"and selection_gene_key.csv ({len(key_rows)} genes)")


if __name__ == "__main__":
    generate()
