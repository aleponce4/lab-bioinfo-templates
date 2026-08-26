"""
Simulate genome coverage depth data for template 12:
  - data/gene_coords.csv     — genome annotation
  - data/coverage.csv        — per-position depth per sample
  - data/manifest.csv        — sample metadata
"""

import csv
import math
import os
import random

SEED = 42
random.seed(SEED)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TEMPLATE_DIR = os.path.dirname(SCRIPT_DIR)
DATA_DIR = os.path.join(TEMPLATE_DIR, "data")
os.makedirs(DATA_DIR, exist_ok=True)

GENOME_LEN = 11447
GENES = {
    "nsp1": (1, 1800), "nsp2": (1801, 4500), "nsp3": (4501, 5700),
    "nsp4": (5701, 7500), "capsid": (7501, 8300), "E3": (8301, 9000),
    "E2": (9001, 10200), "6K": (10201, 10400), "E1": (10401, 11400)
}

DPI_REPS = {1: 3, 3: 3, 5: 3}
# Invented base depth scales per timepoint — placeholder numbers, not measurements.
# They rise with DPI because viral load rises, so a fixed sequencing effort yields more
# on-target reads later. Kept within one order of magnitude: real timepoints differ a few
# fold, not a few hundred fold.
DPI_BASE_DEPTH = {1: 8000, 3: 22000, 5: 45000}
# Coverage is autocorrelated along the genome (fragments span many positions), so depth is
# drawn as a smooth multiplicative field plus per-position noise rather than independent
# draws at every base.
FIELD_SCALE = 250.0   # nt; correlation length of the smooth component
FIELD_SD = 0.35       # log-space amplitude of the smooth component
POS_NOISE_SD = 0.12   # log-space per-position jitter
AMPLICON_DROPS = [(2000, 2200), (4300, 4450), (7000, 7150)]
# Width (nt) of the coverage taper at each genome terminus.
TERMINAL_TAPER = 30.0


def generate():
    rows = []
    for dpi_val, n_reps in DPI_REPS.items():
        base_depth = DPI_BASE_DEPTH[dpi_val]
        for rep in range(1, n_reps + 1):
            sample_id = f"DPI{dpi_val}_R{rep}_Lung"

            # Smooth multiplicative field: a random value every FIELD_SCALE nt, linearly
            # interpolated between anchors. This gives coverage its characteristic gentle
            # undulation instead of independent per-base noise.
            n_anchors = int(GENOME_LEN / FIELD_SCALE) + 2
            anchors = [random.gauss(0.0, FIELD_SD) for _ in range(n_anchors)]

            # Simulate genome coverage with synthetic terminal drop-offs, subgenomic elevation, amplicon dropouts
            for pos in range(1, GENOME_LEN + 1):
                # 3' structural gene subgenomic mRNA elevation (26S RNA)
                subgenomic_mult = 1.8 if pos > 7500 else 1.0

                a_idx = (pos - 1) / FIELD_SCALE
                lo = int(a_idx)
                frac = a_idx - lo
                field = anchors[lo] * (1.0 - frac) + anchors[lo + 1] * frac

                depth = (base_depth * subgenomic_mult
                         * math.exp(field + random.gauss(0.0, POS_NOISE_SD)))
                
                # Amplicon dropout regions
                for ds, de in AMPLICON_DROPS:
                    if ds <= pos <= de:
                        depth *= random.uniform(0.08, 0.35)
                        
                # Terminal drop-offs at genome 5' and 3' ends
                if pos <= TERMINAL_TAPER:
                    depth *= (pos / TERMINAL_TAPER) ** 1.5
                elif pos > GENOME_LEN - TERMINAL_TAPER:
                    depth *= ((GENOME_LEN - pos + 1) / TERMINAL_TAPER) ** 1.5
                    
                depth = max(0, int(round(depth)))
                rows.append({
                    "CHROM": "PathogenX",
                    "Position": pos,
                    "Depth": depth,
                    "sample_id": sample_id,
                    "dpi": dpi_val,
                    "replicate": rep
                })

    with open(os.path.join(DATA_DIR, "coverage.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys(), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"  Wrote data/coverage.csv ({len(rows)} rows)")

    with open(os.path.join(DATA_DIR, "gene_coords.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "start", "end"])
        for name, (s, e) in GENES.items():
            w.writerow([name, s, e])

    with open(os.path.join(DATA_DIR, "manifest.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample_id", "dpi", "replicate"])
        for dpi_val, n_reps in DPI_REPS.items():
            for rep in range(1, n_reps + 1):
                w.writerow([f"DPI{dpi_val}_R{rep}_Lung", dpi_val, rep])


if __name__ == "__main__":
    generate()
