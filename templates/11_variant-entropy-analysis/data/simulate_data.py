"""
Simulate intra-host variant data for template 09:
  - data/gene_coords.csv      — genome annotation (gene start/end)
  - data/loFreq_variants.csv  — LoFreq VCF-derived variant calls
  - data/iVar_variants.csv    — iVar VCF-derived variant calls
  - data/manifest.csv         — sample × DPI × tissue metadata
"""

import csv, random, math, os

SEED = 42
random.seed(SEED)

DATA_DIR = os.path.join(os.path.dirname(__file__) if "__file__" in dir() else ".", "..")
os.makedirs(os.path.dirname(os.path.abspath(__file__)), exist_ok=True)

GENOME_LEN = 11447
GENES = {
    "nsp1": (1, 1800), "nsp2": (1801, 4500), "nsp3": (4501, 5700),
    "nsp4": (5701, 7500), "capsid": (7501, 8300), "E3": (8301, 9000),
    "E2": (9001, 10200), "6K": (10201, 10400), "E1": (10401, 11400)
}
NUCLEOTIDES = ["A", "C", "G", "T"]

DPI_MEAN_VARIANTS = {1: 45, 3: 120, 5: 180}
REPLICATES = 3
CALLERS = ["LoFreq", "iVar"]
FREQ_LO = 0.01
FREQ_HI = 0.95
DEPTH_MEAN = 8000
DEPTH_SD = 2000


def gene_at(pos):
    for name, (s, e) in GENES.items():
        if s <= pos <= e:
            return name
    return "intergenic"


def consequence(pos, gene):
    """Assign missense, synonymous, or intergenic based on position within codon."""
    return "missense_variant" if (pos % 3) != 0 else "synonymous_variant"


def amino_acid_change(ref, alt, pos, is_missense):
    if not is_missense:
        return "synonymous"
    codon_start = ((pos - 1) // 3) * 3 + 1
    return f"p.Val{codon_start // 3 + 1}Ile" if random.random() < 0.5 else f"p.Leu{codon_start // 3 + 1}Phe"


def generate():
    variants = []
    for dpi_val, n_variants in DPI_MEAN_VARIANTS.items():
        for rep in range(1, REPLICATES + 1):
            for caller in CALLERS:
                sample_id = f"DPI{dpi_val}_R{rep}_Lung"
                n = int(random.gauss(n_variants, 15))
                n = max(10, min(n, 300))

                used_positions = set()
                for _ in range(n):
                    while True:
                        pos = random.randint(1, GENOME_LEN)
                        if pos not in used_positions:
                            used_positions.add(pos)
                            break

                    ref = random.choice(NUCLEOTIDES)
                    alt = random.choice([n for n in NUCLEOTIDES if n != ref])
                    af = 10 ** random.uniform(math.log10(FREQ_LO), math.log10(FREQ_HI))
                    dp = max(500, int(random.gauss(DEPTH_MEAN, DEPTH_SD)))
                    gene = gene_at(pos)
                    cons = consequence(pos, gene)
                    is_missense = "missense" in cons
                    aa = amino_acid_change(ref, alt, pos, is_missense)

                    variants.append({
                        "sample_id": sample_id,
                        "caller": caller,
                        "dpi": dpi_val,
                        "replicate": rep,
                        "tissue": "Lung",
                        "CHROM": "VEEV",
                        "POS": pos,
                        "REF": ref,
                        "ALT": alt,
                        "AF": round(af, 6),
                        "DP": dp,
                        "consequence": cons,
                        "gene": gene,
                        "amino_acid_change": aa,
                        "DNA_change": f"c.{pos}{ref}>{alt}"
                    })

    # Write gene coordinates
    with open(os.path.join(os.path.dirname(__file__), "gene_coords.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "start", "end"])
        for name, (s, e) in GENES.items():
            w.writerow([name, s, e])

    # Write manifest
    with open(os.path.join(os.path.dirname(__file__), "manifest.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample_id", "dpi", "replicate", "tissue"])
        for dpi_val in DPI_MEAN_VARIANTS:
            for rep in range(1, REPLICATES + 1):
                w.writerow([f"DPI{dpi_val}_R{rep}_Lung", dpi_val, rep, "Lung"])

    # Write per-caller variant files
    for caller in CALLERS:
        caller_vars = [v for v in variants if v["caller"] == caller]
        # Remove ~20% of iVar variants to simulate caller differences
        if caller == "iVar":
            caller_vars = [v for v in caller_vars if random.random() > 0.20]
        # Add 5 unique LoFreq positions
        if caller == "LoFreq":
            for _ in range(5):
                p = random.randint(1, GENOME_LEN)
                af = 10 ** random.uniform(math.log10(FREQ_LO), math.log10(FREQ_HI))
                g = gene_at(p)
                caller_vars.append({
                    "sample_id": f"DPI{random.choice(list(DPI_MEAN_VARIANTS.keys()))}_R{random.randint(1,3)}_Lung",
                    "caller": "LoFreq", "dpi": 3, "replicate": 1, "tissue": "Lung",
                    "CHROM": "VEEV", "POS": p, "REF": "A", "ALT": "G", "AF": round(af, 6),
                    "DP": 5000, "consequence": "missense_variant", "gene": g,
                    "amino_acid_change": "p.Thr100Ile", "DNA_change": f"c.{p}A>G"
                })

        fname = f"{'LoFreq' if caller == 'LoFreq' else 'iVar'}_variants.csv"
        with open(os.path.join(os.path.dirname(__file__), fname), "w", newline="") as f:
            if caller_vars:
                w = csv.DictWriter(f, fieldnames=caller_vars[0].keys())
                w.writeheader()
                w.writerows(caller_vars)
        print(f"  Wrote {fname} ({len(caller_vars)} variants)")

    print(f"Generated {len(variants)} total variants across {len(DPI_MEAN_VARIANTS)} DPIs × {REPLICATES} reps × {len(CALLERS)} callers")


if __name__ == "__main__":
    generate()
