"""
Simulate intra-host variant data for template 11:
  - data/gene_coords.csv      — genome annotation (gene start/end)
  - data/LoFreq_variants.csv  — LoFreq VCF-derived variant calls
  - data/iVar_variants.csv    — iVar VCF-derived variant calls
  - data/manifest.csv         — sample x DPI x tissue metadata

Both callers are derived from ONE shared truth set of variant positions per sample,
so the caller-overlap Venn diagram in the template compares two overlapping call
sets the way two callers run on the same reads actually would. Drawing the two
position sets independently produces ~1% concordance, which is not a realistic
input for that figure.

Run from the template folder:
    python data/simulate_data.py
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
NUCLEOTIDES = ["A", "C", "G", "T"]

DPI_MEAN_VARIANTS = {1: 45, 3: 120, 5: 180}
REPLICATES = 3
TISSUE = "Lung"
FREQ_LO = 0.01
FREQ_HI = 0.95
DEPTH_MEAN = 8000
DEPTH_SD = 2000

# How each caller derives its call set from the shared per-sample truth set.
#   dropout_rate      — fraction of true variants the caller misses
#   extra_calls       — caller-specific calls not in the truth set
#   af_jitter_sigma   — lognormal spread applied to the true allele frequency
CALLER_BEHAVIOUR = {
    "LoFreq": {"dropout_rate": 0.05, "extra_calls": 5, "af_jitter_sigma": 0.08},
    "iVar":   {"dropout_rate": 0.20, "extra_calls": 0, "af_jitter_sigma": 0.12},
}
CALLERS = list(CALLER_BEHAVIOUR)


def gene_at(pos):
    for name, (s, e) in GENES.items():
        if s <= pos <= e:
            return name
    return "intergenic"


def consequence(pos, gene):
    """Assign intergenic, missense, or synonymous based on gene and codon position."""
    if gene == "intergenic":
        return "intergenic_variant"
    return "missense_variant" if (pos % 3) != 0 else "synonymous_variant"


def amino_acid_change(pos, cons):
    if cons != "missense_variant":
        return "synonymous" if cons == "synonymous_variant" else "none"
    codon_start = ((pos - 1) // 3) * 3 + 1
    return f"p.Val{codon_start // 3 + 1}Ile" if random.random() < 0.5 else f"p.Leu{codon_start // 3 + 1}Phe"


def make_truth_variant(sample_id, dpi_val, rep, pos):
    """One true intra-host variant, independent of which caller finds it."""
    ref = random.choice(NUCLEOTIDES)
    alt = random.choice([n for n in NUCLEOTIDES if n != ref])
    af = 10 ** random.uniform(math.log10(FREQ_LO), math.log10(FREQ_HI))
    gene = gene_at(pos)
    cons = consequence(pos, gene)
    return {
        "sample_id": sample_id,
        "dpi": dpi_val,
        "replicate": rep,
        "tissue": TISSUE,
        "CHROM": "PathogenX",
        "POS": pos,
        "REF": ref,
        "ALT": alt,
        "AF": af,
        "consequence": cons,
        "gene": gene,
        "amino_acid_change": amino_acid_change(pos, cons),
        "DNA_change": f"c.{pos}{ref}>{alt}",
    }


def as_call(truth, caller, af_jitter_sigma):
    """One caller's observation of a variant: same locus, slightly different AF/DP."""
    af = truth["AF"] * random.lognormvariate(0.0, af_jitter_sigma)
    af = min(FREQ_HI, max(FREQ_LO, af))
    return {
        "sample_id": truth["sample_id"],
        "caller": caller,
        "dpi": truth["dpi"],
        "replicate": truth["replicate"],
        "tissue": truth["tissue"],
        "CHROM": truth["CHROM"],
        "POS": truth["POS"],
        "REF": truth["REF"],
        "ALT": truth["ALT"],
        "AF": round(af, 6),
        "DP": max(500, int(random.gauss(DEPTH_MEAN, DEPTH_SD))),
        "consequence": truth["consequence"],
        "gene": truth["gene"],
        "amino_acid_change": truth["amino_acid_change"],
        "DNA_change": truth["DNA_change"],
    }


def generate():
    calls = {caller: [] for caller in CALLERS}
    n_truth_total = 0

    for dpi_val, n_variants in DPI_MEAN_VARIANTS.items():
        for rep in range(1, REPLICATES + 1):
            sample_id = f"DPI{dpi_val}_R{rep}_{TISSUE}"

            # ── One shared truth set of positions for this sample ──────────────
            n = max(10, min(int(random.gauss(n_variants, 15)), 300))
            truth_positions = random.sample(range(1, GENOME_LEN + 1), n)
            truth = [make_truth_variant(sample_id, dpi_val, rep, pos)
                     for pos in sorted(truth_positions)]
            n_truth_total += len(truth)

            # ── Each caller sees a subset of it, plus its own extra calls ──────
            for caller, behaviour in CALLER_BEHAVIOUR.items():
                detected = [t for t in truth
                            if random.random() >= behaviour["dropout_rate"]]
                sample_calls = [as_call(t, caller, behaviour["af_jitter_sigma"])
                                for t in detected]

                # Caller-specific calls at positions no other caller reports.
                taken = set(truth_positions)
                for _ in range(behaviour["extra_calls"]):
                    pos = random.randint(1, GENOME_LEN)
                    if pos in taken:
                        continue
                    taken.add(pos)
                    extra = make_truth_variant(sample_id, dpi_val, rep, pos)
                    sample_calls.append(
                        as_call(extra, caller, behaviour["af_jitter_sigma"])
                    )

                sample_calls.sort(key=lambda c: c["POS"])
                calls[caller].extend(sample_calls)

    # Write gene coordinates
    with open(os.path.join(DATA_DIR, "gene_coords.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "start", "end"])
        for name, (s, e) in GENES.items():
            w.writerow([name, s, e])

    # Write manifest
    with open(os.path.join(DATA_DIR, "manifest.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample_id", "dpi", "replicate", "tissue"])
        for dpi_val in DPI_MEAN_VARIANTS:
            for rep in range(1, REPLICATES + 1):
                w.writerow([f"DPI{dpi_val}_R{rep}_{TISSUE}", dpi_val, rep, TISSUE])

    # Write per-caller variant files
    for caller in CALLERS:
        rows = calls[caller]
        fname = f"{caller}_variants.csv"
        with open(os.path.join(DATA_DIR, fname), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        print(f"  Wrote {fname} ({len(rows)} variants)")

    print(f"Generated {n_truth_total} true variants across "
          f"{len(DPI_MEAN_VARIANTS)} DPIs x {REPLICATES} reps; "
          f"call sets derived per caller.")


if __name__ == "__main__":
    generate()
