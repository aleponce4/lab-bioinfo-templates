"""
simulate_data.py
Generates synthetic VCF files for Template 08 (VCF Mutation Analysis).

Output: data/aa_change_vcf/*.vcf
VCF filename format: aa_change_{ID}-{SOURCE}_{SEGMENT}_{THRESHOLD}.vcf

Run from the template folder:
    python data/simulate_data.py

Three samples × three tissue types × two segments × three thresholds = 54 VCF files.
Synthetic mutations are randomly seeded for reproducibility.
"""

import os
import random
import string
import math

random.seed(42)

# ── Output directory ──────────────────────────────────────────────────────────
OUT_DIR = os.path.join(os.path.dirname(__file__), "aa_change_vcf")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Configuration mirrors template config ─────────────────────────────────────
SAMPLE_IDS = ["SMPL001", "SMPL002", "SMPL003",
               "SMPL004", "SMPL005", "SMPL006"]
TISSUE_TYPES = ["Lung", "Saliva", "Urine"]
SEGMENTS = ["S", "M"]
VAF_THRESHOLDS = ["1%", "2%", "5%"]

# Segment lengths (nt) and CDS boundaries
SEGMENT_CONFIG = {
    "S": {"length": 1902, "cds_start": 250, "cds_end": 1329, "aa_length": 360},
    "M": {"length": 3675, "cds_start":  52, "cds_end": 3461, "aa_length": 1135},
}

DEPTH_MIN = 500   # minimum coverage to pass QC filter

NUCLEOTIDES = ["A", "T", "G", "C"]
AA_CHANGES = [
    "Lys4Arg", "Glu17Gly", "Val52Ile", "Thr91Ala", "Met120Val",
    "Ser200Pro", "Leu245Phe", "Asn310Asp", "Arg412Ser", "His500Tyr",
    "N/A",  # synonymous — included ~40% of the time
]

# ── VCF template ──────────────────────────────────────────────────────────────
VCF_HEADER = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=AA_CHANGE,Number=1,Type=String,Description="Amino acid change">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
"""

CHROM = "Segment"


def random_nucleotide(exclude=""):
    choices = [n for n in NUCLEOTIDES if n != exclude]
    return random.choice(choices)


def make_vcf_records(segment, vaf_str):
    """Generate a plausible set of variant records for one VCF file."""
    cfg = SEGMENT_CONFIG[segment]
    seg_len = cfg["length"]
    cds_start = cfg["cds_start"]
    cds_end = cfg["cds_end"]

    # Parse numeric VAF floor
    vaf_floor = float(vaf_str.rstrip("%")) / 100.0

    # Number of variants varies by threshold (more at lower thresholds)
    n_variants = {0.05: random.randint(8, 18),
                  0.02: random.randint(12, 25),
                  0.01: random.randint(18, 35)}[vaf_floor]

    records = []
    positions_used = set()
    for _ in range(n_variants):
        # Pick a unique position
        for _attempt in range(100):
            pos = random.randint(1, seg_len)
            if pos not in positions_used:
                positions_used.add(pos)
                break

        ref = random_nucleotide()
        alt = random_nucleotide(exclude=ref)

        depth = random.randint(DEPTH_MIN, DEPTH_MIN * 6)
        # Allele frequency bounded between vaf_floor and vaf_floor * 4 (capped at 0.5)
        vaf = min(0.5, random.uniform(vaf_floor, vaf_floor * 4))
        alt_depth = max(1, int(round(depth * vaf)))
        ref_depth = depth - alt_depth

        qual = round(random.uniform(30, 60), 1)

        # Amino acid change — only in CDS
        if cds_start <= pos <= cds_end:
            aa_change = random.choices(
                AA_CHANGES,
                weights=[1] * (len(AA_CHANGES) - 1) + [2],  # slight bias toward N/A
            )[0]
        else:
            aa_change = "N/A"

        info = f"AA_CHANGE={aa_change}"
        fmt = "GT:AD:DP"
        sample = f"0/1:{ref_depth},{alt_depth}:{depth}"

        records.append(
            f"{CHROM}\t{pos}\t.\t{ref}\t{alt}\t{qual}\tPASS\t{info}\t{fmt}\t{sample}\n"
        )

    # Sort by position
    records.sort(key=lambda r: int(r.split("\t")[1]))
    return records


# ── Write VCF files ───────────────────────────────────────────────────────────
file_count = 0
for sample_id in SAMPLE_IDS:
    for source in TISSUE_TYPES:
        for segment in SEGMENTS:
            for threshold in VAF_THRESHOLDS:
                threshold_str = threshold.replace("%", "pct")
                fname = f"aa_change_{sample_id}-{source}_{segment}_{threshold_str}.vcf"
                fpath = os.path.join(OUT_DIR, fname)

                with open(fpath, "w") as fh:
                    fh.write(VCF_HEADER)
                    for rec in make_vcf_records(segment, threshold):
                        fh.write(rec)

                file_count += 1

print(f"Wrote {file_count} synthetic VCF files to {OUT_DIR}/")
print("Each file follows the naming convention:")
print("  aa_change_{SampleID}-{Source}_{Segment}_{Threshold}.vcf")
print("Example: aa_change_SMPL001-Lung_S_1pct.vcf")
