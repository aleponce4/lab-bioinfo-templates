# simulate_data.R
# Generates synthetic data for template 05 (Phylo-Geographic):
#   - data/S_aligned.fasta   (S segment alignment, 10 sequences)
#   - data/M_aligned.fasta   (M segment alignment, 10 sequences)
#   - data/coordinates.csv   (sample locations)
#   - data/metadata.csv      (clade assignments)
#   - data/S_tree.nwk        (Newick tree from simulated topology)
# Run: Rscript data/simulate_data.R

set.seed(42)

library(ape)
library(Biostrings)

# ── Parameters ────────────────────────────────────────────────────────────────
n_taxa    <- 10
clades    <- c("Clade_A", "Clade_B", "Clade_C", "Clade_D")
clade_asgn <- sample(clades, n_taxa, replace = TRUE)
# Keep tip labels as plain "Sample_N" — the template regex extracts this exact pattern
taxa_names <- paste0("Sample_", seq_len(n_taxa))

# ── Simulate FASTA ─────────────────────────────────────────────────────────────
sim_fasta <- function(n, len, names) {
  seqs <- sapply(seq_len(n), function(i) {
    paste(sample(c("A","T","G","C"), len, replace = TRUE), collapse = "")
  })
  names(seqs) <- names
  DNAStringSet(seqs)
}

writeXStringSet(sim_fasta(n_taxa, 1200, taxa_names), "data/S_aligned.fasta")
writeXStringSet(sim_fasta(n_taxa, 3000, taxa_names), "data/M_aligned.fasta")
message("Wrote FASTA files")

# ── Simulate tree ─────────────────────────────────────────────────────────────
tree <- rtree(n_taxa, tip.label = taxa_names)
write.tree(tree, file = "data/S_tree.nwk")
message("Wrote data/S_tree.nwk")

# ── Simulate coordinates ───────────────────────────────────────────────────────
# Centered roughly around a fictional study area
coords <- data.frame(
  SampleID  = taxa_names,
  Longitude = runif(n_taxa, -70, -65),
  Latitude  = runif(n_taxa, -28, -22)
)
write.csv(coords, "data/coordinates.csv", row.names = FALSE)

# ── Simulate metadata ──────────────────────────────────────────────────────────
metadata <- data.frame(
  SampleID = taxa_names,
  Clade    = clade_asgn,
  Year     = sample(2015:2023, n_taxa, replace = TRUE),
  Host     = sample(c("Rodent_A", "Rodent_B"), n_taxa, replace = TRUE)
)
write.csv(metadata, "data/metadata.csv", row.names = FALSE)
message("Wrote coordinate and metadata files")
