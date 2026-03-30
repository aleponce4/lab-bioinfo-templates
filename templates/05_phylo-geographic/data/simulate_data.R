# simulate_data.R
# Generates synthetic data for template 05 (Phylo-Geographic):
#   - data/S_aligned.fasta   (S segment alignment, 30 sequences)
#   - data/M_aligned.fasta   (M segment alignment, 30 sequences)
#   - data/coordinates.csv   (sample locations)
#   - data/metadata.csv      (clade assignments)
#   - data/S_tree.nwk        (Newick tree)
# Run: Rscript data/simulate_data.R

set.seed(42)

library(ape)
library(Biostrings)

# ── Parameters ────────────────────────────────────────────────────────────────
# 30 samples across 8 sites arranged along a geographic gradient.
# Sites 1-4 host Clade_A/B (north); sites 5-8 host Clade_C/D (south).
# This creates a clear isolation-by-distance signal so the correlogram
# shows significant r at short distances and drops off at longer ones.

n_taxa    <- 30
taxa_names <- paste0("Sample_", seq_len(n_taxa))

# 8 sites along a north-south gradient, ~150 km total span
site_lookup <- data.frame(
  Site      = paste0("Site_", 1:8),
  Longitude = c(-69.5, -68.8, -68.2, -67.6, -67.0, -66.3, -65.7, -65.0),
  Latitude  = c(-22.2, -22.9, -23.5, -24.1, -24.8, -25.4, -26.1, -26.7)
)

# Assign samples to sites (3-4 per site) with clades matching geography
sample_plan <- data.frame(
  SampleID = taxa_names,
  Site = c(
    rep("Site_1", 4), rep("Site_2", 4), rep("Site_3", 3), rep("Site_4", 4),
    rep("Site_5", 4), rep("Site_6", 4), rep("Site_7", 3), rep("Site_8", 4)
  ),
  Clade = c(
    rep("Clade_A", 4), rep("Clade_A", 3), "Clade_B",        # sites 1-2
    rep("Clade_B", 3),                                        # site 3
    rep("Clade_B", 2), rep("Clade_A", 2),                    # site 4
    rep("Clade_C", 4), rep("Clade_C", 3), "Clade_D",        # sites 5-6
    rep("Clade_D", 3),                                        # site 7
    rep("Clade_D", 2), rep("Clade_C", 2)                     # site 8
  ),
  stringsAsFactors = FALSE
)

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

# ── Simulate tree with geographic signal ──────────────────────────────────────
# Build a tree where branch lengths reflect the clade/site structure:
# short branches within a site, longer between sites, longest between north/south.
# Achieved by making a balanced topology with scaled branch lengths.

# Intra-site branch length ~ 0.003, inter-site ~ 0.015, north/south ~ 0.08
make_site_clade <- function(ids, bl_intra = 0.003) {
  n <- length(ids)
  if (n == 1) return(list(tip = ids[1]))
  tr <- rtree(n, tip.label = ids, br = NULL)
  tr$edge.length <- runif(nrow(tr$edge), bl_intra * 0.5, bl_intra * 1.5)
  tr
}

# Build per-site subtrees then join
sites <- unique(sample_plan$Site)
subtrees <- lapply(sites, function(s) {
  ids <- sample_plan$SampleID[sample_plan$Site == s]
  make_site_clade(ids)
})

# Join sites 1-4 (north clade) and 5-8 (south clade), then join both
join_trees <- function(t1, t2, bl = 0.015) {
  tr <- bind.tree(t1, t2, where = "root")
  # scale the root edge
  tr$edge.length[tr$edge.length > 0.03] <- runif(
    sum(tr$edge.length > 0.03), 0.012, 0.02)
  tr
}

north <- Reduce(function(a, b) join_trees(a, b, 0.015), subtrees[1:4])
south <- Reduce(function(a, b) join_trees(a, b, 0.015), subtrees[5:8])

# Add a long inter-clade branch separating north and south
full_tree <- bind.tree(north, south, where = "root")
full_tree$edge.length[is.na(full_tree$edge.length)] <- 0.005
# Lengthen a few edges to represent the north-south split
long_edges <- which(full_tree$edge.length > 0.018)
if (length(long_edges) > 0)
  full_tree$edge.length[long_edges[1]] <- 0.08

write.tree(full_tree, file = "data/S_tree.nwk")
message("Wrote data/S_tree.nwk")

# ── Coordinates ───────────────────────────────────────────────────────────────
coords_joined <- merge(
  sample_plan[c("SampleID", "Site")],
  site_lookup,
  by = "Site",
  sort = FALSE
)
coords <- coords_joined[match(taxa_names, coords_joined$SampleID),
                        c("SampleID", "Longitude", "Latitude")]
write.csv(coords, "data/coordinates.csv", row.names = FALSE)

# ── Metadata ──────────────────────────────────────────────────────────────────
metadata <- data.frame(
  SampleID = sample_plan$SampleID,
  Clade    = sample_plan$Clade,
  Year     = sample(2015:2023, n_taxa, replace = TRUE),
  Host     = sample(c("Rodent_A", "Rodent_B"), n_taxa, replace = TRUE)
)
write.csv(metadata, "data/metadata.csv", row.names = FALSE)
message("Wrote coordinate and metadata files")
