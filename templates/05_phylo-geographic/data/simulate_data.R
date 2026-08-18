# simulate_data.R
# Generates synthetic data for template 05 (Phylo-Geographic):
#   - data/S_aligned.fasta   (S segment alignment, 30 sequences)
#   - data/M_aligned.fasta   (M segment alignment, 30 sequences)
#   - data/coordinates.csv   (sample locations)
#   - data/metadata.csv      (clade assignments)
#   - data/S_tree.nwk        (Newick tree)
#
# Design goal: the three artifacts must be MUTUALLY CONSISTENT.
#   geography  ->  tree topology & branch lengths  ->  sequence alignments
# so that the isolation-by-distance signal the template measures is really
# present in the data rather than being an artifact of random noise.
#
# Run: Rscript data/simulate_data.R

set.seed(42)

# Ensure working directory is set to template directory
if (!file.exists("template.qmd")) {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    setwd(dirname(script_dir))
  }
}
dir.create("data", showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(ape)
  library(Biostrings)
})

# ── Parameters ────────────────────────────────────────────────────────────────
# 30 samples across 8 sites arranged along a north-south geographic gradient.
# Sites 1-4 host Clade_A/B (north); sites 5-8 host Clade_C/D (south).

n_taxa     <- 30
taxa_names <- paste0("Sample_", seq_len(n_taxa))

TREE_DEPTH   <- 0.03   # root-to-tip divergence (substitutions/site) at site level
INTRA_SITE   <- c(0.0008, 0.0030)   # within-site terminal branch length range
SUBST_RATE   <- 1.0    # branch lengths are already in substitutions/site
S_LEN        <- 1200
M_LEN        <- 3000
AMBIG_RATE   <- 0.004  # fraction of positions masked as ambiguous, per sequence
N_LOWERCASE  <- 3      # sequences emitted in lowercase, as MAFFT does

# 8 sites in a realistic 2D scatter (not a straight line), ~600 km span
site_lookup <- data.frame(
  Site      = paste0("Site_", 1:8),
  Longitude = c(-69.8, -68.2, -70.1, -67.5, -68.9, -66.8, -69.3, -65.5),
  Latitude  = c(-21.5, -22.8, -24.0, -23.2, -25.6, -25.1, -27.0, -26.3),
  stringsAsFactors = FALSE
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

# ── 1. Site-level tree derived from the geography ─────────────────────────────
# Average-linkage clustering of the inter-site distances gives an ultrametric
# tree whose branching order IS the geographic structure. Because sites 1-4 and
# 5-8 form two spatially separated groups, the north/south split falls out
# automatically instead of having to be grafted on afterwards.
km_per_deg <- 111
site_xy <- cbind(
  x = site_lookup$Longitude * km_per_deg *
        cos(mean(site_lookup$Latitude) * pi / 180),
  y = site_lookup$Latitude * km_per_deg
)
rownames(site_xy) <- site_lookup$Site

site_tree <- as.phylo(hclust(dist(site_xy), method = "average"))
# Rescale so root-to-tip depth equals TREE_DEPTH substitutions/site.
site_tree$edge.length <- site_tree$edge.length *
  (TREE_DEPTH / max(node.depth.edgelength(site_tree)))

# Sanity check: the north (1-4) and south (5-8) sites must be reciprocally
# monophyletic, i.e. the clade split the comments promise actually exists.
north_sites <- paste0("Site_", 1:4)
south_sites <- paste0("Site_", 5:8)
stopifnot(is.monophyletic(site_tree, north_sites),
          is.monophyletic(site_tree, south_sites))

# ── 2. Expand each site into its sampled individuals ──────────────────────────
# bind.tree() with `where = <tip index>` NESTS the subtree at that tip. Using
# where = "root" instead would graft it beside the root and destroy the
# structure built above.
make_site_subtree <- function(ids) {
  bl <- runif(length(ids), INTRA_SITE[1], INTRA_SITE[2])
  if (length(ids) == 1) {
    tr <- read.tree(text = sprintf("(%s:%.6f);", ids, bl))
  } else {
    tr <- read.tree(text = sprintf(
      "(%s);", paste(sprintf("%s:%.6f", ids, bl), collapse = ",")))
  }
  tr$root.edge <- 0
  tr
}

full_tree <- site_tree
for (s in site_lookup$Site) {
  ids <- sample_plan$SampleID[sample_plan$Site == s]
  idx <- which(full_tree$tip.label == s)
  full_tree <- bind.tree(full_tree, make_site_subtree(ids), where = idx)
}

full_tree <- reorder(full_tree, "cladewise")
stopifnot(setequal(full_tree$tip.label, taxa_names))
stopifnot(all(full_tree$edge.length >= 0))

# Every sample from one site must still group together.
for (s in site_lookup$Site) {
  ids <- sample_plan$SampleID[sample_plan$Site == s]
  if (length(ids) > 1) stopifnot(is.monophyletic(full_tree, ids))
}

write.tree(full_tree, file = "data/S_tree.nwk")
message("Wrote data/S_tree.nwk  (", Ntip(full_tree), " tips)")

# ── 3. Simulate alignments ALONG that tree ────────────────────────────────────
# Jukes-Cantor substitution model. Sequences generated this way carry genuine
# phylogenetic signal: two tips that are close on the tree (and therefore close
# geographically) share more sites than two distant tips. Drawing independent
# random sequences instead — as an earlier version of this script did — makes
# the alignment unrelated to both the tree and the map.
BASES <- c("A", "C", "G", "T")

simulate_alignment <- function(tree, n_sites, rate = SUBST_RATE) {
  tree <- reorder(tree, "cladewise")   # parents always precede children
  n_node_total <- Ntip(tree) + Nnode(tree)
  root <- Ntip(tree) + 1L

  states <- vector("list", n_node_total)
  states[[root]] <- sample.int(4L, n_sites, replace = TRUE)

  for (e in seq_len(nrow(tree$edge))) {
    parent <- tree$edge[e, 1]
    child  <- tree$edge[e, 2]
    s      <- states[[parent]]
    bl     <- tree$edge.length[e]

    # P(site differs from parent) under Jukes-Cantor
    p_change <- 0.75 * (1 - exp(-4 / 3 * rate * bl))
    hits <- which(runif(n_sites) < p_change)
    if (length(hits) > 0) {
      # move to one of the three OTHER bases, uniformly
      s[hits] <- (s[hits] + sample.int(3L, length(hits), replace = TRUE) - 1L) %% 4L + 1L
    }
    states[[child]] <- s
  }

  mat <- do.call(rbind, states[seq_len(Ntip(tree))])
  out <- apply(mat, 1, function(r) paste(BASES[r], collapse = ""))
  names(out) <- tree$tip.label
  out
}

# Inject ambiguous bases so the template's alignment-QC table is meaningful.
# Real alignments carry N plus the other IUPAC degeneracy codes, and MAFFT
# writes its output in lowercase — both are exercised here.
IUPAC_AMBIG <- c("N", "N", "N", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V")

add_ambiguity <- function(seqs, rate = AMBIG_RATE, n_lower = N_LOWERCASE) {
  ids <- names(seqs)
  out <- vapply(seq_along(seqs), function(i) {
    chars <- strsplit(seqs[[i]], "", fixed = TRUE)[[1]]
    # sequence-specific quality: a few sequences are much worse than the rest
    r <- rate * ifelse(i %% 7 == 0, 6, ifelse(i %% 3 == 0, 1, 0.15))
    k <- rbinom(1, length(chars), min(r, 1))
    if (k > 0) {
      pos <- sample.int(length(chars), k)
      chars[pos] <- sample(IUPAC_AMBIG, k, replace = TRUE)
    }
    paste(chars, collapse = "")
  }, character(1))
  names(out) <- ids

  # MAFFT-style lowercase output for a subset of sequences
  if (n_lower > 0) {
    lower_idx <- sample.int(length(out), min(n_lower, length(out)))
    out[lower_idx] <- tolower(out[lower_idx])
  }
  out
}

seqs_S <- add_ambiguity(simulate_alignment(full_tree, S_LEN))
seqs_M <- add_ambiguity(simulate_alignment(full_tree, M_LEN))

writeXStringSet(DNAStringSet(seqs_S), "data/S_aligned.fasta")
writeXStringSet(DNAStringSet(seqs_M), "data/M_aligned.fasta")

count_ambig <- function(x) sum(vapply(x, function(s)
  sum(!strsplit(toupper(s), "", fixed = TRUE)[[1]] %in% BASES), integer(1)))
message("Wrote FASTA files  (S: ", S_LEN, " bp, M: ", M_LEN, " bp; ",
        count_ambig(seqs_S), " + ", count_ambig(seqs_M), " ambiguous bases)")
stopifnot(count_ambig(seqs_S) > 0, count_ambig(seqs_M) > 0)

# ── 4. Coordinates ────────────────────────────────────────────────────────────
coords_joined <- merge(
  sample_plan[c("SampleID", "Site")],
  site_lookup,
  by = "Site",
  sort = FALSE
)
coords <- coords_joined[match(taxa_names, coords_joined$SampleID),
                        c("SampleID", "Longitude", "Latitude")]
write.csv(coords, "data/coordinates.csv", row.names = FALSE)

# ── 5. Metadata ───────────────────────────────────────────────────────────────
metadata <- data.frame(
  SampleID = sample_plan$SampleID,
  Clade    = sample_plan$Clade,
  Year     = sample(2015:2023, n_taxa, replace = TRUE),
  Host     = sample(c("Rodent_A", "Rodent_B"), n_taxa, replace = TRUE)
)
write.csv(metadata, "data/metadata.csv", row.names = FALSE)
message("Wrote coordinate and metadata files")

# ── 6. Report the signal that was actually built in ──────────────────────────
# A quick Mantel-style check so a regression in this script is visible here
# rather than only downstream in the template.
pat <- cophenetic.phylo(full_tree)[taxa_names, taxa_names]
geo_xy <- cbind(
  x = coords$Longitude * km_per_deg * cos(mean(coords$Latitude) * pi / 180),
  y = coords$Latitude * km_per_deg
)
geo <- as.matrix(dist(geo_xy))
r <- cor(pat[lower.tri(pat)], geo[lower.tri(geo)])
message(sprintf("Built-in isolation-by-distance signal: Pearson r = %.3f", r))
stopifnot(r > 0.5)
