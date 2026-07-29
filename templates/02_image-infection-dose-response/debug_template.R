suppressPackageStartupMessages({
  library(tidyverse)
  library(mclust)
  library(minpack.lm)
})

source("data/simulate_data.R")

NUCLEUS_FILE <- "data/nucleus.csv"
VIRUS_FILE   <- "data/virus.csv"
NUCLEUS_DAPI_COL    <- "(nucleus) MeanIntensity CH1"
NUCLEUS_AREA_COL    <- "(nucleus) Area"
VIRUS_INTENSITY_COL <- "(Virus) MeanIntensity CH2"

VIRUS_NAME    <- "Virus"
COMPOUND_NAME <- "Compound"
DOSE_UNIT     <- "\u00b5M"

AREA_MIN        <- 30
AREA_MAX        <- 400
INTENSITY_MIN   <- 300
CIRCULARITY_MIN <- 0.40

DOSES_TO_EXCLUDE <- integer(0)
WELLS_TO_EXCLUDE <- character(0)

PLATE_MAP <- tibble::tribble(
  ~WellName,  ~Treatment,          ~Concentration_uM,
  "C3",  "Mock",             NA,
  "C4",  "Mock",             NA,
  "C5",  "Mock",             NA,
  "C6",  "Mock",             NA,
  "C11", "Mock+Ab",          NA,
  "C7",  "Virus Only",        0,
  "C8",  "Virus Only",        0,
  "D1",  "Toxicity Control",  0.01,
  "D2",  "Toxicity Control",  0.1,
  "D3",  "Toxicity Control",  0.5,
  "D4",  "Toxicity Control",  1,
  "D5",  "Toxicity Control",  2,
  "D6",  "Toxicity Control",  4,
  "D7",  "Toxicity Control",  8,
  "D8",  "Toxicity Control",  12,
  "D9",  "Toxicity Control",  16,
  "D10", "Toxicity Control",  20,
  "D11", "Toxicity Control",  30,
  "D12", "Toxicity Control",  50,
  "E1",  "Virus + Compound",  0.01,
  "E2",  "Virus + Compound",  0.1,
  "E3",  "Virus + Compound",  0.5,
  "E4",  "Virus + Compound",  1,
  "E5",  "Virus + Compound",  2,
  "E6",  "Virus + Compound",  4,
  "E7",  "Virus + Compound",  8,
  "E8",  "Virus + Compound",  12,
  "E9",  "Virus + Compound",  16,
  "E10", "Virus + Compound",  20,
  "E11", "Virus + Compound",  30,
  "E12", "Virus + Compound",  50,
  "F1",  "Virus + Compound",  0.01,
  "F2",  "Virus + Compound",  0.1,
  "F3",  "Virus + Compound",  0.5,
  "F4",  "Virus + Compound",  1,
  "F5",  "Virus + Compound",  2,
  "F6",  "Virus + Compound",  4,
  "F7",  "Virus + Compound",  8,
  "F8",  "Virus + Compound",  12,
  "F9",  "Virus + Compound",  16,
  "F10", "Virus + Compound",  20,
  "F11", "Virus + Compound",  30,
  "F12", "Virus + Compound",  50
)

df_nuclei <- read_csv(NUCLEUS_FILE, show_col_types = FALSE)
df_virus  <- read_csv(VIRUS_FILE,   show_col_types = FALSE)

nuclei_clean <- df_nuclei %>% dplyr::select(WellName, FieldIndex, ObjectNumber, dapi_intensity = all_of(NUCLEUS_DAPI_COL), Area = all_of(NUCLEUS_AREA_COL)) %>% mutate(Sphericity = 1.0)
virus_clean  <- df_virus  %>% dplyr::select(WellName, FieldIndex, ObjectNumber, virus_intensity = all_of(VIRUS_INTENSITY_COL))

master_df <- nuclei_clean %>% left_join(virus_clean, by = c("WellName", "FieldIndex", "ObjectNumber")) %>%
  mutate(WellName_std = str_replace(str_to_upper(WellName), "^([A-H])-?0*([1-9][0-9]?).*$", "\\1\\2"))

master_df_annotated <- master_df %>% left_join(PLATE_MAP, by = c("WellName_std" = "WellName")) %>% filter(!is.na(Treatment))

master_df_gated <- master_df_annotated %>% filter(Area > AREA_MIN & Area < AREA_MAX & dapi_intensity > INTENSITY_MIN) %>% mutate(Plate = "P1")

mock_cells <- master_df_gated %>% filter(Treatment %in% c("Mock", "Mock+Ab")) %>% mutate(xlog = log10(virus_intensity + 1))
L_p_val    <- median(mock_cells$xlog, na.rm = TRUE)

vo_cells   <- master_df_gated %>% filter(Treatment == "Virus Only") %>% mutate(xlog = log10(virus_intensity + 1))
H_p_val    <- quantile(vo_cells$xlog, 0.75, na.rm = TRUE)

master_df_norm <- master_df_gated %>% mutate(xlog = log10(virus_intensity + 1), L_p = L_p_val, H_p = H_p_val, span = H_p - L_p, x_norm = (xlog - L_p)/span)

threshold_final <- 0.5
total_counts <- master_df_norm %>% dplyr::count(Treatment, Concentration_uM, WellName_std, name = "total_cell_count")
positive_counts <- master_df_norm %>% filter(x_norm > threshold_final) %>% dplyr::count(Treatment, Concentration_uM, WellName_std, name = "positive_cell_count")

percent_summary <- total_counts %>% left_join(positive_counts, by = c("Treatment", "Concentration_uM", "WellName_std")) %>%
  mutate(positive_cell_count = replace_na(positive_cell_count, 0), percent_positive = (positive_cell_count/total_cell_count)*100)

vo_anchor   <- mean(percent_summary$percent_positive[percent_summary$Treatment == "Virus Only"], na.rm = TRUE)
mock_anchor <- mean(percent_summary$percent_positive[percent_summary$Treatment %in% c("Mock", "Mock+Ab")], na.rm = TRUE)

summary_rel <- percent_summary %>% mutate(percent_infection_rel = 100 * (percent_positive - mock_anchor)/(vo_anchor - mock_anchor))

drc_data_rel <- summary_rel %>% filter(Treatment == "Virus + Compound")
pos_conc  <- drc_data_rel$Concentration_uM[drc_data_rel$Concentration_uM > 0 & is.finite(drc_data_rel$Concentration_uM)]
min_pos   <- if (length(pos_conc) > 0) min(pos_conc) else 0.01
zero_tick <- min_pos / 10

fit_data_wide <- summary_rel %>%
  filter(Treatment == "Virus Only" | Treatment == "Virus + Compound") %>%
  mutate(
    dose_fit   = ifelse(Treatment == "Virus Only", zero_tick, Concentration_uM),
    prop_raw   = positive_cell_count / total_cell_count,
    prop_clamp = pmin(pmax(prop_raw, 1e-3), 1 - 1e-3),
    bin_weight = total_cell_count / (prop_clamp * (1 - prop_clamp))
  )

fit_rel <- tryCatch(
  minpack.lm::nlsLM(
    percent_infection_rel ~ bottom + (top - bottom) / (1 + (dose_fit / ec50) ^ hill),
    data    = fit_data_wide,
    start   = list(ec50 = 8, hill = 0.8, bottom = 0, top = 100),
    lower   = c(ec50 = 1e-4, hill = 0.05, bottom = -20, top = 80),
    upper   = c(ec50 = 1e3,  hill = 50, bottom = 50, top = 120),
    weights = fit_data_wide$bin_weight,
    control = minpack.lm::nls.lm.control(maxiter = 500)
  ),
  error   = function(e) NULL
)

cat("=== FIT_REL COEFFICIENTS ===\n")
print(if (!is.null(fit_rel)) coef(fit_rel) else "FIT FAILED")

pts_sum <- summary_rel %>%
  filter(Treatment %in% c("Virus Only", "Virus + Compound")) %>%
  mutate(x_plot = ifelse(Treatment == "Virus Only", zero_tick, Concentration_uM)) %>%
  group_by(x_plot) %>%
  summarise(n  = n(),
            y  = mean(percent_infection_rel, na.rm = TRUE),
            se = ifelse(n > 1, sd(percent_infection_rel, na.rm = TRUE) / sqrt(n), 0),
            .groups = "drop") %>%
  mutate(y    = ifelse(x_plot == zero_tick, 100, y),
         ymin = pmax(0,   y - se),
         ymax = pmin(100, y + se))

cat("=== PTS_SUM CONTENT ===\n")
print(pts_sum)

max_pos <- if (length(pos_conc) > 0) max(pos_conc) else 100
newx    <- exp(seq(log(zero_tick), log(max_pos), length.out = 300))

pred_df_rel <- tibble(
  Concentration_uM = newx,
  y = coef(fit_rel)["bottom"] + (coef(fit_rel)["top"] - coef(fit_rel)["bottom"]) / (1 + (Concentration_uM / coef(fit_rel)["ec50"])^coef(fit_rel)["hill"])
)

dose_breaks     <- sort(unique(c(zero_tick, drc_data_rel$Concentration_uM)))
dose_labels_vec <- c("0", as.character(dose_breaks[-1]))

p_dr_rel <- ggplot() +
  geom_line(data = pred_df_rel, aes(x = Concentration_uM, y = y), color = "steelblue", linewidth = 0.9) +
  geom_errorbar(data = pts_sum, aes(x = x_plot, ymin = ymin, ymax = ymax), width = 0.06, linewidth = 0.5, colour = "steelblue4") +
  geom_point(data = pts_sum, aes(x = x_plot, y = y), shape = 21, size = 2.6, stroke = 0.4, fill = "steelblue4", colour = "steelblue4") +
  scale_x_log10(breaks = dose_breaks, labels = dose_labels_vec, minor_breaks = NULL) +
  coord_cartesian(ylim = c(0, 105), clip = "off") +
  theme_bw()

ggsave("test_plot.png", p_dr_rel, width = 5, height = 4)
cat("Wrote test_plot.png successfully!\n")
