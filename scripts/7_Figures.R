#===============================================================================================
# INITIAL SETUP
# ===============================================================================================

# Clear environment
rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(grid)
  library(pROC)
  library(survival)
  library(scales)
  library(ggtext)
  library(forestplot)
  library(patchwork)
  library(cowplot)
  
})



# Load custom utility functions
source(here("scripts", "0_functions.R"))


# =========================================================
# 1. Load datasets
# =========================================================
load(here("results", "PIMA_Metrics_T1_T9_chimera.Rdata"))
load(here("results", "PIMA_Metrics_T1_T9_synthpop.Rdata"))
load(here("results", "PIMA_Metrics_T1_T9_ctgan.Rdata"))

load(here("results", "AIDS_Metrics_T1_T9_chimera.Rdata"))
load(here("results", "AIDS_Metrics_T1_T9_synthpop.Rdata"))
load(here("results", "AIDS_Metrics_T1_T9_ctgan.Rdata"))

load(here("results", "REIN_Metrics_T1_T9_chimera.Rdata"))
load(here("results", "REIN_Metrics_T1_T9_synthpop.Rdata"))
load(here("results", "REIN_Metrics_T1_T9_ctgan.Rdata"))


# =========================================================
# 2. NPJ-inspired color palette
# =========================================================
npj_palette <- c(
  "CHIMERA"  = "#00A087",
  "SYNTHPOP" = "#E64B35",
  "CTGAN"    = "#4DBBD5"
)

best_point_fill   <- "#F2C14E"   # Fill color for the best-of-50 diamond marker
best_point_colour <- "#2B2B2B"   # Border color for the best-of-50 marker
ref_line_colour   <- "#8C8C8C"   # Color for the reference (target) dashed line
grid_colour       <- "#D9E0E6"   # Major grid line color
strip_fill        <- "#F3F6F8"   # Facet strip background color
border_colour     <- "#D0D7DE"   # Panel border color
text_colour       <- "#1F2933"   # Global text color


# =========================================================
# 3. Metric metadata
# =========================================================
metric_info <- tibble(
  metric = paste0("T", 1:9),
  metric_label = c(
    "T1: JSD ↓",
    "T2: CFD ↓",
    "T3: DiscPred (0.5)",
    "T4: SDiff ↓",
    "T5: ΔAUC TSTR (0)",
    "T6: ΔAUC TSRTR (0)",
    "T7: MIR-RM ↑",
    "T8: MIR-Holdout (0.5)",
    "T9: AIR (F1) ↓"
  ),
  target = c(0, 0, 0.5, 0, 0, 0, 1, 0.5, 0),   # Ideal reference value for each metric
  metric_group = c(
    "Fidelity", "Fidelity", "Fidelity",
    "Utility",  "Utility",  "Utility",
    "Privacy",  "Privacy",  "Privacy"
  )
)


# =========================================================
# 4. Data preparation — long format
# =========================================================
SYNTH_ORDER <- c("CHIMERA", "SYNTHPOP", "CTGAN")

# Convert one wide-format metrics data frame to long format,
# attach metric metadata, and apply factor levels.
prepare_one_synth <- function(df, dataset_name, synthesizer_name) {
  
  df |>
    dplyr::mutate(
      syn_id       = dataset,
      replicate_id = dplyr::row_number(),
      dataset_name = dataset_name,
      synthesizer  = toupper(synthesizer_name)
    ) |>
    dplyr::select(
      dataset_name,
      synthesizer,
      syn_id,
      replicate_id,
      dplyr::starts_with("T")
    ) |>
    tidyr::pivot_longer(
      cols      = dplyr::starts_with("T"),
      names_to  = "metric",
      values_to = "value"
    ) |>
    dplyr::left_join(metric_info, by = "metric") |>
    dplyr::mutate(
      synthesizer  = factor(synthesizer, levels = SYNTH_ORDER),
      metric       = factor(metric, levels = paste0("T", 1:9)),
      metric_label = factor(metric_label, levels = metric_info$metric_label),
      metric_group = factor(
        metric_group,
        levels = c("Fidelity", "Utility", "Privacy")
      )
    ) |>
    dplyr::filter(is.finite(value))
}

# =========================================================
# 5. Best-of-50 selection — identify the run with minimum L_total
# =========================================================
get_best_syn <- function(df, dataset_name, synthesizer_name) {
  
  df |>
    dplyr::mutate(
      dataset_name = dataset_name,
      synthesizer  = synthesizer_name,
      best_syn_id  = dataset
    ) |>
    dplyr::filter(is.finite(L_total)) |>
    dplyr::slice_min(
      order_by = L_total,
      n = 1,
      with_ties = FALSE
    ) |>
    dplyr::select(
      dataset_name,
      synthesizer,
      best_syn_id,
      L_total
    )
}

# =========================================================
# 6. Build the best-of-50 reference table
# =========================================================
best_of_50 <- bind_rows(
  get_best_syn(AIDS_Metrics_T1_T9_chimera,  "AIDS", "CHIMERA"),
  get_best_syn(AIDS_Metrics_T1_T9_synthpop, "AIDS", "SYNTHPOP"),
  get_best_syn(AIDS_Metrics_T1_T9_ctgan,    "AIDS", "CTGAN"),
  
  get_best_syn(PIMA_Metrics_T1_T9_chimera,  "PIMA", "CHIMERA"),
  get_best_syn(PIMA_Metrics_T1_T9_synthpop, "PIMA", "SYNTHPOP"),
  get_best_syn(PIMA_Metrics_T1_T9_ctgan,    "PIMA", "CTGAN"),
  
  get_best_syn(REIN_Metrics_T1_T9_chimera,  "REIN", "CHIMERA"),
  get_best_syn(REIN_Metrics_T1_T9_synthpop, "REIN", "SYNTHPOP"),
  get_best_syn(REIN_Metrics_T1_T9_ctgan,    "REIN", "CTGAN")
)

print(best_of_50)


# =========================================================
# 7. Build the full long-format plotting dataset
#    and flag the best run for each synthesizer × dataset
# =========================================================
plot_df <- bind_rows(
  prepare_one_synth(AIDS_Metrics_T1_T9_chimera,  "AIDS", "CHIMERA"),
  prepare_one_synth(AIDS_Metrics_T1_T9_synthpop, "AIDS", "SYNTHPOP"),
  prepare_one_synth(AIDS_Metrics_T1_T9_ctgan,    "AIDS", "CTGAN"),
  
  prepare_one_synth(PIMA_Metrics_T1_T9_chimera,  "PIMA", "CHIMERA"),
  prepare_one_synth(PIMA_Metrics_T1_T9_synthpop, "PIMA", "SYNTHPOP"),
  prepare_one_synth(PIMA_Metrics_T1_T9_ctgan,    "PIMA", "CTGAN"),
  
  prepare_one_synth(REIN_Metrics_T1_T9_chimera,  "REIN", "CHIMERA"),
  prepare_one_synth(REIN_Metrics_T1_T9_synthpop, "REIN", "SYNTHPOP"),
  prepare_one_synth(REIN_Metrics_T1_T9_ctgan,    "REIN", "CTGAN")
) |>
  dplyr::left_join(
    best_of_50 |> dplyr::select(dataset_name, synthesizer, best_syn_id),
    by = c("dataset_name", "synthesizer")
  ) |>
  dplyr::mutate(
    is_best    = syn_id == best_syn_id,
    synthesizer = factor(synthesizer, levels = SYNTH_ORDER)
  )

# Quick verification: one best run per synthesizer × dataset combination
plot_df |>
  filter(is_best) |>
  distinct(dataset_name, synthesizer, syn_id, best_syn_id) |>
  arrange(dataset_name, synthesizer) |>
  print()


# =========================================================
# 8. Reference lines — one per metric (ideal target value)
# =========================================================
ref_lines <- metric_info |>
  mutate(
    metric       = factor(metric, levels = paste0("T", 1:9)),
    metric_label = factor(metric_label, levels = metric_info$metric_label),
    metric_group = factor(metric_group, levels = c("Fidelity", "Utility", "Privacy"))
  )


# =========================================================
# 9. Publication-ready ggplot theme
# =========================================================
theme_npj_box <- function(base_size = 11, base_family = "Helvetica") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      text         = element_text(colour = text_colour),
      plot.title   = element_text(face = "bold", size = rel(1.2), hjust = 0),
      plot.subtitle = element_text(size = rel(0.95), hjust = 0, margin = margin(b = 8)),
      
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = rel(1.0), margin = margin(r = 8)),
      axis.text.x  = element_text(face = "bold", size = rel(0.9), colour = text_colour),
      axis.text.y  = element_text(size = rel(0.85), colour = text_colour),
      
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.3),
      
      strip.background = element_rect(fill = strip_fill, colour = NA),
      strip.text       = element_text(face = "bold", size = rel(0.9), colour = text_colour),
      
      panel.border  = element_rect(fill = NA, colour = border_colour, linewidth = 0.5),
      panel.spacing = unit(0.8, "lines"),
      
      legend.position = "none",
      plot.margin     = margin(8, 10, 8, 10)
    )
}


# =========================================================
# 10. Plot one row of metrics (one metric group per row)
# =========================================================
make_metric_row <- function(data, dataset_selected, group_selected, y_lab) {
  
  df_plot       <- data |> filter(dataset_name == dataset_selected, metric_group == group_selected)
  ref_lines_row <- ref_lines |> filter(metric_group == group_selected)
  
  ggplot(df_plot, aes(x = synthesizer, y = value, fill = synthesizer)) +
    
    # Dashed reference line at the ideal target value
    geom_hline(
      data        = ref_lines_row,
      aes(yintercept = target),
      inherit.aes = FALSE,
      linewidth   = 0.35,
      linetype    = "22",
      colour      = ref_line_colour
    ) +
    
    # Boxplot (outliers hidden to reduce clutter)
    geom_boxplot(
      width         = 0.62,
      outlier.shape = NA,
      linewidth     = 0.45,
      colour        = "#3A3A3A",
      alpha         = 0.95
    ) +
    
    # Diamond marker for the best-of-50 run
    geom_point(
      data        = df_plot |> filter(is_best),
      aes(x = synthesizer, y = value),
      inherit.aes = FALSE,
      shape       = 23,
      size        = 2.8,
      stroke      = 0.45,
      fill        = best_point_fill,
      colour      = best_point_colour
    ) +
    
    facet_wrap(~ metric_label, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = npj_palette) +
    scale_x_discrete(labels = SYNTH_ORDER) +
    labs(y = y_lab) +
    coord_cartesian(clip = "off") +
    theme_npj_box() +
    theme(plot.title = element_blank())
}


# =========================================================
# 11. Assemble the full 3-row figure for one dataset
# =========================================================
make_metric_boxplot <- function(data, dataset_selected) {
  p1 <- make_metric_row(data, dataset_selected, "Fidelity", "Fidelity")
  p2 <- make_metric_row(data, dataset_selected, "Utility",  "Utility")
  p3 <- make_metric_row(data, dataset_selected, "Privacy",  "Privacy")
  
  (p1 / p2 / p3) + plot_annotation()
}


# =========================================================
# 12. Generate figures for each dataset
# =========================================================
p_AIDS <- make_metric_boxplot(plot_df, "AIDS")
p_PIMA <- make_metric_boxplot(plot_df, "PIMA")
p_REIN <- make_metric_boxplot(plot_df, "REIN")

p_AIDS
p_PIMA
p_REIN


# =========================================================
# 13. Save figures (300 mm × 200 mm)
# =========================================================
fig_width_in  <- 300 / 25.4
fig_height_in <- 200 / 25.4

ggsave(here("figures", "FigS2_DistrBestofM_PIMA.pdf"), p_PIMA,
       width = fig_width_in, height = fig_height_in, device = cairo_pdf, bg = "white")

ggsave(here("figures", "FigS3_DistrBestofM_AIDS.pdf"), p_AIDS,
       width = fig_width_in, height = fig_height_in, device = cairo_pdf, bg = "white")

ggsave(here("figures", "FigS4_DistrBestofM_REIN.pdf"), p_REIN,
       width = fig_width_in, height = fig_height_in, device = cairo_pdf, bg = "white")


ggsave(here("figures", "FigS2_DistrBestofM_PIMA.png"), p_PIMA,
       width = fig_width_in, height = fig_height_in, dpi=350, bg = "white")

ggsave(here("figures", "FigS3_DistrBestofM_AIDS.png"), p_AIDS,
       width = fig_width_in, height = fig_height_in, dpi=350, bg = "white")

ggsave(here("figures", "FigS4_DistrBestofM_REIN.png"), p_REIN,
       width = fig_width_in, height = fig_height_in, dpi=350, bg = "white")


# ============================================================
# ABLATION STUDY — CHIMERA vs. CHIMERA without matching vs. MICE
# ============================================================

# =========================================================
# 1. Load datasets
# =========================================================
load(here("results", "PIMA_Metrics_T1_T9_chimera.Rdata"))
load(here("results", "PIMA_Metrics_T1_T9_chimera_no_matching.Rdata"))
load(here("results", "PIMA_Metrics_T1_T9_mice.Rdata"))

load(here("results", "AIDS_Metrics_T1_T9_chimera.Rdata"))
load(here("results", "AIDS_Metrics_T1_T9_chimera_no_matching.Rdata"))
load(here("results", "AIDS_Metrics_T1_T9_mice.Rdata"))

load(here("results", "REIN_Metrics_T1_T9_chimera.Rdata"))
load(here("results", "REIN_Metrics_T1_T9_chimera_no_matching.Rdata"))
load(here("results", "REIN_Metrics_T1_T9_mice.Rdata"))


# =========================================================
# 2. NPJ-inspired color palette (ablation study)
# =========================================================
npj_palette <- c(
  "MICE"                = "#4DBBD5",
  "CHIMERA NO MATCHING" = "#E64B35",
  "CHIMERA"             = "#00A087"
)

best_point_fill   <- "#F2C14E"
best_point_colour <- "#2B2B2B"
ref_line_colour   <- "#8C8C8C"
grid_colour       <- "#D9E0E6"
strip_fill        <- "#F3F6F8"
border_colour     <- "#D0D7DE"
text_colour       <- "#1F2933"


# =========================================================
# 3. Metric metadata (same structure as Section 3 above)
# =========================================================
metric_info <- tibble(
  metric = paste0("T", 1:9),
  metric_label = c(
    "T1: JSD ↓", "T2: CFD ↓", "T3: DiscPred (0.5)",
    "T4: SDiff ↓", "T5: ΔAUC TSTR (0)", "T6: ΔAUC TSRTR (0)",
    "T7: MIR-RM ↑", "T8: MIR-Holdout (0.5)", "T9: AIR (F1) ↓"
  ),
  target       = c(0, 0, 0.5, 0, 0, 0, 1, 0.5, 0),
  metric_group = c(
    "Fidelity", "Fidelity", "Fidelity",
    "Utility",  "Utility",  "Utility",
    "Privacy",  "Privacy",  "Privacy"
  )
)


# =========================================================
# 4. Data preparation (long format)
# =========================================================
SYNTH_ORDER <- c("MICE", "CHIMERA NO MATCHING", "CHIMERA")

# =========================================================
# 6. Build the best-of-50 reference table (ablation)
# =========================================================
best_of_50 <- bind_rows(
  get_best_syn(AIDS_Metrics_T1_T9_chimera,             "AIDS", "CHIMERA"),
  get_best_syn(AIDS_Metrics_T1_T9_chimera_no_matching, "AIDS", "CHIMERA NO MATCHING"),
  get_best_syn(AIDS_Metrics_T1_T9_mice,                "AIDS", "MICE"),
  
  get_best_syn(PIMA_Metrics_T1_T9_chimera,             "PIMA", "CHIMERA"),
  get_best_syn(PIMA_Metrics_T1_T9_chimera_no_matching, "PIMA", "CHIMERA NO MATCHING"),
  get_best_syn(PIMA_Metrics_T1_T9_mice,                "PIMA", "MICE"),
  
  get_best_syn(REIN_Metrics_T1_T9_chimera,             "REIN", "CHIMERA"),
  get_best_syn(REIN_Metrics_T1_T9_chimera_no_matching, "REIN", "CHIMERA NO MATCHING"),
  get_best_syn(REIN_Metrics_T1_T9_mice,                "REIN", "MICE")
)

print(best_of_50)


# =========================================================
# 7. Build the full long-format plotting dataset (ablation)
# =========================================================
plot_df <- bind_rows(
  prepare_one_synth(AIDS_Metrics_T1_T9_chimera,             "AIDS", "CHIMERA"),
  prepare_one_synth(AIDS_Metrics_T1_T9_chimera_no_matching, "AIDS", "CHIMERA NO MATCHING"),
  prepare_one_synth(AIDS_Metrics_T1_T9_mice,                "AIDS", "MICE"),
  
  prepare_one_synth(PIMA_Metrics_T1_T9_chimera,             "PIMA", "CHIMERA"),
  prepare_one_synth(PIMA_Metrics_T1_T9_chimera_no_matching, "PIMA", "CHIMERA NO MATCHING"),
  prepare_one_synth(PIMA_Metrics_T1_T9_mice,                "PIMA", "MICE"),
  
  prepare_one_synth(REIN_Metrics_T1_T9_chimera,             "REIN", "CHIMERA"),
  prepare_one_synth(REIN_Metrics_T1_T9_chimera_no_matching, "REIN", "CHIMERA NO MATCHING"),
  prepare_one_synth(REIN_Metrics_T1_T9_mice,                "REIN", "MICE")
) |>
  dplyr::left_join(
    best_of_50 |> dplyr::select(dataset_name, synthesizer, best_syn_id),
    by = c("dataset_name", "synthesizer")
  ) |>
  dplyr::mutate(
    is_best     = syn_id == best_syn_id,
    synthesizer = factor(synthesizer, levels = SYNTH_ORDER)
  )


# =========================================================
# 8. Reference lines (ablation)
# =========================================================
ref_lines <- metric_info |>
  mutate(
    metric       = factor(metric, levels = paste0("T", 1:9)),
    metric_label = factor(metric_label, levels = metric_info$metric_label),
    metric_group = factor(metric_group, levels = c("Fidelity", "Utility", "Privacy"))
  )



# =========================================================
# 12. Generate figures (ablation study)
# =========================================================
p_AIDS <- make_metric_boxplot(plot_df, "AIDS")
p_PIMA <- make_metric_boxplot(plot_df, "PIMA")
p_REIN <- make_metric_boxplot(plot_df, "REIN")

p_AIDS
p_PIMA
p_REIN


# =========================================================
# 13. Save figures (ablation study)
# =========================================================
fig_width_in  <- 300 / 25.4
fig_height_in <- 200 / 25.4

ggsave(here("figures", "FigS5_AblationStudy_PIMA.pdf"),
       p_PIMA, width = fig_width_in, height = fig_height_in, device = cairo_pdf, bg = "white")

ggsave(here("figures", "FigS6_AblationStudy_AIDS.pdf"),
       p_AIDS, width = fig_width_in, height = fig_height_in, device = cairo_pdf, bg = "white")

ggsave(here("figures", "FigS7_AblationStudy_REIN.pdf"),
       p_REIN, width = fig_width_in, height = fig_height_in, device = cairo_pdf, bg = "white")

ggsave(here("figures", "FigS5_AblationStudy_PIMA.png"),
       p_PIMA, width = fig_width_in, height = fig_height_in, dpi = 350, bg = "white")

ggsave(here("figures", "FigS6_AblationStudy_AIDS.png"),
       p_AIDS, width = fig_width_in, height = fig_height_in, dpi = 350, bg = "white")

ggsave(here("figures", "FigS7_AblationStudy_REIN.png"),
       p_REIN, width = fig_width_in, height = fig_height_in, dpi = 350, bg = "white")


# ============================================================
# UNIVARIATE DISTRIBUTION — CHIMERA vs. SYNTHPOP vs. CTGAN
# ============================================================

# =========================================================
# PIMA — Load data
# =========================================================
load(here("data", "PIMA.Rdata"))
load(here("data", "PIMA_CHIMERA.Rdata"))
load(here("data", "PIMA_SYNTHPOP.Rdata"))
load(here("data", "PIMA_CTGAN.Rdata"))


# =========================================================
# NPJ color palette (includes real data)
# =========================================================
npj_palette <- c(
  "REAL"     = "#000000",
  "CHIMERA"  = "#00A087",
  "SYNTHPOP" = "#E64B35",
  "CTGAN"    = "#4DBBD5"
)

# Line types: REAL uses dashed to stand out from synthetic curves
npj_lty <- c(
  "REAL"     = "dashed",
  "CHIMERA"  = "solid",
  "SYNTHPOP" = "dashed",
  "CTGAN"    = "solid"
)

# Line widths: REAL is thicker to ensure visual prominence
npj_lwd <- c(
  "REAL"     = 1.2,
  "CHIMERA"  = 0.6,
  "SYNTHPOP" = 0.6,
  "CTGAN"    = 0.6
)


# =========================================================
# Units of measurement for numeric variables (used in axis labels)
# =========================================================
variable_units <- list(
  Plasma_glucose_concentration = "mg/dL",
  Diastolic_blood_pressure     = "mmHg",
  Triceps_skin_fold_thickness  = "mm",
  Two_hour_serum_insulin       = "µU/mL",
  Body_mass_index              = "kg/m\u00b2",
  Age                          = "years"
)

# Helper: build HTML-formatted strip labels with variable name + unit
format_label <- function(value) {
  vapply(value, function(v) {
    var_clean <- gsub("_", " ", v)
    unit      <- variable_units[[v]]
    if (!is.null(unit)) {
      paste0("<b>", var_clean, "</b><br><span style='font-size:6.5pt'>(", unit, ")</span>")
    } else {
      paste0("<b>", var_clean, "</b>")
    }
  }, character(1), USE.NAMES = FALSE)
}


# =========================================================
# NPJ theme for univariate distribution panels
# =========================================================
theme_npj <- function(base_size = 7) {
  theme_minimal(base_size = base_size, base_family = "Helvetica") +
    theme(
      text       = element_text(colour = "#1F2933"),
      plot.title = element_text(size = rel(1.05), hjust = 0, face = "bold", margin = margin(b = 4)),
      
      # Axes
      axis.title        = element_blank(),
      axis.text         = element_text(size = rel(0.85), colour = "#3A3A3A"),
      axis.text.x       = element_text(size = rel(0.80)),
      axis.ticks        = element_line(colour = "#D0D7DE", linewidth = 0.25),
      axis.ticks.length = unit(1.5, "pt"),
      
      # Grid
      panel.grid.major.y = element_line(colour = "#D9E0E6", linewidth = 0.2),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      
      # Facet strips — element_markdown() required to render inline HTML labels
      strip.background = element_rect(fill = "#F3F6F8", colour = NA),
      strip.text       = element_markdown(
        size = rel(0.90), colour = "#1F2933", halign = 0.5,
        margin = margin(t = 3, b = 3), lineheight = 1.35
      ),
      
      # Panel border
      panel.border  = element_rect(fill = NA, colour = "#D0D7DE", linewidth = 0.3),
      panel.spacing = unit(0.55, "lines"),
      
      # Legend
      legend.position  = "top",
      legend.direction = "horizontal",
      legend.title     = element_blank(),
      legend.text      = element_text(size = rel(0.95)),
      legend.key.size  = unit(8,  "pt"),
      legend.key.width = unit(12, "pt"),
      legend.spacing.x = unit(6,  "pt"),
      legend.margin    = margin(t = 2, b = 2),
      
      plot.margin = margin(4, 6, 4, 6)
    )
}


# =========================================================
# Merge real and synthetic PIMA datasets
# =========================================================
df_pima <- bind_rows(
  PIMA          |> mutate(Dataset = "REAL"),
  PIMA_CHIMERA  |> mutate(Dataset = "CHIMERA"),
  PIMA_SYNTHPOP |> mutate(Dataset = "SYNTHPOP"),
  PIMA_CTGAN    |> mutate(Dataset = "CTGAN")
) |>
  mutate(Dataset = factor(Dataset, levels = c("REAL", "CHIMERA", "SYNTHPOP", "CTGAN")))


# =========================================================
# Panel A — Numeric variable density plots
# =========================================================
df_long_num <- df_pima |>
  select(where(is.numeric), Dataset) |>
  pivot_longer(-Dataset, names_to = "variable", values_to = "value") |>
  filter(is.finite(value))

n_num  <- df_long_num |> distinct(variable) |> nrow()
ncol_A <- min(4L, n_num)   # Maximum 4 columns

p_num <- ggplot(df_long_num, aes(x = value, colour = Dataset, linetype = Dataset, linewidth = Dataset)) +
  geom_density(adjust = 1.1) +
  facet_wrap(~ variable, scales = "free", ncol = ncol_A, labeller = labeller(variable = format_label)) +
  scale_colour_manual(values = npj_palette) +
  scale_linetype_manual(values = npj_lty) +
  scale_linewidth_manual(values = npj_lwd) +
  scale_x_continuous(labels = scales::label_number(accuracy = NULL, big.mark = "\u202f")) +
  labs(title = "a \u2014 Univariate numerical distributions") +
  theme_npj() +
  guides(colour = "none", linetype = "none", linewidth = "none") +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())


# =========================================================
# Panel B — Categorical variable bar charts (proportions)
# =========================================================
cat_vars <- df_pima |>
  select(where(~ is.factor(.) | is.character(.))) |>
  names()

df_cat_prop <- df_pima |>
  select(all_of(cat_vars), Dataset) |>
  pivot_longer(-Dataset, names_to = "variable", values_to = "value") |>
  filter(!is.na(value)) |>
  mutate(value = as.factor(value)) |>
  group_by(Dataset, variable, value) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(Dataset, variable) |>
  mutate(prop = n / sum(n)) |>
  ungroup()

n_cat  <- df_cat_prop |> distinct(variable) |> nrow()
ncol_B <- min(1L, n_cat)   # One column for PIMA (single categorical variable)

p_cat <- ggplot(df_cat_prop, aes(x = value, y = prop, fill = Dataset)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.68, colour = "#2C2C2C", linewidth = 0.20) +
  facet_wrap(~ variable, scales = "free_x", ncol = ncol_B, labeller = labeller(variable = format_label)) +
  scale_fill_manual(values = npj_palette) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.05))) +
  labs(title = "b \u2014 Univariate categorical distributions") +
  guides(fill = guide_legend(override.aes = list(alpha = 1, colour = "#2C2C2C", linewidth = 0.2)), colour = "none") +
  theme_npj() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = rel(0.78)))


# =========================================================
# Assemble panels A and B with patchwork
# =========================================================
height_ratio <- c(ceiling(n_num / ncol_A), ceiling(n_cat / ncol_B))

p_final <- p_num / p_cat +
  plot_layout(heights = height_ratio, guides = "collect") &
  theme(legend.position = "bottom")

p_final


# =========================================================
# Save — Nature standard width (~180 mm)
# =========================================================
w_in <- 180 / 25.4   # 7.087 in
h_in <- 220 / 25.4   # 8.661 in

ggsave(
  filename = here("figures", "Fig1_UnivariateDistr_PIMA.pdf"),
  plot = p_final, width = w_in, height = h_in, units = "in", device = cairo_pdf, bg = "white"
)

ggsave(
  filename = here("figures", "Fig1_UnivariateDistr_PIMA.png"),
  plot = p_final, width = w_in, height = h_in, units = "in", dpi = 350, bg = "white"
)


# ============================================================
# AIDS — Univariate distribution comparison
# ============================================================

# =========================================================
# Load data
# =========================================================
load(here("data", "AIDS.Rdata"))
load(here("data", "AIDS_CHIMERA.Rdata"))
load(here("data", "AIDS_SYNTHPOP.Rdata"))
load(here("data", "AIDS_CTGAN.Rdata"))

# Recode Censorship as a labeled factor (0 = Censored, 1 = Observed)
for (df_name in c("AIDS", "AIDS_CHIMERA", "AIDS_SYNTHPOP", "AIDS_CTGAN")) {
  df <- get(df_name)
  df$Censored <- factor(df$Censored, levels = c(0, 1), labels = c("Censored", "Observed"))
  assign(df_name, df)
}

# Units of measurement for AIDS numeric variables
variable_units <- list(
  Age             = "years",
  Weight          = "kg",
  Karnofsky_score = "%",
  CD4_baseline    = "cells/mm³",
  CD8_baseline    = "cells/mm³"
)

# Merge datasets
df_aids <- bind_rows(
  AIDS          |> mutate(Dataset = "REAL"),
  AIDS_CHIMERA  |> mutate(Dataset = "CHIMERA"),
  AIDS_SYNTHPOP |> mutate(Dataset = "SYNTHPOP"),
  AIDS_CTGAN    |> mutate(Dataset = "CTGAN")
) |>
  mutate(Dataset = factor(Dataset, levels = c("REAL", "CHIMERA", "SYNTHPOP", "CTGAN")))

# Exclude time-to-event variables from distribution panels
df_aids_bis <- df_aids |> select(-c("Censored", "times"))

# Panel A — Numeric distributions
df_long_num <- df_aids_bis |>
  select(where(is.numeric), Dataset) |>
  pivot_longer(-Dataset, names_to = "variable", values_to = "value") |>
  filter(is.finite(value))

n_num  <- df_long_num |> distinct(variable) |> nrow()
ncol_A <- min(4L, n_num)

p_num <- ggplot(df_long_num, aes(x = value, colour = Dataset, linetype = Dataset, linewidth = Dataset)) +
  geom_density(adjust = 1.1) +
  facet_wrap(~ variable, scales = "free", ncol = ncol_A, labeller = labeller(variable = format_label)) +
  scale_colour_manual(values = npj_palette) +
  scale_linetype_manual(values = npj_lty) +
  scale_linewidth_manual(values = npj_lwd) +
  scale_x_continuous(labels = scales::label_number(accuracy = NULL, big.mark = "\u202f")) +
  labs(title = "a \u2014 Univariate numerical distributions") +
  theme_npj() +
  guides(colour = "none", linetype = "none", linewidth = "none") +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Panel B — Categorical distributions
cat_vars <- df_aids_bis |> select(where(~ is.factor(.) | is.character(.))) |> names()

df_cat_prop <- df_aids_bis |>
  select(all_of(cat_vars), Dataset) |>
  pivot_longer(-Dataset, names_to = "variable", values_to = "value") |>
  filter(!is.na(value)) |>
  mutate(value = as.factor(value)) |>
  group_by(Dataset, variable, value) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(Dataset, variable) |>
  mutate(prop = n / sum(n)) |>
  ungroup()

n_cat  <- df_cat_prop |> distinct(variable) |> nrow()
ncol_B <- min(4L, n_cat)

p_cat <- ggplot(df_cat_prop, aes(x = value, y = prop, fill = Dataset)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.68, colour = "#2C2C2C", linewidth = 0.20) +
  facet_wrap(~ variable, scales = "free_x", ncol = ncol_B, labeller = labeller(variable = format_label)) +
  scale_fill_manual(values = npj_palette) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.05))) +
  labs(title = "b \u2014 Univariate categorical distributions") +
  guides(fill = guide_legend(override.aes = list(alpha = 1, colour = "#2C2C2C", linewidth = 0.2)), colour = "none") +
  theme_npj() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = rel(0.78)))

# Assemble and save
height_ratio <- c(ceiling(n_num / ncol_A), ceiling(n_cat / ncol_B))
p_final <- p_num / p_cat + plot_layout(heights = height_ratio, guides = "collect") & theme(legend.position = "bottom")
p_final

ggsave(filename=here("figures", "Fig1_UnivariateDistr_AIDS.pdf"),
       p_final, width = w_in, height = h_in, units = "in", device = cairo_pdf, bg = "white")

ggsave(filename=here("figures", "Fig1_UnivariateDistr_AIDS.png"),
       p_final, width = w_in, height = h_in, units = "in", dpi = 350, bg = "white")


# ============================================================
# REIN — Univariate distribution comparison
# ============================================================

# Load data
load(here("data", "REIN.Rdata"))
load(here("data", "REIN_CHIMERA.Rdata"))
load(here("data", "REIN_SYNTHPOP.Rdata"))
load(here("data", "REIN_CTGAN.Rdata"))

# Units of measurement for REIN numeric variables
variable_units <- list(
  Serum_albumin_level = "g/dL",
  Body_mass_index     = "kg/m\u00b2",
  eGFR                = "mL/min/1.73m\u00b2",
  Age                 = "years"
)

# Merge datasets
df_rein <- bind_rows(
  REIN          |> mutate(Dataset = "REAL"),
  REIN_CHIMERA  |> mutate(Dataset = "CHIMERA"),
  REIN_SYNTHPOP |> mutate(Dataset = "SYNTHPOP"),
  REIN_CTGAN    |> mutate(Dataset = "CTGAN")
) |>
  mutate(Dataset = factor(Dataset, levels = c("REAL", "CHIMERA", "SYNTHPOP", "CTGAN")))

# Exclude time-to-event variables from distribution panels
df_rein_bis <- df_rein |> select(-c("Censored", "times"))

# Panel A — Numeric distributions
df_long_num <- df_rein_bis |>
  select(where(is.numeric), Dataset) |>
  pivot_longer(-Dataset, names_to = "variable", values_to = "value") |>
  filter(is.finite(value))

n_num  <- df_long_num |> distinct(variable) |> nrow()
ncol_A <- min(4L, n_num)

p_num <- ggplot(df_long_num, aes(x = value, colour = Dataset, linetype = Dataset, linewidth = Dataset)) +
  geom_density(adjust = 1.1) +
  facet_wrap(~ variable, scales = "free", ncol = ncol_A, labeller = labeller(variable = format_label)) +
  scale_colour_manual(values = npj_palette) +
  scale_linetype_manual(values = npj_lty) +
  scale_linewidth_manual(values = npj_lwd) +
  scale_x_continuous(labels = scales::label_number(accuracy = NULL, big.mark = "\u202f")) +
  labs(title = "a \u2014 Univariate numerical distributions") +
  theme_npj() +
  guides(colour = "none", linetype = "none", linewidth = "none") +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Panel B — Categorical distributions
cat_vars <- df_rein_bis |> select(where(~ is.factor(.) | is.character(.))) |> names()

df_cat_prop <- df_rein_bis |>
  select(all_of(cat_vars), Dataset) |>
  pivot_longer(-Dataset, names_to = "variable", values_to = "value") |>
  filter(!is.na(value)) |>
  mutate(value = as.factor(value)) |>
  group_by(Dataset, variable, value) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(Dataset, variable) |>
  mutate(prop = n / sum(n)) |>
  ungroup()

n_cat  <- df_cat_prop |> distinct(variable) |> nrow()
ncol_B <- min(4L, n_cat)

p_cat <- ggplot(df_cat_prop, aes(x = value, y = prop, fill = Dataset)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.68, colour = "#2C2C2C", linewidth = 0.20) +
  facet_wrap(~ variable, scales = "free_x", ncol = ncol_B, labeller = labeller(variable = format_label)) +
  scale_fill_manual(values = npj_palette) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.05))) +
  labs(title = "b \u2014 Univariate categorical distributions") +
  guides(fill = guide_legend(override.aes = list(alpha = 1, colour = "#2C2C2C", linewidth = 0.2)), colour = "none") +
  theme_npj() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = rel(0.78)))

# Assemble and save
height_ratio <- c(ceiling(n_num / ncol_A), ceiling(n_cat / ncol_B))
p_final <- p_num / p_cat + plot_layout(heights = height_ratio, guides = "collect") & theme(legend.position = "bottom")
p_final

ggsave(filename=here("figures", "Fig1_UnivariateDistr.pdf"),
       p_final, width = w_in, height = h_in, units = "in", device = cairo_pdf, bg = "white")

ggsave(filename=here("figures", "Fig1_UnivariateDistr.png"),
       p_final, width = w_in, height = h_in, units = "in", dpi = 350, bg = "white")


# ============================================================
# PIMA — ROC curve and forest plot (logistic regression)
# ============================================================

# ── Shared graphic settings ────────────────────────────────
npj_palette <- c(
  "REAL"     = "#000000",
  "CHIMERA"  = "#00A087",
  "SYNTHPOP" = "#E64B35",
  "CTGAN"    = "#4DBBD5"
)

npj_linetypes <- c(
  "REAL"     = "solid",
  "CHIMERA"  = "solid",
  "SYNTHPOP" = "solid",
  "CTGAN"    = "solid"
)

# Updated NPJ theme for larger panels (base_size = 12)
theme_npj <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      text          = element_text(colour = "#000000"),
      plot.title    = element_text(size = rel(1.05), hjust = 0, margin = margin(b = 5)),
      axis.title    = element_text(size = rel(1.0), family = "sans"),
      axis.text     = element_text(size = rel(1),   family = "sans"),
      axis.ticks        = element_line(linewidth = 0.25),
      axis.ticks.length = unit(2, "pt"),
      panel.grid.major  = element_line(colour = "#E8ECF0", linewidth = 0.3),
      panel.grid.minor  = element_blank(),
      panel.border      = element_rect(fill = NA, colour = "#000000", linewidth = 0.5),
      legend.position        = "inside",
      legend.position.inside = c(0.95, 0.05),
      legend.justification   = c(1, 0),
      legend.background = element_rect(
        fill = scales::alpha("white", 0.95), colour = "#B0B0B0", linewidth = 0.3
      ),
      legend.title      = element_blank(),
      legend.text       = element_text(size = 10, family = "sans"),
      legend.key.width  = unit(0.7, "cm"),
      legend.key.height = unit(0.30, "cm"),
      legend.spacing.y  = unit(0.05, "cm"),
      legend.margin     = margin(3, 5, 3, 5),
      plot.margin       = margin(5, 8, 5, 5)
    )
}


# =========================================================
# ROC curve — PIMA
# =========================================================
load(here("results", "PIMA_tbl_OR.Rdata"))
load(here("results", "pima.roc.original.Rdata"))
load(here("results", "pima.roc.chimera.Rdata"))
load(here("results", "pima.roc.synthpop.Rdata"))
load(here("results", "pima.roc.ctgan.Rdata"))

# Helper: compute AUC and 95% CI for one ROC object
compute_auc_label <- function(roc_obj, name) {
  auc_val <- as.numeric(auc(roc_obj))
  ci_val  <- as.numeric(ci.auc(roc_obj, conf.level = 0.95))
  tibble(Dataset = name, AUC = auc_val, CI_lo = ci_val[1], CI_hi = ci_val[3])
}

auc_tbl <- bind_rows(
  compute_auc_label(pima.roc.original, "REAL"),
  compute_auc_label(pima.roc.chimera,  "CHIMERA"),
  compute_auc_label(pima.roc.synthpop, "SYNTHPOP"),
  compute_auc_label(pima.roc.ctgan,    "CTGAN")
) |>
  mutate(
    Dataset = factor(Dataset, levels = c("REAL", "CHIMERA", "SYNTHPOP", "CTGAN")),
    label   = sprintf("%-8s  AUC = %.3f  [95%%\u202fCI: %.3f\u2013%.3f]", Dataset, AUC, CI_lo, CI_hi)
  )

# Build long-format ROC data for ggplot
roc_list <- list(
  REAL     = pima.roc.original,
  CHIMERA  = pima.roc.chimera,
  SYNTHPOP = pima.roc.synthpop,
  CTGAN    = pima.roc.ctgan
)

df_roc <- lapply(names(roc_list), function(nm) {
  r <- roc_list[[nm]]
  tibble(Dataset = nm, Specificity = r$specificities, Sensitivity = r$sensitivities)
}) |>
  bind_rows() |>
  mutate(
    Dataset = factor(Dataset, levels = c("REAL", "CHIMERA", "SYNTHPOP", "CTGAN")),
    FPR     = 1 - Specificity   # Conventional x-axis: 1 - Specificity
  )

# ROC plot (AUC values added as inline text annotations)
plt_roc <- ggplot(df_roc, aes(x = FPR, y = Sensitivity, colour = Dataset, linetype = Dataset)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.4, colour = "#000000") +
  geom_line(linewidth = 0.65) +
  scale_colour_manual(values = npj_palette) +
  scale_linetype_manual(values = npj_linetypes) +
  scale_x_continuous(name = "1 \u2212 Specificity", limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     labels = scales::number_format(accuracy = 0.1), expand = expansion(mult = 0.01)) +
  scale_y_continuous(name = "Sensitivity", limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     labels = scales::number_format(accuracy = 0.1), expand = expansion(mult = 0.01)) +
  labs(title = "") +
  theme_npj() +
  theme(legend.position = "none") +
  # AUC annotations (manually positioned to avoid legend overlap)
  annotate("text", x = 0.59, y = 0.295, label = "REAL",     hjust = 0, size = 5, colour = npj_palette["REAL"]) +
  annotate("text", x = 0.70, y = 0.295, label = "AUC: 0.839 [95% CI: 0.811\u20130.868]", hjust = 0, size = 5) +
  annotate("text", x = 0.59, y = 0.262, label = "CHIMERA",  hjust = 0, size = 5, colour = npj_palette["CHIMERA"]) +
  annotate("text", x = 0.70, y = 0.262, label = "AUC: 0.839 [95% CI: 0.808\u20130.869]", hjust = 0, size = 5) +
  annotate("text", x = 0.59, y = 0.229, label = "SYNTHPOP", hjust = 0, size = 5, colour = npj_palette["SYNTHPOP"]) +
  annotate("text", x = 0.70, y = 0.229, label = "AUC: 0.852 [95% CI: 0.824\u20130.881]", hjust = 0, size = 5) +
  annotate("text", x = 0.59, y = 0.196, label = "CTGAN",    hjust = 0, size = 5, colour = npj_palette["CTGAN"]) +
  annotate("text", x = 0.70, y = 0.196, label = "AUC: 0.712 [95% CI: 0.676\u20130.748]", hjust = 0, size = 5)


plt_roc

# =========================================================
# Forest plot — odds ratios from logistic regression (PIMA)
# =========================================================
tbl_df <- PIMA_tbl_OR

# Helper: format OR or HR with 95% CI as a string
fmt <- function(est, lo, hi) sprintf("%.3f [%.3f; %.3f]", est, lo, hi)

# Build the matrix of point estimates and confidence intervals
# DATA
RRforest_combined <- cbind.data.frame(
  OR_original    = c(NA, tbl_df$OR_original),
  lower_original = c(NA, tbl_df$IC95_low_original),
  upper_original = c(NA, tbl_df$IC95_high_original),
  
  OR_chimera     = c(NA, tbl_df$OR_chimera),
  lower_chimera  = c(NA, tbl_df$IC95_low_chimera),
  upper_chimera  = c(NA, tbl_df$IC95_high_chimera),
  
  OR_synthpop    = c(NA, tbl_df$OR_synthpop),
  lower_synthpop = c(NA, tbl_df$IC95_low_synthpop),
  upper_synthpop = c(NA, tbl_df$IC95_high_synthpop),
  
  OR_ctgan       = c(NA, tbl_df$OR_ctgan),
  lower_ctgan    = c(NA, tbl_df$IC95_low_ctgan),
  upper_ctgan    = c(NA, tbl_df$IC95_high_ctgan)
)

# COVARIATES
covariates_text <- c(
  "Covariate",
  "Age (10 years)",
  "Body mass index",
  "Diabetes pedigree function",
  "Diastolic blood pressure",
  "Number of pregnancies",
  "Plasma glucose",
  "Triceps skinfold thickness",
  "Two hour serum insulin"
)

# OR TEXT 
OR_text <- c(
  "Odds Ratio [95% CI]",
  paste(
    sprintf("%-14s : %s", "REAL",
            fmt(tbl_df$OR_original,
                tbl_df$IC95_low_original,
                tbl_df$IC95_high_original)),
    
    sprintf("%-10s : %s", "CHIMERA",
            fmt(tbl_df$OR_chimera,
                tbl_df$IC95_low_chimera,
                tbl_df$IC95_high_chimera)),
    
    sprintf("%-1s : %s", "SYNTHPOP",
            fmt(tbl_df$OR_synthpop,
                tbl_df$IC95_low_synthpop,
                tbl_df$IC95_high_synthpop)),
    
    sprintf("%-12s : %s", "CTGAN",
            fmt(tbl_df$OR_ctgan,
                tbl_df$IC95_low_ctgan,
                tbl_df$IC95_high_ctgan)),
    
    sep = "\n"
  )
)

# SDIFF TEXT (même format)
SD_text <- c(
  "SDiff",
  paste(
    sprintf("%s %s", "",""),  # ✔ ligne vide pour alignement
    sprintf("%s %.3f", "",  tbl_df$SMD_log1),
    sprintf("%s %.3f", "", tbl_df$SMD_log2),
    sprintf("%s %.3f", "",    tbl_df$SMD_log3),
    sep = "\n"
  )
)

# TABLE TEXT
tabletext <- cbind.data.frame(covariates_text, OR_text, SD_text)


# FOREST PLOT
FP.plot <- forestplot(
  tabletext,
  
  mean = cbind(
    RRforest_combined$OR_original,
    RRforest_combined$OR_chimera,
    RRforest_combined$OR_synthpop,
    RRforest_combined$OR_ctgan
  ),
  
  lower = cbind(
    RRforest_combined$lower_original,
    RRforest_combined$lower_chimera,
    RRforest_combined$lower_synthpop,
    RRforest_combined$lower_ctgan
  ),
  
  upper = cbind(
    RRforest_combined$upper_original,
    RRforest_combined$upper_chimera,
    RRforest_combined$upper_synthpop,
    RRforest_combined$upper_ctgan
  ),
  
  zero = 1,
  Xlog = TRUE,
  xticks = c(0.75, 1, 1.25, 1.5, 1.75, 2),
  graph.pos = 2,
  align = c("l", "l", "r"),
  
  is.summary = c(TRUE, rep(FALSE, nrow(RRforest_combined) - 1)),
  
  graphwidth = unit(14.5, "cm"),
  boxsize = 0.1,
  lwd.ci=1.5,
  
  colgap = unit(2, "mm"),
  line.margin = unit(2.5, "mm"),
  
  txt_gp = fpTxtGp(
    label = gpar(
      fontfamily = "sans",
      cex = 1,
      lineheight = 1,
      fontface = "plain"
    ),
    ticks = gpar(
      cex = 1,
      fontface = "plain"
    ),
    xlab = gpar(
      cex = 1,
      fontface = "plain"
    ),
    summary = gpar(
      fontface = "plain"
    )
  ),
  
  xlab = "Odds Ratio [95% CI]",
  
  col = fpColors(
    box = c("#000000","#00A087","#E64B35","#4DBBD5"),
    line = c("#000000","#00A087","#E64B35","#4DBBD5"),
    summary = "#000000"
  )
) |>
  fp_add_lines(
    h_2 = gpar(col = "#000000", lwd = 1)
  ) |>
  fp_set_zebra_style("#E3EDF5")

FP.plot

## Combine AUC + Forestplot ----
g_forest <- grid.grabExpr(print(FP.plot))

plotAUC_forest_PIMA <- plt_roc / wrap_elements(full = g_forest) +
  plot_layout(heights = c(0.7, 1), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(family = "sans", face = "bold", size = 16),
        plot.tag.position = c(0.01, 0.99))

ggsave(
  filename = here("figures", "Fig2_AUC_Forest_PIMA.pdf"),
  plot     = plotAUC_forest_PIMA,
  width    = 300 / 25.4,
  height   = (200 / 25.4) * 2,
  units    = "in", scale = 1
)


# ============================================================
# AIDS — Kaplan–Meier, time-dependent AUC, and forest plot
# ============================================================

# ── Kaplan–Meier curves ──────────────────────────────────────
load(here("results", "km_aids.Rdata"))

p_km <- ggplot(km_aids, aes(x = time, y = surv, colour = group, linetype = group, fill = group)) +
  geom_step(linewidth = 1) +
  scale_colour_manual(values = npj_palette, name = "Dataset") +
  scale_linetype_manual(values = npj_linetypes, name = "Dataset") +
  scale_fill_manual(values = npj_palette, guide = "none") +
  scale_x_continuous(name = "Follow-up time (days)", limits = c(0, NA),
                     breaks = seq(0, 1400, by = 200), expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(name = "Survival probability", limits = c(0.50, 1.00),
                     breaks = seq(0.5, 1.0, by = 0.1),
                     labels = scales::number_format(accuracy = 0.1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "") +
  guides(colour = "none", linetype = "none") +
  theme_npj() +
  theme(legend.key.width = unit(1.4, "cm"), legend.key.height = unit(0.5, "cm"))



# ── Time-dependent AUC ───────────────────────────────────────
load(here("results", "auc_all_aids.Rdata"))

auc_all <- auc_all_aids |>
  mutate(source = factor(source, levels = c("REAL", "CHIMERA", "SYNTHPOP", "CTGAN")))

# Compute integrated AUC (iAUC) and mean 95% CI per method
auc_means <- auc_all |>
  group_by(source) |>
  summarise(mean_auc = mean(auc, na.rm = TRUE), mean_lower = mean(lower, na.rm = TRUE),
            mean_upper = mean(upper, na.rm = TRUE)) |>
  ungroup() |>
  mutate(source = factor(source, levels = c("REAL", "CHIMERA", "SYNTHPOP", "CTGAN"))) |>
  arrange(source)

p_auc <- ggplot(auc_all, aes(x = time, y = auc, colour = source, linetype = source, fill = source)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = npj_palette, name = NULL) +
  scale_linetype_manual(values = npj_linetypes, name = NULL) +
  scale_fill_manual(values = npj_palette, guide = "none") +
  scale_x_continuous(name = "Follow-up time (days)", breaks = seq(0, 1100, by = 100),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(name = "Time-dependent AUC", limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.2), labels = number_format(accuracy = 0.1),
                     expand = expansion(mult = c(0, 0.02))) +
  # iAUC annotations
  annotate("text", x = 300, y = 0.295, label = "REAL",     hjust = 0, size = 4, colour = npj_palette["REAL"]) +
  annotate("text", x = 550, y = 0.295, label = "iAUC: 0.672 [95% CI: 0.658\u20130.682]", hjust = 0, size = 4) +
  annotate("text", x = 300, y = 0.262, label = "CHIMERA",  hjust = 0, size = 4, colour = npj_palette["CHIMERA"]) +
  annotate("text", x = 550, y = 0.262, label = "iAUC: 0.687 [95% CI: 0.676\u20130.695]", hjust = 0, size = 4) +
  annotate("text", x = 300, y = 0.229, label = "SYNTHPOP", hjust = 0, size = 4, colour = npj_palette["SYNTHPOP"]) +
  annotate("text", x = 550, y = 0.229, label = "iAUC: 0.645 [95% CI: 0.633\u20130.653]", hjust = 0, size = 4) +
  annotate("text", x = 300, y = 0.196, label = "CTGAN",    hjust = 0, size = 4, colour = npj_palette["CTGAN"]) +
  annotate("text", x = 550, y = 0.196, label = "iAUC: 0.617 [95% CI: 0.609\u20130.623]", hjust = 0, size = 4) +
  labs(title = "") + theme_npj() + theme(legend.position = "none")


# ── Forest plot — hazard ratios from Cox model (AIDS) ────────
load(here("results", "tbl_AIDS_HR.Rdata"))

# Reorder covariates for the forest plot
order_cov <- c("Treatment_discontinuationYes", "Treatment_indicatorZDV only",
               "SexMale", "Age", "Karnofsky_score",
               "Intravenous_drug_useYes", "Prior_opportunistic_infectionsYes", "CD4_baseline")

tbl_df <- tbl_AIDS_HR[match(order_cov, tbl_AIDS_HR$Covariates), ]

RRforest_combined <- cbind.data.frame(
  HR_original = c(NA, tbl_df$HR_original), lower_original = c(NA, tbl_df$IC95_low),  upper_original = c(NA, tbl_df$IC95_high),
  HR_chimera  = c(NA, tbl_df$HR_chimera),  lower_chimera  = c(NA, tbl_df$IC95_low_chimera), upper_chimera = c(NA, tbl_df$IC95_high_chimera),
  HR_synthpop = c(NA, tbl_df$HR_synthpop), lower_synthpop = c(NA, tbl_df$IC95_low_synthpop), upper_synthpop = c(NA, tbl_df$IC95_high_synthpop),
  HR_ctgan    = c(NA, tbl_df$HR_ctgan),    lower_ctgan    = c(NA, tbl_df$IC95_low_ctgan),    upper_ctgan    = c(NA, tbl_df$IC95_high_ctgan)
)

covariates_text <- c(
  "Covariate",
  "Treatment discontinuation (Yes)",
  "Treatment indicator (ZDV only)",
  "Sex (Male)", "Age (10 years)", 
  "Karnofsky score (5 %)",
  "Intravenous drug use (Yes)", 
  "Prior opportunistic infections (Yes)", 
  "CD4 baseline (100 Cells/mm³)"
)

HR_text <- c(
  "Hazard Ratio [95% CI]",
  paste(
    sprintf("%-14s : %s", "REAL",     fmt(tbl_df$HR_original, tbl_df$IC95_low,          tbl_df$IC95_high)),
    sprintf("%-10s : %s", "CHIMERA",  fmt(tbl_df$HR_chimera,  tbl_df$IC95_low_chimera,  tbl_df$IC95_high_chimera)),
    sprintf("%-1s : %s",  "SYNTHPOP", fmt(tbl_df$HR_synthpop, tbl_df$IC95_low_synthpop, tbl_df$IC95_high_synthpop)),
    sprintf("%-12s : %s", "CTGAN",    fmt(tbl_df$HR_ctgan,    tbl_df$IC95_low_ctgan,    tbl_df$IC95_high_ctgan)),
    sep = "\n"
  )
)

SD_text <- c(
  "SDiff",
  paste(
    sprintf("%s %s",   "", ""),
    sprintf("%s %.3f", "", tbl_df$SMD_log1),
    sprintf("%s %.3f", "", tbl_df$SMD_log2),
    sprintf("%s %.3f", "", tbl_df$SMD_log3),
    sep = "\n"
  )
)

tabletext <- cbind.data.frame(covariates_text, HR_text, SD_text)

FP.plot <- forestplot(
  tabletext,
  mean  = cbind(RRforest_combined$HR_original, RRforest_combined$HR_chimera,
                RRforest_combined$HR_synthpop,  RRforest_combined$HR_ctgan),
  lower = cbind(RRforest_combined$lower_original, RRforest_combined$lower_chimera,
                RRforest_combined$lower_synthpop,  RRforest_combined$lower_ctgan),
  upper = cbind(RRforest_combined$upper_original, RRforest_combined$upper_chimera,
                RRforest_combined$upper_synthpop,  RRforest_combined$upper_ctgan),
  zero = 1, Xlog = TRUE, xticks = c(0.75, 1, 1.25, 1.5, 1.75, 2),
  graph.pos = 2, align = c("l", "l", "r"),
  is.summary = c(TRUE, rep(FALSE, nrow(RRforest_combined) - 1)),
  graphwidth = unit(14.5, "cm"), boxsize = 0.1, lwd.ci = 1.5,
  colgap = unit(2, "mm"), line.margin = unit(2.5, "mm"),
  txt_gp = fpTxtGp(
    label   = gpar(fontfamily = "sans", cex = 1, lineheight = 1, fontface = "plain"),
    ticks   = gpar(cex = 1, fontface = "plain"),
    xlab    = gpar(cex = 1, fontface = "plain"),
    summary = gpar(fontface = "plain")
  ),
  xlab = "Hazard Ratio [95% CI]",
  col  = fpColors(box = c("#000000", "#00A087", "#E64B35", "#4DBBD5"),
                  line = c("#000000", "#00A087", "#E64B35", "#4DBBD5"), summary = "#000000")
) |>
  fp_add_lines(h_2 = gpar(col = "#000000", lwd = 1)) |>
  fp_set_zebra_style("#E3EDF5")

# Combine KM + time-dependent AUC + forest plot — AIDS
g_forest <- grid.grabExpr(print(FP.plot))

plot_KM_AUC_forest_AIDS <- (p_km | p_auc) / wrap_elements(full = g_forest) +
  plot_layout(heights = c(0.7, 1), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(family = "sans", face = "bold", size = 16),
        plot.tag.position = c(0.01, 0.99))

ggsave(
  filename = here("figures", "Fig3_KM_AUC_Forest_AIDS.pdf"),
  plot     = plot_KM_AUC_forest_AIDS,
  width    = 300 / 25.4, height = (200 / 25.4) * 2, units = "in", scale = 1
)


# ============================================================
# REIN — Kaplan–Meier, time-dependent AUC, and forest plot
# ============================================================

# ── Kaplan–Meier ────────────────────────────────────────────
load(here("results", "km_rein.Rdata"))

p_km <- ggplot(km_rein, aes(x = time, y = surv, colour = group, linetype = group, fill = group)) +
  geom_step(linewidth = 1) +
  scale_colour_manual(values = npj_palette, name = "Dataset") +
  scale_linetype_manual(values = npj_linetypes, name = "Dataset") +
  scale_fill_manual(values = npj_palette, guide = "none") +
  scale_x_continuous(name = "Follow-up time (days)", limits = c(0, NA),
                     breaks = seq(0, 4000, by = 500), expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(name = "Survival probability", limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.1),
                     labels = scales::number_format(accuracy = 0.1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "") + guides(colour = "none", linetype = "none") + theme_npj() +
  theme(legend.key.width = unit(1.4, "cm"), legend.key.height = unit(0.5, "cm"))


# ── Time-dependent AUC ───────────────────────────────────────
load(here("results", "auc_all_rein.Rdata"))

auc_all <- auc_all_rein |>
  mutate(source = factor(source, levels = c("REAL", "CHIMERA", "SYNTHPOP", "CTGAN")))

auc_means <- auc_all |>
  group_by(source) |>
  summarise(mean_auc = mean(auc, na.rm = TRUE), mean_lower = mean(lower, na.rm = TRUE),
            mean_upper = mean(upper, na.rm = TRUE)) |>
  ungroup() |>
  mutate(source = factor(source, levels = c("REAL", "CHIMERA", "SYNTHPOP", "CTGAN"))) |>
  arrange(source)

p_auc <- ggplot(auc_all, aes(x = time, y = auc, colour = source, linetype = source, fill = source)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = npj_palette, name = NULL) +
  scale_linetype_manual(values = npj_linetypes, name = NULL) +
  scale_fill_manual(values = npj_palette, guide = "none") +
  scale_x_continuous(name = "Follow-up time (days)", breaks = seq(0, 4000, by = 500),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(name = "Time-dependent AUC", limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.2), labels = number_format(accuracy = 0.1),
                     expand = expansion(mult = c(0, 0.02))) +
  annotate("text", x = 1200, y = 0.295, label = "REAL",     hjust = 0, size = 4, colour = npj_palette["REAL"]) +
  annotate("text", x = 2000, y = 0.295, label = "iAUC: 0.694 [95% CI: 0.692\u20130.696]", hjust = 0, size = 4) +
  annotate("text", x = 1200, y = 0.262, label = "CHIMERA",  hjust = 0, size = 4, colour = npj_palette["CHIMERA"]) +
  annotate("text", x = 2000, y = 0.262, label = "iAUC: 0.662 [95% CI: 0.660\u20130.664]", hjust = 0, size = 4) +
  annotate("text", x = 1200, y = 0.229, label = "SYNTHPOP", hjust = 0, size = 4, colour = npj_palette["SYNTHPOP"]) +
  annotate("text", x = 2000, y = 0.229, label = "iAUC: 0.679 [95% CI: 0.677\u20130.681]", hjust = 0, size = 4) +
  annotate("text", x = 1200, y = 0.196, label = "CTGAN",    hjust = 0, size = 4, colour = npj_palette["CTGAN"]) +
  annotate("text", x = 2000, y = 0.196, label = "iAUC: 0.688 [95% CI: 0.686\u20130.689]", hjust = 0, size = 4) +
  labs(title = "") + theme_npj() + theme(legend.position = "none")


# ── Forest plot — hazard ratios from Cox model (REIN) ────────
load(here("results", "tbl_REIN_HR.Rdata"))
tbl_df <- tbl_REIN_HR[1:9, ]

RRforest_combined <- cbind.data.frame(
  HR_original = c(NA, tbl_df$HR_original), lower_original = c(NA, tbl_df$IC95_low),  upper_original = c(NA, tbl_df$IC95_high),
  HR_chimera  = c(NA, tbl_df$HR_chimera),  lower_chimera  = c(NA, tbl_df$IC95_low_chimera), upper_chimera = c(NA, tbl_df$IC95_high_chimera),
  HR_synthpop = c(NA, tbl_df$HR_synthpop), lower_synthpop = c(NA, tbl_df$IC95_low_synthpop), upper_synthpop = c(NA, tbl_df$IC95_high_synthpop),
  HR_ctgan    = c(NA, tbl_df$HR_ctgan),    lower_ctgan    = c(NA, tbl_df$IC95_low_ctgan),    upper_ctgan    = c(NA, tbl_df$IC95_high_ctgan)
)

covariates_text <- c(
  "Covariate",
  "Age (10 years)", 
  "Body mass index (kg/m\u00b2)",
  "Cardiac arrhythmia (Yes)", 
  "Chronic respiratory failure (Yes)",
  "Cirrhosis (Yes)", 
  "Coronary artery disease (Yes)",
  "Diabetes (Yes)", 
  "Heart failure (NYHA I\u2013II)", 
  "Heart failure (NYHA III\u2013IV)"
)

HR_text <- c(
  "Hazard Ratio [95% CI]",
  paste(
    sprintf("%-14s : %s", "REAL",     fmt(tbl_df$HR_original, tbl_df$IC95_low,          tbl_df$IC95_high)),
    sprintf("%-10s : %s", "CHIMERA",  fmt(tbl_df$HR_chimera,  tbl_df$IC95_low_chimera,  tbl_df$IC95_high_chimera)),
    sprintf("%-1s : %s",  "SYNTHPOP", fmt(tbl_df$HR_synthpop, tbl_df$IC95_low_synthpop, tbl_df$IC95_high_synthpop)),
    sprintf("%-12s : %s", "CTGAN",    fmt(tbl_df$HR_ctgan,    tbl_df$IC95_low_ctgan,    tbl_df$IC95_high_ctgan)),
    sep = "\n"
  )
)

SD_text <- c(
  "SDiff",
  paste(
    sprintf("%s %s",   "", ""),
    sprintf("%s %.3f", "", tbl_df$SMD_log1),
    sprintf("%s %.3f", "", tbl_df$SMD_log2),
    sprintf("%s %.3f", "", tbl_df$SMD_log3),
    sep = "\n"
  )
)

tabletext <- cbind.data.frame(covariates_text, HR_text, SD_text)

FP.plot <- forestplot(
  tabletext,
  mean  = cbind(RRforest_combined$HR_original, RRforest_combined$HR_chimera,
                RRforest_combined$HR_synthpop,  RRforest_combined$HR_ctgan),
  lower = cbind(RRforest_combined$lower_original, RRforest_combined$lower_chimera,
                RRforest_combined$lower_synthpop,  RRforest_combined$lower_ctgan),
  upper = cbind(RRforest_combined$upper_original, RRforest_combined$upper_chimera,
                RRforest_combined$upper_synthpop,  RRforest_combined$upper_ctgan),
  zero = 1, Xlog = TRUE, xticks = c(0.75, 1, 1.25, 1.5, 1.75, 2),
  graph.pos = 2, align = c("l", "l", "r"),
  is.summary = c(TRUE, rep(FALSE, nrow(RRforest_combined) - 1)),
  graphwidth = unit(14.5, "cm"), boxsize = 0.1, lwd.ci = 1.5,
  colgap = unit(2, "mm"), line.margin = unit(2.5, "mm"),
  txt_gp = fpTxtGp(
    label   = gpar(fontfamily = "sans", cex = 1, lineheight = 1, fontface = "plain"),
    ticks   = gpar(cex = 1, fontface = "plain"),
    xlab    = gpar(cex = 1, fontface = "plain"),
    summary = gpar(fontface = "plain")
  ),
  xlab = "Hazard Ratio [95% CI]",
  col  = fpColors(box = c("#000000", "#00A087", "#E64B35", "#4DBBD5"),
                  line = c("#000000", "#00A087", "#E64B35", "#4DBBD5"), summary = "#000000")
) |>
  fp_add_lines(h_2 = gpar(col = "#000000", lwd = 1)) |>
  fp_set_zebra_style("#E3EDF5")

# Combine KM + time-dependent AUC + forest plot — REIN
g_forest <- grid.grabExpr(print(FP.plot))

plot_KM_AUC_forest_REIN <- (p_km | p_auc) / wrap_elements(full = g_forest) +
  plot_layout(heights = c(0.7, 1), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(family = "sans", face = "bold", size = 16),
        plot.tag.position = c(0.01, 0.99))

ggsave(
  filename = here( "figures", "Fig4_KM_AUC_Forest_REIN.pdf"),
  plot     = plot_KM_AUC_forest_REIN,
  width    = 300 / 25.4, height = (200 / 25.4) * 2, units = "in", scale = 1
)


# ============================================================
# REIN — Clinical score: selection performance boxplots
# ============================================================

# Load selection metrics (sensitivity, specificity, kappa) and final-model AUC
load(here("results", "chimera_selection_metrics_df.Rdata"))
load(here("results", "synthpop_selection_metrics_df.Rdata"))
load(here("results", "ctgan_selection_metrics_df.Rdata"))
load(here("results", "auc_final_model_rein_syn.Rdata"))

# Append final-model AUC as metric T4
chimera_selection_metrics_df$T4  <- auc_final_model_rein_syn |> filter(source == "CHIMERA")  |> pull(auc)
synthpop_selection_metrics_df$T4 <- auc_final_model_rein_syn |> filter(source == "SYNTHPOP") |> pull(auc)
ctgan_selection_metrics_df$T4    <- auc_final_model_rein_syn |> filter(source == "CTGAN")    |> pull(auc)


# ── Palette and aesthetic settings ───────────────────────────
npj_palette <- c("CHIMERA" = "#00A087", "SYNTHPOP" = "#E64B35", "CTGAN" = "#4DBBD5")
best_point_fill   <- "#F2C14E"
best_point_colour <- "#2B2B2B"
ref_line_colour   <- "#8C8C8C"
grid_colour       <- "#D9E0E6"
strip_fill        <- "#F3F6F8"
border_colour     <- "#D0D7DE"
text_colour       <- "#1F2933"
SYNTH_ORDER       <- c("CHIMERA", "SYNTHPOP", "CTGAN")

# ── Metric metadata for score selection ──────────────────────
metric_info <- tibble(
  metric       = paste0("T", 1:4),
  metric_label = factor(c("Sensitivity", "Specificity", "Cohen's kappa", "AUC"),
                        levels = c("Sensitivity", "Specificity", "Cohen's kappa", "AUC")),
  metric_group = "Selection",
  target       = c(1, 1, 1, 1)   # All metrics are maximized
)

ref_lines <- metric_info

# ── Data preparation ──────────────────────────────────────────
prepare_one_synth <- function(df, dataset_name, synthesizer_name) {
  df |>
    mutate(syn_id = row_number(), replicate_id = row_number(),
           dataset_name = dataset_name, synthesizer = toupper(synthesizer_name)) |>
    select(dataset_name, synthesizer, syn_id, replicate_id, starts_with("T")) |>
    pivot_longer(cols = starts_with("T"), names_to = "metric", values_to = "value") |>
    left_join(metric_info, by = "metric") |>
    mutate(
      synthesizer  = factor(synthesizer, levels = SYNTH_ORDER),
      metric_label = factor(metric_label, levels = c("Sensitivity", "Specificity", "Cohen's kappa", "AUC")),
      metric_group = metric_group
    ) |>
    filter(is.finite(value))
}

# Select the best run based on minimum Euclidean distance to the ideal vector (1, 1, 1, 1)
get_best_syn <- function(df, dataset_name, synthesizer_name) {
  df |>
    mutate(syn_id = row_number(), dataset_name = dataset_name, synthesizer = toupper(synthesizer_name)) |>
    slice_min(order_by = distance, n = 1, with_ties = FALSE) |>
    transmute(dataset_name, synthesizer, best_syn_id = syn_id, distance)
}

best_of_50 <- bind_rows(
  get_best_syn(chimera_selection_metrics_df,  "REIN", "CHIMERA"),
  get_best_syn(synthpop_selection_metrics_df, "REIN", "SYNTHPOP"),
  get_best_syn(ctgan_selection_metrics_df,    "REIN", "CTGAN")
)
print(best_of_50)

plot_df <- bind_rows(
  prepare_one_synth(chimera_selection_metrics_df,  "REIN", "CHIMERA"),
  prepare_one_synth(synthpop_selection_metrics_df, "REIN", "SYNTHPOP"),
  prepare_one_synth(ctgan_selection_metrics_df,    "REIN", "CTGAN")
) |>
  left_join(best_of_50 |> select(dataset_name, synthesizer, best_syn_id),
            by = c("dataset_name", "synthesizer")) |>
  mutate(is_best = syn_id == best_syn_id, synthesizer = factor(synthesizer, levels = SYNTH_ORDER))

# ── Theme ─────────────────────────────────────────────────────
theme_npj_box <- function(base_size = 11, base_family = "Helvetica") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      text          = element_text(colour = text_colour),
      plot.title    = element_text(face = "bold", size = rel(1.2), hjust = 0),
      plot.subtitle = element_text(size = rel(0.95), hjust = 0, margin = margin(b = 8)),
      axis.title.x  = element_blank(),
      axis.title.y  = element_text(size = rel(1.0), margin = margin(r = 8)),
      axis.text.x   = element_text(face = "bold", size = rel(0.9), colour = text_colour),
      axis.text.y   = element_text(size = rel(0.85), colour = text_colour),
      panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.3),
      strip.background   = element_rect(fill = strip_fill, colour = NA),
      strip.text         = element_text(face = "bold", size = rel(0.9), colour = text_colour),
      panel.border       = element_rect(fill = NA, colour = border_colour, linewidth = 0.5),
      panel.spacing      = unit(0.8, "lines"),
      legend.position    = "none",
      plot.margin        = margin(8, 10, 8, 10)
    )
}

# ── Plot ──────────────────────────────────────────────────────
make_metric_row <- function(data, dataset_selected, y_lab) {
  df_plot <- data |> filter(dataset_name == dataset_selected)
  ggplot(df_plot, aes(x = synthesizer, y = value, fill = synthesizer)) +
    geom_hline(data = ref_lines, aes(yintercept = target), inherit.aes = FALSE,
               linetype = "dashed", colour = ref_line_colour, linewidth = 0.4) +
    geom_boxplot(width = 0.6, outlier.shape = NA, colour = "#3A3A3A") +
    geom_point(data = df_plot |> filter(is_best), shape = 23, size = 2.8, stroke = 0.4,
               fill = best_point_fill, colour = best_point_colour) +
    facet_wrap(~ metric_label, ncol = 2) +
    scale_fill_manual(values = npj_palette) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    labs(y = y_lab) + coord_cartesian(clip = "off") + theme_npj_box()
}

p_REIN <- make_metric_row(plot_df, "REIN", "")
p_REIN

ggsave(filename=here("figures", "FigS8_ScoreDistrBestofM_REIN.pdf"),
       p_REIN, width = 300 / 25.4, height = (200 / 25.4) / 1.5, device = cairo_pdf, bg = "white")

ggsave(filename=here("figures", "FigS8_ScoreDistrBestofM_REIN.png"),
       p_REIN, width = 300 / 25.4, height = (200 / 25.4) / 1.5, dpi = 350, bg = "white")


# ============================================================
# REIN — Final model: ROC curve and calibration curves
# ============================================================

npj_palette <- c(
  "REAL"     = "#000000",
  "CHIMERA"  = "#00A087",
  "SYNTHPOP" = "#E64B35",
  "CTGAN"    = "#4DBBD5"
)

# ── ROC curve — REIN final model ─────────────────────────────
load(here("results", "rein.roc.original.Rdata"))
load(here("results", "rein.roc.chimera.Rdata"))
load(here("results", "rein.roc.synthpop.Rdata"))
load(here("results", "rein.roc.ctgan.Rdata"))

roc_list <- list(
  REAL     = rein.roc.original,
  CHIMERA  = rein.roc.chimera,
  SYNTHPOP = rein.roc.synthpop,
  CTGAN    = rein.roc.ctgan
)

df_roc <- lapply(names(roc_list), function(nm) {
  r <- roc_list[[nm]]
  tibble(Dataset = nm, Specificity = r$specificities, Sensitivity = r$sensitivities)
}) |>
  bind_rows() |>
  mutate(Dataset = factor(Dataset, levels = c("REAL", "CHIMERA", "SYNTHPOP", "CTGAN")),
         FPR = 1 - Specificity)

plt_roc <- ggplot(df_roc, aes(x = FPR, y = Sensitivity, colour = Dataset, linetype = Dataset)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.4, colour = "#000000") +
  geom_line(linewidth = 0.65) +
  scale_colour_manual(values = npj_palette) +
  scale_linetype_manual(values = npj_linetypes) +
  scale_x_continuous(name = "1 \u2212 Specificity", limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     labels = scales::number_format(accuracy = 0.1), expand = expansion(mult = 0.01)) +
  scale_y_continuous(name = "Sensitivity", limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     labels = scales::number_format(accuracy = 0.1), expand = expansion(mult = 0.01)) +
  labs(title = "") + theme_npj() + theme(legend.position = "none") +
  annotate("text", x = 0.29, y = 0.262, label = "REAL",     hjust = 0, size = 4, colour = npj_palette["REAL"]) +
  annotate("text", x = 0.49, y = 0.262, label = "AUC: 0.741 [95% CI: 0.724\u20130.757]", hjust = 0, size = 4) +
  annotate("text", x = 0.29, y = 0.229, label = "CHIMERA",  hjust = 0, size = 4, colour = npj_palette["CHIMERA"]) +
  annotate("text", x = 0.49, y = 0.229, label = "AUC: 0.720 [95% CI: 0.703\u20130.738]", hjust = 0, size = 4) +
  annotate("text", x = 0.29, y = 0.196, label = "SYNTHPOP", hjust = 0, size = 4, colour = npj_palette["SYNTHPOP"]) +
  annotate("text", x = 0.49, y = 0.196, label = "AUC: 0.734 [95% CI: 0.717\u20130.751]", hjust = 0, size = 4) +
  annotate("text", x = 0.29, y = 0.163, label = "CTGAN",    hjust = 0, size = 4, colour = npj_palette["CTGAN"]) +
  annotate("text", x = 0.49, y = 0.163, label = "AUC: 0.718 [95% CI: 0.701\u20130.736]", hjust = 0, size = 4)


# ── Calibration curves ────────────────────────────────────────
load(here("results", "res.original.cal.Rdata"))
load(here("results", "res.chimera.cal.Rdata"))
load(here("results", "res.synthpop.cal.Rdata"))
load(here("results", "res.ctgan.cal.Rdata"))

# Reference: calibration curve from real data (imputation = 0)
res.original.cal.curve <- res.original.cal$calibration_curve
res.original.cal.curve$imputation <- 0

# Helper to build a calibration plot for one synthesizer
# - Shows the best synthetic run (color) and the real data curve (black)
# - Includes 95% CI ribbons for both
make_calib_plot <- function(synth_curves, best_id, colour_hex, group_label) {
  df <- rbind.data.frame(res.original.cal.curve, synth_curves) |>
    mutate(groupe = case_when(
      imputation == 0       ~ "REAL",
      imputation == best_id ~ paste0("Best ", group_label),
      TRUE                  ~ paste0("Other ", group_label)
    )) |>
    filter(groupe %in% c("REAL", paste0("Best ", group_label)))
  
  ggplot(df, aes(x = pHat, y = obs, group = interaction(imputation, groupe))) +
    geom_ribbon(data = subset(df, groupe != "REAL"),
                aes(ymin = ymin, ymax = ymax), fill = colour_hex, alpha = 0.10) +
    geom_line(data  = subset(df, groupe != "REAL"), color = colour_hex, linewidth = 1) +
    geom_ribbon(data = subset(df, groupe == "REAL"),
                aes(ymin = ymin, ymax = ymax), fill = "black", alpha = 0.10) +
    geom_line(data  = subset(df, groupe == "REAL"), color = "black", linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    coord_cartesian(xlim = c(0, 0.6), ylim = c(0, 1)) +
    labs(x = "Predicted probability", y = "Observed probability") +
    theme_npj(base_size = 12)
}

plt_chimera  <- make_calib_plot(res.chimera.cal$calibration_curves,  4,  "#00A087", "chimera")
plt_synthpop <- make_calib_plot(res.synthpop.cal$calibration_curves, 33, "#E64B35", "synthpop")
plt_ctgan    <- make_calib_plot(res.ctgan.cal$calibration_curves,    32, "#4DBBD5", "ctgan")

# Combine ROC + calibration panels (main figure)
plot_ROCFinalModel_Cal_REIN <- (plt_roc | plt_chimera) / (plt_synthpop | plt_ctgan) +
  plot_layout(heights = c(0.7, 0.7), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(family = "sans", face = "bold", size = 16),
        plot.tag.position = c(0.01, 0.99))

ggsave(
  filename = here("figures", "Fig5_AUC_Calib_REIN.pdf"),
  plot     = plot_ROCFinalModel_Cal_REIN,
  width    = 300 / 25.4, height = (150 / 25.4) * 2, units = "in", scale = 1
)


# ── Supplementary: all calibration curves (best + other runs) ─
make_calib_plot_all <- function(synth_curves, best_id, colour_hex, group_label) {
  df <- rbind.data.frame(res.original.cal.curve, synth_curves) |>
    mutate(groupe = case_when(
      imputation == 0       ~ "REAL",
      imputation == best_id ~ paste0("Best ", group_label),
      TRUE                  ~ paste0("Other ", group_label)
    ))
  
  ggplot(df, aes(x = pHat, y = obs, group = interaction(imputation, groupe))) +
    # All other synthetic runs in light grey
    geom_line(data = subset(df, groupe == paste0("Other ", group_label)),
              color = "grey80", linewidth = 0.3) +
    # Best synthetic run
    geom_line(data = subset(df, groupe == paste0("Best ", group_label)),
              color = colour_hex, linewidth = 1) +
    # Real data curve
    geom_line(data = subset(df, groupe == "REAL"), color = "black", linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    coord_cartesian(xlim = c(0, 0.6), ylim = c(0, 1)) +
    labs(x = "Predicted probability", y = "Observed probability") +
    theme_npj(base_size = 12)
}

plt_chimera  <- make_calib_plot_all(res.chimera.cal$calibration_curves,  4,  "#00A087", "chimera")
plt_synthpop <- make_calib_plot_all(res.synthpop.cal$calibration_curves, 33, "#E64B35", "synthpop")
plt_ctgan    <- make_calib_plot_all(res.ctgan.cal$calibration_curves,    32, "#4DBBD5", "ctgan")

plot_ROCFinalModel_Cal_REIN_supp <- (plt_chimera | plt_synthpop | plt_ctgan) +
  plot_layout(heights = c(0.5), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(family = "sans", face = "bold", size = 16),
        plot.tag.position = c(0.01, 0.99))

ggsave(
  filename = here("figures", "FigS9_AUC_Calib_REIN.pdf"),
  plot     = plot_ROCFinalModel_Cal_REIN_supp,
  width    = 300 / 25.4, height = 130 / 25.4, units = "in", scale = 1
)
