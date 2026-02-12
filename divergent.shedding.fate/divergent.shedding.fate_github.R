library(tidyverse)    
library(patchwork)
library(ggridges)
library(ComplexUpset)
library(pheatmap)
library(corrplot)
library(entropy)
library(kableExtra)
library(stringdist)
library(purrr)
library(readxl)
library(readr)
library(reshape2)
library(grid)
library(gridExtra)

##### Fig 1: Kidney late-shed vs late non-shed classification + panels A-C

library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)


##### Fig 1 panel
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(patchwork)

urine_bc_levels <- read_tsv("urine_bc_levels.tsv", show_col_types = FALSE)
tissue_bc_levels <- read_tsv("tissue_bc_levels.tsv", show_col_types = FALSE)
cutoff_99pct_stock_barcodes <- read_tsv("cutoff_99pct_stock_barcodes.tsv", show_col_types = FALSE)

SAVE_DIR <- "../plots/urine_tissue_manuscript5"


LAST_K <- 2

if (!exists("is_kidney")) {
  is_kidney <- function(x) grepl("kidney", x, ignore.case = TRUE)
}

LABEL_SHED    <- "Late-shed"
LABEL_NONSHED <- "Late non-shed"
status_levels <- c(LABEL_SHED, LABEL_NONSHED)

scale_def_colors <- scale_fill_manual(values = c(
  "Late-shed" = "#F8766D",      # ggplot default coral
  "Late non-shed" = "#00BFC4"   # ggplot default teal
))


last_k_days <- urine_bc_levels %>%
  distinct(animal, Days_pi) %>%
  group_by(animal) %>%
  arrange(Days_pi, .by_group = TRUE) %>%
  slice_tail(n = LAST_K) %>%
  mutate(in_lastK = TRUE) %>%
  ungroup()

present_lastK <- urine_bc_levels %>%
  inner_join(last_k_days, by = c("animal", "Days_pi")) %>%
  group_by(animal, barcode) %>%
  summarise(present_lastK = any(bc_level > 0, na.rm = TRUE), .groups = "drop") %>%
  filter(present_lastK) %>%
  select(animal, barcode)

kidney_barcodes <- tissue_bc_levels %>%
  filter(is_kidney(organ)) %>%
  group_by(animal, barcode) %>%
  summarise(kidney_present = any(bc_level > 0, na.rm = TRUE), .groups = "drop") %>%
  filter(kidney_present) %>%
  select(animal, barcode)

kidney_classified <- kidney_barcodes %>%
  left_join(present_lastK %>% mutate(in_lastK = TRUE), by = c("animal", "barcode")) %>%
  mutate(
    in_lastK = if_else(is.na(in_lastK), FALSE, in_lastK),
    status   = if_else(in_lastK, LABEL_SHED, LABEL_NONSHED),
    status   = factor(status, levels = status_levels)
  )


panel_title_theme <- theme(
  plot.title = element_text(face = "plain", size = 10, hjust = 0.5),
  plot.title.position = "plot"
)


kidney_counts <- kidney_classified %>%
  count(animal, status, name = "n_barcodes")

p_A <- ggplot(kidney_counts, aes(x = animal, y = n_barcodes, fill = status)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  labs(
    title = "Kidney barcode counts by late urine detection",
    x = "Mouse",
    y = "Number of unique kidney barcodes"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  ) +
  panel_title_theme


kidney_total_by_status <- tissue_bc_levels %>%
  filter(is_kidney(organ)) %>%
  group_by(animal, barcode) %>%
  summarise(kidney_bc_level = sum(bc_level, na.rm = TRUE), .groups = "drop") %>%
  inner_join(kidney_classified %>% select(animal, barcode, status),
             by = c("animal", "barcode")) %>%
  group_by(animal, status) %>%
  summarise(total_bc_level = sum(kidney_bc_level, na.rm = TRUE), .groups = "drop") %>%
  mutate(status = factor(status, levels = status_levels))

p_B <- ggplot(kidney_total_by_status,
              aes(x = animal, y = total_bc_level, fill = status)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  scale_y_log10() +
  labs(
    title = "Total kidney barcode level by late urine detection",
    x = "Mouse",
    y = "Total kidney barcode level (log10)",
    fill = ""   # the ONLY color legend (blank title)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  ) +
  panel_title_theme


kidney_bc_levels_labeled <- tissue_bc_levels %>%
  filter(is_kidney(organ)) %>%
  group_by(animal, barcode) %>%
  summarise(bc_level = sum(bc_level, na.rm = TRUE), .groups = "drop") %>%
  inner_join(kidney_classified, by = c("animal", "barcode")) %>%
  filter(bc_level > 0)

p_C <- ggplot(kidney_bc_levels_labeled,
              aes(x = status, y = bc_level, fill = status)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9) +
  facet_wrap(~ animal, scales = "free_y") +
  scale_y_log10() +
  labs(
    title = "Kidney barcode level distributions by late urine detection",
    x = NULL,
    y = "Kidney barcode level (log10)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 25, hjust = 1)
  ) +
  panel_title_theme


# fig1_panel <- (p_A | p_B) / p_C +
#   plot_layout(heights = c(1, 1.15)) +
#   plot_annotation(tag_levels = "A") &
#   theme(
#     plot.tag = element_text(face = "bold", size = 12),
#     plot.tag.position = c(0.0, 1.02)
#   )

fig1_panel <- (p_A | p_B) / p_C +
  plot_layout(heights = c(1, 1.15)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 12),
    plot.tag.position = c(0.0, 1.02),
    legend.position = "top",
    legend.justification = "center",
    legend.box.just = "center"
  )&
  guides(fill = guide_legend(title = NULL))


print(fig1_panel)

ggsave(
  filename = "fig1_panel.pdf",
  plot     = fig1_panel,
  path     = SAVE_DIR,
  width    = 7.24, height = 8
)








##### Fig 2 panel: Donuts (A,B) on row 1; Sum overlays (C,D) on row 2

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(purrr)
  library(scales)
  library(grid)
  library(gridExtra)
  library(patchwork)
})



SAVE_DIR <- "../plots/urine_tissue_manuscript5"

TOP_N  <- 50
LAST_K <- 2
REQUIRE_ALL_LASTK <- FALSE  # TRUE = present in *all* of the last K; FALSE = present in *any* of the last K

is_kidney <- function(x) grepl("kidney", x, ignore.case = TRUE)

# If animal_colors is not defined upstream, fall back safely:
if (!exists("animal_colors")) {
  animals <- c("FL","FR","ML","MR")
} else {
  animals <- names(animal_colors)
}


plot_diversity_donut_w_colored_barcodes <- function(bc_table,
                                                    barcode_color_highlights,
                                                    faceting_variable) {
  
  valid_highlights <- barcode_color_highlights[barcode_color_highlights %in% bc_table$barcode]
  if (length(valid_highlights) == 0) {
    stop("None of the barcode highlights are present in the dataset.")
  }
  
  winner_and_bg_colors <- set_names(
    c(scales::hue_pal()(length(valid_highlights)), 'gray75', 'gray60'),
    nm = c(valid_highlights, 'background1', 'background2')
  )
  
  bc_table <- bc_table %>%
    group_by(sample) %>%
    arrange(desc(frac_weighted_count), .by_group = TRUE) %>%
    mutate(
      ymax = cumsum(frac_weighted_count),
      ymin = c(0, head(ymax, n = -1)),
      bc_color = ifelse(barcode %in% valid_highlights,
                        barcode,
                        ifelse(row_number() %% 2 == 0, 'background1', 'background2'))
    ) %>%
    ungroup()
  
  ggplot(bc_table, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3.2, fill = bc_color)) +
    scale_fill_manual(values = winner_and_bg_colors, labels = names(winner_and_bg_colors)) +
    geom_rect(show.legend = FALSE) +
    coord_polar(theta = 'y') +
    xlim(c(2, 4)) +
    facet_wrap(vars({{faceting_variable}})) +
    theme_classic() +
    theme(
      line = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    ) +
    guides(fill = guide_legend('Barcode'))
}


wrap_with_title <- function(grob_obj, title_text) {
  gridExtra::arrangeGrob(
    grob_obj,
    top = grid::textGrob(title_text, gp = grid::gpar(fontsize = 10), just = "center")
  )
}


last_k_days <- urine_bc_levels %>%
  distinct(animal, Days_pi) %>%
  group_by(animal) %>%
  arrange(Days_pi, .by_group = TRUE) %>%
  slice_tail(n = LAST_K) %>%
  mutate(in_lastK = TRUE) %>%
  ungroup()

present_lastK <- urine_bc_levels %>%
  inner_join(last_k_days, by = c("animal","Days_pi")) %>%
  group_by(animal, barcode) %>%
  summarise(
    present_lastK = if (REQUIRE_ALL_LASTK) all(bc_level > 0, na.rm = TRUE) else any(bc_level > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(present_lastK)

kidney_rank <- tissue_bc_levels %>%
  filter(is_kidney(organ)) %>%
  group_by(animal, barcode) %>%
  summarise(kidney_load = sum(bc_level, na.rm = TRUE), .groups = "drop") %>%
  inner_join(present_lastK, by = c("animal","barcode")) %>%
  group_by(animal) %>%
  arrange(desc(kidney_load), .by_group = TRUE) %>%
  mutate(kidney_rank = row_number()) %>%
  slice_head(n = TOP_N) %>%
  ungroup()

winner_barcodes_by_animal <- kidney_rank %>%
  group_by(animal) %>%
  summarise(barcodes = list(barcode), .groups = "drop") %>%
  deframe()

urine_bc_kidTop <- urine_bc_levels %>%
  semi_join(kidney_rank, by = c("animal","barcode"))


donut_plot_list <- lapply(animals, function(a) {
  bc_table <- urine_bc_levels %>% filter(animal == a)
  highlights <- winner_barcodes_by_animal[[a]]
  if (is.null(highlights) || length(highlights) == 0) {
    return(ggplot() + theme_void() + labs(title = a))
  }
  plot_diversity_donut_w_colored_barcodes(
    bc_table = bc_table,
    barcode_color_highlights = highlights,
    faceting_variable = Days_pi
  ) + labs(title = a)
})

donut_combined <- gridExtra::arrangeGrob(grobs = donut_plot_list, ncol = 2)


sum_selected <- urine_bc_kidTop %>%
  group_by(animal, Days_pi, mean_ul_urine) %>%
  summarise(sum_selected_barcodes = sum(bc_level, na.rm = TRUE), .groups = "drop")

sum_overlay_plots <- sum_selected %>%
  split(.$animal) %>%
  imap(~ ggplot(.x, aes(x = Days_pi)) +
         geom_area(aes(y = mean_ul_urine), fill = "grey85") +
         geom_area(aes(y = sum_selected_barcodes), fill = "navajowhite3") +
         labs(title = .y, x = "", y = "") +
         theme_minimal() +
         theme(legend.position = "none"))

overlay_combined <- gridExtra::arrangeGrob(
  sum_overlay_plots[["FL"]],
  sum_overlay_plots[["FR"]],
  sum_overlay_plots[["ML"]],
  sum_overlay_plots[["MR"]],
  nrow = 4, ncol = 1,
  bottom = grid::textGrob("Days Post-Injection", gp = grid::gpar(fontsize = 11)),
  left   = grid::textGrob("Barcode level", rot = 90, gp = grid::gpar(fontsize = 11))
)


absent_lastK <- urine_bc_levels %>%
  distinct(animal, barcode) %>%
  anti_join(present_lastK, by = c("animal","barcode"))

kidney_rank_notLastK <- tissue_bc_levels %>%
  filter(is_kidney(organ)) %>%
  group_by(animal, barcode) %>%
  summarise(kidney_load = sum(bc_level, na.rm = TRUE), .groups = "drop") %>%
  inner_join(absent_lastK, by = c("animal","barcode")) %>%
  group_by(animal) %>%
  arrange(desc(kidney_load), .by_group = TRUE) %>%
  mutate(kidney_rank = row_number()) %>%
  slice_head(n = TOP_N) %>%
  ungroup()

winner_barcodes_by_animal_notLastK <- kidney_rank_notLastK %>%
  group_by(animal) %>%
  summarise(barcodes = list(barcode), .groups = "drop") %>%
  deframe()

urine_bc_kidTop_notLastK <- urine_bc_levels %>%
  semi_join(kidney_rank_notLastK, by = c("animal","barcode"))

# --- Donuts: late non-shed (B)
donut_plot_list_notLastK <- lapply(animals, function(a) {
  bc_table <- urine_bc_levels %>% filter(animal == a)
  highlights <- winner_barcodes_by_animal_notLastK[[a]]
  
  if (is.null(highlights) || length(highlights) == 0) {
    return(ggplot() + theme_void() + labs(title = a))
  }
  
  highlights_in_urine <- intersect(highlights, bc_table$barcode)
  if (length(highlights_in_urine) == 0) {
    return(ggplot() + theme_void() + labs(title = a))
  }
  
  plot_diversity_donut_w_colored_barcodes(
    bc_table = bc_table,
    barcode_color_highlights = highlights_in_urine,
    faceting_variable = Days_pi
  ) + labs(title = a)
})

donut_combined_notLastK <- gridExtra::arrangeGrob(grobs = donut_plot_list_notLastK, ncol = 2)

sum_selected_notLastK <- urine_bc_kidTop_notLastK %>%
  group_by(animal, Days_pi, mean_ul_urine) %>%
  summarise(sum_selected_barcodes = sum(bc_level, na.rm = TRUE), .groups = "drop")

sum_overlay_plots_notLastK <- sum_selected_notLastK %>%
  split(.$animal) %>%
  imap(~ ggplot(.x, aes(x = Days_pi)) +
         geom_area(aes(y = mean_ul_urine), fill = "grey85") +
         geom_area(aes(y = sum_selected_barcodes), fill = "navajowhite3") +
         labs(title = .y, x = "", y = "") +
         theme_minimal() +
         theme(legend.position = "none"))

overlay_combined_notLastK <- gridExtra::arrangeGrob(
  sum_overlay_plots_notLastK[["FL"]],
  sum_overlay_plots_notLastK[["FR"]],
  sum_overlay_plots_notLastK[["ML"]],
  sum_overlay_plots_notLastK[["MR"]],
  nrow = 4, ncol = 1,
  bottom = grid::textGrob("Days Post-Injection", gp = grid::gpar(fontsize = 11)),
  left   = grid::textGrob("Barcode level", rot = 90, gp = grid::gpar(fontsize = 11))
)


A_grob <- wrap_with_title(donut_combined,
                          "Donut plots: top kidney barcodes (Late-shed)")
B_grob <- wrap_with_title(donut_combined_notLastK,
                          "Donut plots: top kidney barcodes (Late non-shed)")
C_grob <- wrap_with_title(overlay_combined,
                          "Sum-overlay: urine vs selected barcodes (Late-shed)")
D_grob <- wrap_with_title(overlay_combined_notLastK,
                          "Sum-overlay: urine vs selected barcodes (Late non-shed)")


pA <- patchwork::wrap_elements(full = A_grob)
pB <- patchwork::wrap_elements(full = B_grob)
pC <- patchwork::wrap_elements(full = C_grob)
pD <- patchwork::wrap_elements(full = D_grob)


fig2_panel <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1, 1), widths = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 12),
    plot.tag.position = c(0.0, 1.02)
  )

print(fig2_panel)

ggsave(
  filename = "fig2_panel.pdf",
  plot     = fig2_panel,
  path     = SAVE_DIR,
  width    = 16, height = 32
)





##### Fig 3: clone-resolved urine trajectories (cap=25), stacked A/B/C (one column)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(purrr)
  library(grid)
  library(gridExtra)
  library(scales)
  library(tibble)
})


TOP_N    <- 50
CAP_VAL  <- 25
LAST_K   <- 2
REQUIRE_ALL_LASTK <- FALSE

SAVE_DIR <- "../plots/urine_tissue_manuscript5"
dir.create(SAVE_DIR, showWarnings = FALSE, recursive = TRUE)

TOP_LABEL <- paste0("Top-", TOP_N)


is_kidney <- function(x) grepl("kidney|renal", x, ignore.case = TRUE)

random_topN_control <- function(urine_tbl, selected_tbl, n = TOP_N, seed = 1) {
  set.seed(seed)
  avail <- urine_tbl %>% distinct(animal, barcode)
  pool  <- anti_join(avail,
                     selected_tbl %>% distinct(animal, barcode),
                     by = c("animal","barcode"))
  pool %>%
    group_by(animal) %>%
    mutate(.rand = runif(n())) %>%
    arrange(.rand, .by_group = TRUE) %>%
    mutate(.rank = row_number()) %>%
    filter(.rank <= n) %>%
    ungroup() %>%
    select(animal, barcode)
}

plot_barcode_lines_capped_for_mouse <- function(urine_tbl, selected_tbl, mouse_id,
                                                title_suffix = "", cap = 25) {
  sel_bcs <- selected_tbl %>%
    filter(animal == mouse_id) %>%
    pull(barcode) %>% unique()
  
  if (length(sel_bcs) == 0) {
    return(ggplot() + theme_void() + labs(title = paste(mouse_id, "(no data)")))
  }
  
  bc_ts <- urine_tbl %>%
    semi_join(tibble(animal = mouse_id, barcode = sel_bcs),
              by = c("animal","barcode")) %>%
    group_by(animal, barcode, Days_pi) %>%
    summarise(bc_day_level = sum(bc_level, na.rm = TRUE), .groups = "drop") %>%
    mutate(bc_capped = pmin(bc_day_level, cap))
  
  pal <- setNames(scales::hue_pal()(length(sel_bcs)), sel_bcs)
  
  ggplot(bc_ts, aes(x = Days_pi, y = bc_capped, group = barcode, color = barcode)) +
    geom_line(linewidth = 0.35, alpha = 0.85, show.legend = FALSE) +
    scale_color_manual(values = pal) +
    scale_y_continuous(limits = c(0, cap),
                       name = "Barcode level") +          # <-- no "capped @ 25"
    labs(title = title_suffix, x = "Days post-injection") +
    theme_minimal()
}

safe_arrange <- function(plot_list, nrow = NULL, ncol = NULL) {
  ord <- c("FL","FR","ML","MR")
  have <- intersect(ord, names(plot_list))
  do.call(gridExtra::grid.arrange, c(plot_list[have], nrow = nrow, ncol = ncol))
}

make_panel_grob <- function(selection_tbl, panel_title, cap = CAP_VAL) {
  animals <- if (exists("animal_colors")) names(animal_colors) else c("FL","FR","ML","MR")
  
  grobs <- setNames(vector("list", length(animals)), animals)
  for (a in animals) {
    grobs[[a]] <- plot_barcode_lines_capped_for_mouse(
      urine_tbl    = urine_bc_levels,
      selected_tbl = selection_tbl,
      mouse_id     = a,
      title_suffix = a,
      cap          = cap
    )
  }
  
  #centered panel title above the 2x2 mouse grid
  gridExtra::arrangeGrob(
    safe_arrange(grobs, nrow = 2, ncol = 2),
    top = grid::textGrob(panel_title, gp = grid::gpar(fontsize = 10), just = "center")
  )
}


if (!exists("kidney_rank") || !exists("kidney_rank_notLastK")) {
  
  last_k_days <- urine_bc_levels %>%
    distinct(animal, Days_pi) %>%
    group_by(animal) %>%
    arrange(Days_pi, .by_group = TRUE) %>%
    slice_tail(n = LAST_K) %>%
    mutate(in_lastK = TRUE) %>%
    ungroup()
  
  present_lastK <- urine_bc_levels %>%
    inner_join(last_k_days, by = c("animal","Days_pi")) %>%
    group_by(animal, barcode) %>%
    summarise(
      present_lastK = if (REQUIRE_ALL_LASTK) all(bc_level > 0, na.rm = TRUE) else any(bc_level > 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(present_lastK) %>%
    distinct(animal, barcode)
  
  absent_lastK <- urine_bc_levels %>%
    distinct(animal, barcode) %>%
    anti_join(present_lastK, by = c("animal","barcode"))
  
  # Late-shed: top-N kidney barcodes among those present in last-K urine
  kidney_rank <- tissue_bc_levels %>%
    filter(is_kidney(organ)) %>%
    group_by(animal, barcode) %>%
    summarise(kidney_load = sum(bc_level, na.rm = TRUE), .groups = "drop") %>%
    inner_join(present_lastK, by = c("animal","barcode")) %>%
    group_by(animal) %>%
    arrange(desc(kidney_load), .by_group = TRUE) %>%
    slice_head(n = TOP_N) %>%
    ungroup() %>%
    select(animal, barcode)
  
  # Late non-shed: top-N kidney barcodes among those absent in last-K urine
  kidney_rank_notLastK <- tissue_bc_levels %>%
    filter(is_kidney(organ)) %>%
    group_by(animal, barcode) %>%
    summarise(kidney_load = sum(bc_level, na.rm = TRUE), .groups = "drop") %>%
    inner_join(absent_lastK, by = c("animal","barcode")) %>%
    group_by(animal) %>%
    arrange(desc(kidney_load), .by_group = TRUE) %>%
    slice_head(n = TOP_N) %>%
    ungroup() %>%
    select(animal, barcode)
}

selected_late_shed_tbl    <- kidney_rank %>% distinct(animal, barcode)
selected_late_nonshed_tbl <- kidney_rank_notLastK %>% distinct(animal, barcode)

#random control
selected_union <- bind_rows(selected_late_shed_tbl, selected_late_nonshed_tbl) %>%
  distinct(animal, barcode)

rand_tbl <- random_topN_control(
  urine_tbl    = urine_bc_levels,
  selected_tbl = selected_union,
  n = TOP_N,
  seed = 777
)


panelA_grob <- make_panel_grob(
  selection_tbl = selected_late_shed_tbl,
  panel_title   = paste0("Late-shed (", TOP_LABEL, " by kidney barcode level)"),
  cap           = CAP_VAL
)

panelB_grob <- make_panel_grob(
  selection_tbl = selected_late_nonshed_tbl,
  panel_title   = paste0("Late non-shed (", TOP_LABEL, " by kidney barcode level)"),
  cap           = CAP_VAL
)

panelC_grob <- make_panel_grob(
  selection_tbl = rand_tbl,
  panel_title   = "Random control",
  cap           = CAP_VAL
)


fig3_pdf <- file.path(SAVE_DIR, paste0("fig3_trajectory_panels_cap", CAP_VAL, "_", TOP_LABEL, ".pdf"))

pdf(fig3_pdf, width = 7.24, height = 9.5, onefile = TRUE)
gridExtra::grid.arrange(
  panelA_grob, panelB_grob, panelC_grob,
  ncol = 1,
  heights = c(1, 1, 1),
  top = NULL
)

grid::grid.text("A", x = unit(0.02, "npc"), y = unit(0.99, "npc"),
                gp = grid::gpar(fontsize = 14, fontface = "bold"))
grid::grid.text("B", x = unit(0.02, "npc"), y = unit(0.665, "npc"),
                gp = grid::gpar(fontsize = 14, fontface = "bold"))
grid::grid.text("C", x = unit(0.02, "npc"), y = unit(0.330, "npc"),
                gp = grid::gpar(fontsize = 14, fontface = "bold"))
dev.off()






##cap change

##### Fig 3: clone-resolved urine trajectories (cap=10), stacked A/B/C (one column)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(purrr)
  library(grid)
  library(gridExtra)
  library(scales)
  library(tibble)
})


TOP_N    <- 50
CAP_VAL  <- 5
LAST_K   <- 2
REQUIRE_ALL_LASTK <- FALSE

SAVE_DIR <- "../plots/urine_tissue_manuscript5"
dir.create(SAVE_DIR, showWarnings = FALSE, recursive = TRUE)

TOP_LABEL <- paste0("Top-", TOP_N)


is_kidney <- function(x) grepl("kidney|renal", x, ignore.case = TRUE)

random_topN_control <- function(urine_tbl, selected_tbl, n = TOP_N, seed = 1) {
  set.seed(seed)
  avail <- urine_tbl %>% distinct(animal, barcode)
  pool  <- anti_join(avail,
                     selected_tbl %>% distinct(animal, barcode),
                     by = c("animal","barcode"))
  pool %>%
    group_by(animal) %>%
    mutate(.rand = runif(n())) %>%
    arrange(.rand, .by_group = TRUE) %>%
    mutate(.rank = row_number()) %>%
    filter(.rank <= n) %>%
    ungroup() %>%
    select(animal, barcode)
}

plot_barcode_lines_capped_for_mouse <- function(urine_tbl, selected_tbl, mouse_id,
                                                title_suffix = "", cap = 25) {
  sel_bcs <- selected_tbl %>%
    filter(animal == mouse_id) %>%
    pull(barcode) %>% unique()
  
  if (length(sel_bcs) == 0) {
    return(ggplot() + theme_void() + labs(title = paste(mouse_id, "(no data)")))
  }
  
  bc_ts <- urine_tbl %>%
    semi_join(tibble(animal = mouse_id, barcode = sel_bcs),
              by = c("animal","barcode")) %>%
    group_by(animal, barcode, Days_pi) %>%
    summarise(bc_day_level = sum(bc_level, na.rm = TRUE), .groups = "drop") %>%
    mutate(bc_capped = pmin(bc_day_level, cap))
  
  pal <- setNames(scales::hue_pal()(length(sel_bcs)), sel_bcs)
  
  ggplot(bc_ts, aes(x = Days_pi, y = bc_capped, group = barcode, color = barcode)) +
    geom_line(linewidth = 0.35, alpha = 0.85, show.legend = FALSE) +
    scale_color_manual(values = pal) +
    scale_y_continuous(limits = c(0, cap),
                       name = "Barcode level") +          # <-- no "capped @ 25"
    labs(title = title_suffix, x = "Days post-injection") +
    theme_minimal()
}

safe_arrange <- function(plot_list, nrow = NULL, ncol = NULL) {
  ord <- c("FL","FR","ML","MR")
  have <- intersect(ord, names(plot_list))
  do.call(gridExtra::grid.arrange, c(plot_list[have], nrow = nrow, ncol = ncol))
}

make_panel_grob <- function(selection_tbl, panel_title, cap = CAP_VAL) {
  animals <- if (exists("animal_colors")) names(animal_colors) else c("FL","FR","ML","MR")
  
  grobs <- setNames(vector("list", length(animals)), animals)
  for (a in animals) {
    grobs[[a]] <- plot_barcode_lines_capped_for_mouse(
      urine_tbl    = urine_bc_levels,
      selected_tbl = selection_tbl,
      mouse_id     = a,
      title_suffix = a,
      cap          = cap
    )
  }
  
  
  gridExtra::arrangeGrob(
    safe_arrange(grobs, nrow = 2, ncol = 2),
    top = grid::textGrob(panel_title, gp = grid::gpar(fontsize = 10), just = "center")
  )
}


if (!exists("kidney_rank") || !exists("kidney_rank_notLastK")) {
  
  last_k_days <- urine_bc_levels %>%
    distinct(animal, Days_pi) %>%
    group_by(animal) %>%
    arrange(Days_pi, .by_group = TRUE) %>%
    slice_tail(n = LAST_K) %>%
    mutate(in_lastK = TRUE) %>%
    ungroup()
  
  present_lastK <- urine_bc_levels %>%
    inner_join(last_k_days, by = c("animal","Days_pi")) %>%
    group_by(animal, barcode) %>%
    summarise(
      present_lastK = if (REQUIRE_ALL_LASTK) all(bc_level > 0, na.rm = TRUE) else any(bc_level > 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(present_lastK) %>%
    distinct(animal, barcode)
  
  absent_lastK <- urine_bc_levels %>%
    distinct(animal, barcode) %>%
    anti_join(present_lastK, by = c("animal","barcode"))
  
  # Late-shed: top-N kidney barcodes among those present in last-K urine
  kidney_rank <- tissue_bc_levels %>%
    filter(is_kidney(organ)) %>%
    group_by(animal, barcode) %>%
    summarise(kidney_load = sum(bc_level, na.rm = TRUE), .groups = "drop") %>%
    inner_join(present_lastK, by = c("animal","barcode")) %>%
    group_by(animal) %>%
    arrange(desc(kidney_load), .by_group = TRUE) %>%
    slice_head(n = TOP_N) %>%
    ungroup() %>%
    select(animal, barcode)
  
  # Late non-shed: top-N kidney barcodes among those absent in last-K urine
  kidney_rank_notLastK <- tissue_bc_levels %>%
    filter(is_kidney(organ)) %>%
    group_by(animal, barcode) %>%
    summarise(kidney_load = sum(bc_level, na.rm = TRUE), .groups = "drop") %>%
    inner_join(absent_lastK, by = c("animal","barcode")) %>%
    group_by(animal) %>%
    arrange(desc(kidney_load), .by_group = TRUE) %>%
    slice_head(n = TOP_N) %>%
    ungroup() %>%
    select(animal, barcode)
}

selected_late_shed_tbl    <- kidney_rank %>% distinct(animal, barcode)
selected_late_nonshed_tbl <- kidney_rank_notLastK %>% distinct(animal, barcode)

#random control: exclude both selected sets
selected_union <- bind_rows(selected_late_shed_tbl, selected_late_nonshed_tbl) %>%
  distinct(animal, barcode)

rand_tbl <- random_topN_control(
  urine_tbl    = urine_bc_levels,
  selected_tbl = selected_union,
  n = TOP_N,
  seed = 777
)

panelA_grob <- make_panel_grob(
  selection_tbl = selected_late_shed_tbl,
  panel_title   = paste0("Late-shed (", TOP_LABEL, " by kidney barcode level)"),
  cap           = CAP_VAL
)

panelB_grob <- make_panel_grob(
  selection_tbl = selected_late_nonshed_tbl,
  panel_title   = paste0("Late non-shed (", TOP_LABEL, " by kidney barcode level)"),
  cap           = CAP_VAL
)

panelC_grob <- make_panel_grob(
  selection_tbl = rand_tbl,
  panel_title   = "Random control",
  cap           = CAP_VAL
)


fig3_pdf <- file.path(SAVE_DIR, paste0("fig3s_trajectory_panels_cap", CAP_VAL, "_", TOP_LABEL, ".pdf"))

pdf(fig3_pdf, width = 7.24, height = 9.5, onefile = TRUE)
gridExtra::grid.arrange(
  panelA_grob, panelB_grob, panelC_grob,
  ncol = 1,
  heights = c(1, 1, 1),
  top = NULL
)
grid::grid.text("A", x = unit(0.02, "npc"), y = unit(0.99, "npc"),
                gp = grid::gpar(fontsize = 14, fontface = "bold"))
grid::grid.text("B", x = unit(0.02, "npc"), y = unit(0.665, "npc"),
                gp = grid::gpar(fontsize = 14, fontface = "bold"))
grid::grid.text("C", x = unit(0.02, "npc"), y = unit(0.330, "npc"),
                gp = grid::gpar(fontsize = 14, fontface = "bold"))
dev.off()










suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(scales)
  library(grid)
  library(gridExtra)
  library(patchwork)
})

SAVE_DIR <- "../plots/urine_tissue_manuscript5"
dir.create(SAVE_DIR, showWarnings = FALSE, recursive = TRUE)

TOP_N  <- 50
LAST_K <- 2
REQUIRE_ALL_LASTK <- FALSE  # TRUE = present in *all* of the last K; FALSE = present in *any* of the last K

is_kidney <- function(x) grepl("kidney|renal", x, ignore.case = TRUE)


if (!exists("animal_colors")) {
  animals <- c("FL","FR","ML","MR")
} else {
  animals <- names(animal_colors)
}


plot_diversity_donut_w_colored_barcodes <- function(bc_table,
                                                    barcode_color_highlights,
                                                    faceting_variable) {
  
  valid_highlights <- barcode_color_highlights[barcode_color_highlights %in% bc_table$barcode]
  if (length(valid_highlights) == 0) stop("None of the barcode highlights are present in the dataset.")
  
  winner_and_bg_colors <- set_names(
    c(scales::hue_pal()(length(valid_highlights)), "gray75", "gray60"),
    nm = c(valid_highlights, "background1", "background2")
  )
  
  bc_table <- bc_table %>%
    group_by(sample) %>%
    arrange(desc(frac_weighted_count), .by_group = TRUE) %>%
    mutate(
      ymax = cumsum(frac_weighted_count),
      ymin = c(0, head(ymax, n = -1)),
      bc_color = ifelse(barcode %in% valid_highlights,
                        barcode,
                        ifelse(row_number() %% 2 == 0, "background1", "background2"))
    ) %>%
    ungroup()
  
  ggplot(bc_table, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3.2, fill = bc_color)) +
    scale_fill_manual(values = winner_and_bg_colors, labels = names(winner_and_bg_colors)) +
    geom_rect(show.legend = FALSE) +
    coord_polar(theta = "y") +
    xlim(c(2, 4)) +
    facet_wrap(vars({{faceting_variable}})) +
    theme_classic() +
    theme(
      line = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    ) +
    guides(fill = guide_legend("Barcode"))
}


wrap_with_title <- function(grob_obj, title_text) {
  gridExtra::arrangeGrob(
    grob_obj,
    top = grid::textGrob(title_text, gp = grid::gpar(fontsize = 10), just = "center")
  )
}


last_k_days <- urine_bc_levels %>%
  distinct(animal, Days_pi) %>%
  group_by(animal) %>%
  arrange(Days_pi, .by_group = TRUE) %>%
  slice_tail(n = LAST_K) %>%
  mutate(in_lastK = TRUE) %>%
  ungroup()

present_lastK <- urine_bc_levels %>%
  inner_join(last_k_days, by = c("animal","Days_pi")) %>%
  group_by(animal, barcode) %>%
  summarise(
    present_lastK = if (REQUIRE_ALL_LASTK) all(bc_level > 0, na.rm = TRUE) else any(bc_level > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(present_lastK) %>%
  distinct(animal, barcode)

absent_lastK <- urine_bc_levels %>%
  distinct(animal, barcode) %>%
  anti_join(present_lastK, by = c("animal","barcode"))


kidney_rank <- tissue_bc_levels %>%
  filter(is_kidney(organ)) %>%
  group_by(animal, barcode) %>%
  summarise(kidney_load = sum(bc_level, na.rm = TRUE), .groups = "drop") %>%
  inner_join(present_lastK, by = c("animal","barcode")) %>%
  group_by(animal) %>%
  arrange(desc(kidney_load), .by_group = TRUE) %>%
  slice_head(n = TOP_N) %>%
  ungroup()

winner_barcodes_by_animal <- kidney_rank %>%
  group_by(animal) %>%
  summarise(barcodes = list(barcode), .groups = "drop") %>%
  deframe()

donut_plot_list <- lapply(animals, function(a) {
  bc_table <- urine_bc_levels %>% filter(animal == a)
  highlights <- winner_barcodes_by_animal[[a]]
  if (is.null(highlights) || length(highlights) == 0) return(ggplot() + theme_void() + labs(title = a))
  
  plot_diversity_donut_w_colored_barcodes(
    bc_table = bc_table,
    barcode_color_highlights = highlights,
    faceting_variable = Days_pi
  ) + labs(title = a)
})

donut_combined <- gridExtra::arrangeGrob(grobs = donut_plot_list, ncol = 2)
A_grob <- wrap_with_title(donut_combined, "Urine barcode composition with top-50 kidney barcodes highlighted (Late-shed)")


kidney_rank_notLastK <- tissue_bc_levels %>%
  filter(is_kidney(organ)) %>%
  group_by(animal, barcode) %>%
  summarise(kidney_load = sum(bc_level, na.rm = TRUE), .groups = "drop") %>%
  inner_join(absent_lastK, by = c("animal","barcode")) %>%
  group_by(animal) %>%
  arrange(desc(kidney_load), .by_group = TRUE) %>%
  slice_head(n = TOP_N) %>%
  ungroup()

winner_barcodes_by_animal_notLastK <- kidney_rank_notLastK %>%
  group_by(animal) %>%
  summarise(barcodes = list(barcode), .groups = "drop") %>%
  deframe()

donut_plot_list_notLastK <- lapply(animals, function(a) {
  bc_table <- urine_bc_levels %>% filter(animal == a)
  highlights <- winner_barcodes_by_animal_notLastK[[a]]
  
  if (is.null(highlights) || length(highlights) == 0) return(ggplot() + theme_void() + labs(title = a))
  
  highlights_in_urine <- intersect(highlights, bc_table$barcode)
  if (length(highlights_in_urine) == 0) return(ggplot() + theme_void() + labs(title = a))
  
  plot_diversity_donut_w_colored_barcodes(
    bc_table = bc_table,
    barcode_color_highlights = highlights_in_urine,
    faceting_variable = Days_pi
  ) + labs(title = a)
})

donut_combined_notLastK <- gridExtra::arrangeGrob(grobs = donut_plot_list_notLastK, ncol = 2)
B_grob <- wrap_with_title(donut_combined_notLastK, "Urine barcode composition with top-50 kidney barcodes highlighted (Late non-shed)")


urine_bc_kidTop <- urine_bc_levels %>%
  semi_join(kidney_rank, by = c("animal","barcode"))

sum_selected <- urine_bc_kidTop %>%
  group_by(animal, Days_pi, mean_ul_urine) %>%
  summarise(sum_selected_barcodes = sum(bc_level, na.rm = TRUE), .groups = "drop")

sum_overlay_plots <- sum_selected %>%
  split(.$animal) %>%
  imap(~ ggplot(.x, aes(x = Days_pi)) +
         geom_area(aes(y = mean_ul_urine), fill = "grey85") +
         geom_area(aes(y = sum_selected_barcodes), fill = "navajowhite3") +
         labs(title = .y, x = "", y = "") +
         theme_minimal() +
         theme(legend.position = "none"))

overlay_combined <- gridExtra::arrangeGrob(
  sum_overlay_plots[["FL"]],
  sum_overlay_plots[["FR"]],
  sum_overlay_plots[["ML"]],
  sum_overlay_plots[["MR"]],
  nrow = 4, ncol = 1,
  bottom = grid::textGrob("Days Post-Injection", gp = grid::gpar(fontsize = 11)),
  left   = grid::textGrob("Barcode level", rot = 90, gp = grid::gpar(fontsize = 11))
)
C_grob <- wrap_with_title(overlay_combined, "Summed urine levels of top-50 kidney barcodes over time (Late-shed)")


urine_bc_kidTop_notLastK <- urine_bc_levels %>%
  semi_join(kidney_rank_notLastK, by = c("animal","barcode"))

sum_selected_notLastK <- urine_bc_kidTop_notLastK %>%
  group_by(animal, Days_pi, mean_ul_urine) %>%
  summarise(sum_selected_barcodes = sum(bc_level, na.rm = TRUE), .groups = "drop")

sum_overlay_plots_notLastK <- sum_selected_notLastK %>%
  split(.$animal) %>%
  imap(~ ggplot(.x, aes(x = Days_pi)) +
         geom_area(aes(y = mean_ul_urine), fill = "grey85") +
         geom_area(aes(y = sum_selected_barcodes), fill = "navajowhite3") +
         labs(title = .y, x = "", y = "") +
         theme_minimal() +
         theme(legend.position = "none"))

overlay_combined_notLastK <- gridExtra::arrangeGrob(
  sum_overlay_plots_notLastK[["FL"]],
  sum_overlay_plots_notLastK[["FR"]],
  sum_overlay_plots_notLastK[["ML"]],
  sum_overlay_plots_notLastK[["MR"]],
  nrow = 4, ncol = 1,
  bottom = grid::textGrob("Days Post-Injection", gp = grid::gpar(fontsize = 11)),
  left   = grid::textGrob("Barcode level", rot = 90, gp = grid::gpar(fontsize = 11))
)
D_grob <- wrap_with_title(overlay_combined_notLastK, "Summed urine levels of top-50 kidney barcodes over time (Late non-shed)")


pA <- patchwork::wrap_elements(full = A_grob)
pB <- patchwork::wrap_elements(full = B_grob)
pC <- patchwork::wrap_elements(full = C_grob)
pD <- patchwork::wrap_elements(full = D_grob)

fig2AB <- (pA | pB) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 12),
    plot.tag.position = c(0.0, .99995)
  )

print(fig2AB)

ggsave(
  filename = "Fig2AB_donuts.pdf",
  plot     = fig2AB,
  path     = SAVE_DIR,
  width    = 12, height = 10
)


fig2CD <- (pC | pD) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = list(c("C","D"))) &
  theme(
    plot.tag = element_text(face = "bold", size = 12),
    plot.tag.position = c(0.0, .99995)
  )

print(fig2CD)

ggsave(
  filename = "Fig2CD_sum_overlays.pdf",
  plot     = fig2CD,
  path     = SAVE_DIR,
  width    = 12, height = 7
)













##### Fig 4: Stock-rank scatter (Late-shed on top, Late non-shed on bottom)


suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(forcats)
  library(patchwork)
  library(grid)
})

SAVE_DIR <- "../plots/urine_tissue_manuscript5"
dir.create(SAVE_DIR, showWarnings = FALSE, recursive = TRUE)


stock_rank_tbl <- cutoff_99pct_stock_barcodes %>%
  mutate(stock_rank = min_rank(desc(count))) %>%
  select(barcode, count, stock_rank)

max_stock_rank <- nrow(cutoff_99pct_stock_barcodes)


rank_plot_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "plain", size = 10, hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(linetype = "dotted", color = "gray80", linewidth = 0.3),
    panel.grid.minor.y = element_blank(),
    panel.spacing = unit(2, "lines"),
    strip.background = element_rect(fill = "gray90", color = "gray90"),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8)
  )

pt_size  <- .5
pt_alpha <- 1


kidney_present_with_stock <- kidney_rank %>%
  left_join(stock_rank_tbl, by = "barcode") %>%
  mutate(kidney_rank_f = as_factor(kidney_rank))

p_present <- kidney_present_with_stock %>%
  ggplot(aes(x = fct_relevel(kidney_rank_f, rev), y = stock_rank)) +
  geom_point(aes(color = as_factor(animal)),
             size = pt_size, alpha = pt_alpha, show.legend = FALSE) +
  facet_wrap(vars(animal)) +
  labs(
    title = "Late-shed kidney barcodes: input stock rank (Top 50)",
    x = "Kidney barcode rank (within group)",
    y = "Rank in stock"
  ) +
  scale_y_continuous(limits = c(1, max_stock_rank)) +
  scale_x_discrete(breaks = as.character(seq(0, TOP_N, by = 5))) +
  scale_color_manual(values = animal_colors) +
  coord_flip() +
  rank_plot_theme

kidney_present_median_stock_rank <- kidney_present_with_stock %>%
  group_by(animal) %>%
  summarise(median_stock_rank = median(stock_rank, na.rm = TRUE), .groups = "drop")

print(kidney_present_median_stock_rank)


kidney_notLastK_with_stock <- kidney_rank_notLastK %>%
  left_join(stock_rank_tbl, by = "barcode") %>%
  mutate(kidney_rank_f = as_factor(kidney_rank))

p_notLastK <- kidney_notLastK_with_stock %>%
  ggplot(aes(x = fct_relevel(kidney_rank_f, rev), y = stock_rank)) +
  geom_point(aes(color = as_factor(animal)),
             size = pt_size, alpha = pt_alpha, show.legend = FALSE) +
  facet_wrap(vars(animal)) +
  labs(
    title = "Late non-shed kidney barcodes: input stock rank (Top 50)",
    x = "Kidney barcode rank (within group)",
    y = "Rank in stock"
  ) +
  scale_y_continuous(limits = c(1, max_stock_rank)) +
  scale_x_discrete(breaks = as.character(seq(0, TOP_N, by = 5))) +
  scale_color_manual(values = animal_colors) +
  coord_flip() +
  rank_plot_theme

kidney_notLastK_median_stock_rank <- kidney_notLastK_with_stock %>%
  group_by(animal) %>%
  summarise(median_stock_rank = median(stock_rank, na.rm = TRUE), .groups = "drop")

print(kidney_notLastK_median_stock_rank)


fig4_panel <- p_present / p_notLastK +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 12),
    plot.tag.position = c(0.0, 1)
  )

print(fig4_panel)

ggsave(
  filename = "fig4_stock_rank_panel.pdf",
  plot     = fig4_panel,
  path     = SAVE_DIR,
  width    = 7, height = 9
)












# Fig 5: Peak urine barcode level per barcode (Top-N cohorts), faceted by mouse

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(purrr)
  library(forcats)
})


SMOLDER_TOP_N  <- if (exists("TOP_N")) TOP_N else 50
LOG_EPS        <- 1e-6
SAVE_DIR_LOCAL <- if (exists("SAVE_DIR")) SAVE_DIR else "."
dir.create(SAVE_DIR_LOCAL, showWarnings = FALSE, recursive = TRUE)

LABEL_RANDOM   <- "Random control"
LABEL_SHED     <- "Late-shed"
LABEL_NONSHED  <- "Late non-shed"
cohort_levels  <- c(LABEL_RANDOM, LABEL_SHED, LABEL_NONSHED)


scale_def_colors <- scale_fill_manual(values = c(
  "Late-shed"     = "#F8766D",
  "Late non-shed" = "#00BFC4",
  "Random control"= "#7F7F7F"   # neutral gray, color-blind friendly
))

scale_def_colors_point <- scale_color_manual(values = c(
  "Late-shed"     = "#F8766D",
  "Late non-shed" = "#00BFC4",
  "Random control"= "#7F7F7F"
))


random_topN_from_urine <- function(urine_tbl, exclude_tbl, n = SMOLDER_TOP_N, seed = 123) {
  set.seed(seed)
  avail <- urine_tbl %>% distinct(animal, barcode)
  pool  <- anti_join(avail, exclude_tbl %>% distinct(animal, barcode),
                     by = c("animal","barcode"))
  pool %>%
    group_by(animal) %>%
    mutate(.rand = runif(n())) %>%
    arrange(.rand, .by_group = TRUE) %>%
    mutate(.rank = row_number()) %>%
    filter(.rank <= n) %>%
    ungroup() %>%
    select(animal, barcode)
}


kidney_notLastK_topN <- if (exists("kidney_rank_notLastK_topN")) kidney_rank_notLastK_topN else kidney_rank_notLastK
cohort_kidney_notLastK <- kidney_notLastK_topN %>%
  distinct(animal, barcode) %>%
  mutate(cohort = LABEL_NONSHED)

kidney_present_topN <- if (exists("kidney_rank_topN")) kidney_rank_topN else kidney_rank
cohort_kidney_present <- kidney_present_topN %>%
  distinct(animal, barcode) %>%
  mutate(cohort = LABEL_SHED)

# Random control: exclude BOTH late-shed and late non-shed sets (so it’s a true baseline)
exclude_union <- bind_rows(cohort_kidney_present, cohort_kidney_notLastK) %>%
  distinct(animal, barcode)

random_for_baseline <- random_topN_from_urine(
  urine_tbl   = urine_bc_levels,
  exclude_tbl = exclude_union,
  n           = SMOLDER_TOP_N,
  seed        = 777
) %>%
  mutate(cohort = LABEL_RANDOM)

cohorts_all <- bind_rows(
  random_for_baseline,
  cohort_kidney_present,
  cohort_kidney_notLastK
) %>%
  distinct(animal, barcode, cohort) %>%
  mutate(cohort = factor(cohort, levels = cohort_levels))


max_urine_bc <- urine_bc_levels %>%
  group_by(animal, barcode) %>%
  summarise(max_bc = max(bc_level, na.rm = TRUE), .groups = "drop")

cohorts_with_max <- cohorts_all %>%
  left_join(max_urine_bc, by = c("animal","barcode")) %>%
  mutate(
    max_bc = replace_na(max_bc, 0),
    max_bc_plot = pmax(max_bc, LOG_EPS)
  )


p_max_bc_violin <- ggplot(
  cohorts_with_max,
  aes(x = cohort, y = max_bc_plot, fill = cohort)
) +
  geom_violin(trim = TRUE, scale = "width", color = "grey30", alpha = 0.85) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.95) +
  geom_jitter(color = "black", width = 0.1, height = 0, alpha = 0.5, size = 0.9, show.legend = FALSE) +
  scale_y_log10() +
  scale_def_colors +
  scale_def_colors_point +
  facet_wrap(vars(animal), nrow = 2) +
  labs(
    title = paste0("Peak urine barcode level per barcode (Top ", SMOLDER_TOP_N, " per cohort)"),
    x     = "",
    y     = "Peak urine barcode level across timepoints (log10)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid.minor.y = element_blank()
  )

print(p_max_bc_violin)

ggsave(
  filename = paste0("fig5_peak_urine_level_by_cohort_TopN", SMOLDER_TOP_N, ".pdf"),
  plot     = p_max_bc_violin,
  path     = SAVE_DIR_LOCAL,
  width    = 7.24, height = 5
)








######### Fig 6: Top-50 smolderers enrichment in kidney categories vs random

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(scales)
})


SMOLDER_TOP_N  <- 50
DETECT         <- 1e-100
SAVE_DIR_LOCAL <- if (exists("SAVE_DIR")) SAVE_DIR else "."
dir.create(SAVE_DIR_LOCAL, showWarnings = FALSE, recursive = TRUE)


LABEL_SHED     <- "Late-shed"
LABEL_NONSHED  <- "Late non-shed"
LABEL_RANDOM   <- "Random control"
LABEL_TOP      <- "Top smolderers"

category_levels <- c(LABEL_SHED, LABEL_NONSHED)
group_levels    <- c(LABEL_RANDOM, LABEL_TOP)

FILL_TOP_SHED    <- paste0(LABEL_TOP, " in ", LABEL_SHED)
FILL_TOP_NONSHED <- paste0(LABEL_TOP, " in ", LABEL_NONSHED)

fill_levels <- c(LABEL_RANDOM, FILL_TOP_SHED, FILL_TOP_NONSHED)

scale_def_colors <- scale_fill_manual(
  breaks = fill_levels,
  values = c(
    "Random control"                 = "#7F7F7F",
    "Top smolderers in Late-shed"    = "#F8766D",
    "Top smolderers in Late non-shed"= "#00BFC4"
  )
)


compute_pct_detected_per_barcode <- function(urine_tbl, selected_tbl, detect_thresh = DETECT) {
  denom_tbl <- urine_tbl %>%
    distinct(animal, Days_pi) %>%
    count(animal, name = "n_timepoints")
  
  numer_tbl <- urine_tbl %>%
    semi_join(selected_tbl %>% distinct(animal, barcode), by = c("animal","barcode")) %>%
    group_by(animal, barcode, Days_pi) %>%
    summarise(bc_day_level = sum(bc_level, na.rm = TRUE), .groups = "drop") %>%
    mutate(detected = bc_day_level >= detect_thresh) %>%
    group_by(animal, barcode) %>%
    summarise(n_detected_days = sum(detected, na.rm = TRUE), .groups = "drop")
  
  numer_tbl %>%
    left_join(denom_tbl, by = "animal") %>%
    mutate(pct_detected = 100 * n_detected_days / n_timepoints) %>%
    left_join(selected_tbl %>% distinct(animal, barcode), by = c("animal","barcode"))
}


random_topN_from_urine <- function(urine_tbl, exclude_tbl, n = SMOLDER_TOP_N, seed = 999) {
  set.seed(seed)
  avail <- urine_tbl %>% distinct(animal, barcode)
  pool  <- anti_join(avail, exclude_tbl %>% distinct(animal, barcode),
                     by = c("animal","barcode"))
  pool %>%
    group_by(animal) %>%
    mutate(.rand = runif(n())) %>%
    arrange(.rand, .by_group = TRUE) %>%
    mutate(.rank = row_number()) %>%
    filter(.rank <= n) %>%
    ungroup() %>%
    select(animal, barcode)
}


kidney_present_topN  <- if (exists("kidney_rank_topN"))           kidney_rank_topN          else kidney_rank
kidney_notLastK_topN <- if (exists("kidney_rank_notLastK_topN"))  kidney_rank_notLastK_topN else kidney_rank_notLastK

categories_all <- bind_rows(
  kidney_present_topN  %>% distinct(animal, barcode) %>% mutate(category = LABEL_SHED),
  kidney_notLastK_topN %>% distinct(animal, barcode) %>% mutate(category = LABEL_NONSHED)
) %>%
  distinct(animal, barcode, category)


urine_universe <- urine_bc_levels %>% distinct(animal, barcode)

smolder_all <- compute_pct_detected_per_barcode(
  urine_tbl     = urine_bc_levels,
  selected_tbl  = urine_universe,
  detect_thresh = DETECT
) %>%
  mutate(pct_detected = replace_na(pct_detected, 0))

top_smolderers <- smolder_all %>%
  group_by(animal) %>%
  arrange(desc(pct_detected), .by_group = TRUE) %>%
  slice_head(n = SMOLDER_TOP_N) %>%
  ungroup() %>%
  select(animal, barcode) %>%
  mutate(group = LABEL_TOP)

rand_smolder_ctrl <- random_topN_from_urine(
  urine_tbl   = urine_bc_levels,
  exclude_tbl = top_smolderers,
  n           = SMOLDER_TOP_N,
  seed        = 2025
) %>%
  mutate(group = LABEL_RANDOM)

smolder_vs_random <- bind_rows(top_smolderers, rand_smolder_ctrl) %>%
  distinct(animal, barcode, group)


smolder_cat <- smolder_vs_random %>%
  left_join(categories_all, by = c("animal","barcode")) %>%
  mutate(category = replace_na(category, "None"))


enrich_summary <- smolder_cat %>%
  group_by(animal, category, group) %>%
  summarise(n_in_cat = n(), .groups = "drop_last") %>%
  group_by(animal, group) %>%
  mutate(
    n_total_group = sum(n_in_cat),
    frac_in_cat   = n_in_cat / pmax(n_total_group, 1)
  ) %>%
  ungroup() %>%
  filter(category %in% category_levels) %>%
  mutate(
    category = factor(category, levels = category_levels),
    group    = factor(group, levels = group_levels)
  )


enrich_complete <- enrich_summary %>%
  tidyr::complete(
    animal,
    category,
    group,
    fill = list(frac_in_cat = 0)
  )


enrich_complete <- enrich_complete %>%
  mutate(
    fill_key = case_when(
      group == LABEL_RANDOM ~ LABEL_RANDOM,
      group == LABEL_TOP & category == LABEL_SHED    ~ FILL_TOP_SHED,
      group == LABEL_TOP & category == LABEL_NONSHED ~ FILL_TOP_NONSHED,
      TRUE ~ LABEL_RANDOM
    ),
    fill_key = factor(fill_key, levels = fill_levels)
  )


p_frac_cat <- ggplot(
  enrich_complete,
  aes(x = category, y = frac_in_cat, fill = fill_key)
) +
  geom_col(
    position = position_dodge(width = 0.7, preserve = "single"),
    width = 0.6
  ) +
  facet_wrap(vars(animal), nrow = 2) +
  labs(
    title = paste0("Top-", SMOLDER_TOP_N, " smolderers: kidney category membership vs random control"),
    x     = "",
    y     = "Fraction of barcodes",
    fill  = ""
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_def_colors +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )

print(p_frac_cat)

ggsave(
  filename = paste0("fig6_enrichment_TopSmolder_vs_Random_kidneyCategories_Top", SMOLDER_TOP_N, ".pdf"),
  plot     = p_frac_cat,
  path     = SAVE_DIR_LOCAL,
  width    = 7, height = 5
)






########## Fig 7: Early detection frequency for barcode groups defined by late-phase behavior


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
})

TOP_LATE_N         <- if (exists("TOP_N")) TOP_N else 50
DETECT             <- 1e-116
LATE_FRACTION      <- 0.25
EARLY_SMOLDER_FRAC <- 0.5
SAVE_DIR_LOCAL     <- if (exists("SAVE_DIR")) SAVE_DIR else "."
dir.create(SAVE_DIR_LOCAL, showWarnings = FALSE, recursive = TRUE)


urine_tp <- urine_bc_levels %>%
  distinct(animal, Days_pi) %>%
  group_by(animal) %>%
  arrange(Days_pi, .by_group = TRUE) %>%
  mutate(
    tp_index = row_number(),
    n_tp     = n(),
    frac_pos = tp_index / n_tp,
    is_late  = frac_pos > (1 - LATE_FRACTION),
    is_early = !is_late
  ) %>%
  ungroup()

urine_flagged <- urine_bc_levels %>%
  inner_join(urine_tp %>% select(animal, Days_pi, is_early, is_late),
             by = c("animal","Days_pi"))

denom_tbl <- urine_tp %>%
  group_by(animal) %>%
  summarise(
    n_early_tp = sum(is_early),
    n_late_tp  = sum(is_late),
    .groups = "drop"
  )

bc_day <- urine_flagged %>%
  group_by(animal, barcode, Days_pi, is_early, is_late) %>%
  summarise(bc_day_level = sum(bc_level, na.rm = TRUE), .groups = "drop")

smolder_early_late <- bc_day %>%
  mutate(detected = bc_day_level >= DETECT) %>%
  group_by(animal, barcode) %>%
  summarise(
    n_early_detect  = sum(detected & is_early, na.rm = TRUE),
    n_late_detect   = sum(detected & is_late,  na.rm = TRUE),
    late_total_load = sum(if_else(is_late, bc_day_level, 0), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(denom_tbl, by = "animal") %>%
  mutate(
    pct_detect_early = if_else(n_early_tp > 0, 100 * n_early_detect / n_early_tp, NA_real_),
    pct_detect_late  = if_else(n_late_tp  > 0, 100 * n_late_detect  / n_late_tp,  NA_real_)
  )


LABEL_RANDOM   <- "Random control"
LABEL_HIGH     <- "Top late high-shedders"
LABEL_SMOLDER  <- "Top late smolderers"

topLate_smolder <- smolder_early_late %>%
  group_by(animal) %>%
  arrange(desc(pct_detect_late), desc(late_total_load), .by_group = TRUE) %>%
  slice_head(n = TOP_LATE_N) %>%
  ungroup() %>%
  transmute(animal, barcode, group = LABEL_SMOLDER)

topLate_highshed <- smolder_early_late %>%
  group_by(animal) %>%
  arrange(desc(late_total_load), .by_group = TRUE) %>%
  slice_head(n = TOP_LATE_N) %>%
  ungroup() %>%
  transmute(animal, barcode, group = LABEL_HIGH)

urine_universe <- urine_bc_levels %>% distinct(animal, barcode)

set.seed(2026)
random_late <- urine_universe %>%
  anti_join(bind_rows(topLate_smolder, topLate_highshed),
            by = c("animal","barcode")) %>%
  group_by(animal) %>%
  mutate(.rand = runif(n())) %>%
  arrange(.rand, .by_group = TRUE) %>%
  slice_head(n = TOP_LATE_N) %>%
  ungroup() %>%
  transmute(animal, barcode, group = LABEL_RANDOM)

all_sets <- bind_rows(topLate_smolder, topLate_highshed, random_late) %>%
  left_join(smolder_early_late, by = c("animal","barcode")) %>%
  mutate(
    group = factor(group, levels = c(LABEL_RANDOM, LABEL_HIGH, LABEL_SMOLDER)),
    early_smolder = !is.na(pct_detect_early) & pct_detect_early >= (EARLY_SMOLDER_FRAC * 100)
  )


group_fill <- c(
  "Random control"        = "#7F7F7F",
  "Top late high-shedders"= "#E69F00",
  "Top late smolderers"   = "#0072B2"
)


p_early_pct <- all_sets %>%
  filter(!is.na(pct_detect_early)) %>%
  ggplot(aes(x = group, y = pct_detect_early, fill = group)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.7) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(color = "black", width = 0.08, height = 0, size = 0.9, alpha = 0.35) +
  facet_wrap(vars(animal)) +
  labs(
    title = paste0("Early detection frequency of late-defined barcode groups (Top ", TOP_LATE_N, ")"),
    x = "",
    y = "Early urine timepoints with detection (%)"
  ) +
  scale_fill_manual(values = group_fill) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.major.x = element_blank()
  )

print(p_early_pct)

ggsave(
  filename = paste0("fig7_early_detection_by_lateDefinedGroups_Top", TOP_LATE_N, ".pdf"),
  plot     = p_early_pct,
  path     = SAVE_DIR_LOCAL,
  width    = 7.24, height = 7
)



urine_tp

late_day_cutoff_tbl <- urine_tp %>%
  group_by(animal) %>%
  summarise(
    max_day = max(Days_pi),
    late_day_cutoff = max_day * (1 - LATE_FRACTION),
    .groups = "drop"
  )

print(late_day_cutoff_tbl)







# alt version (time-based late window): uses day cutoff = max_day*(1-0.25)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
})

TOP_LATE_N_alt         <- if (exists("TOP_N_alt")) TOP_N_alt else 50
DETECT_alt             <- 1e-116
LATE_FRACTION_alt      <- 0.25
EARLY_SMOLDER_FRAC_alt <- 0.5
SAVE_DIR_alt           <- SAVE_DIR
SAVE_DIR

late_day_cutoff_tbl_alt <- urine_bc_levels %>%
  distinct(animal, Days_pi) %>%
  group_by(animal) %>%
  summarise(max_day_alt = max(Days_pi, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    late_day_cutoff_alt = max_day_alt * (1 - LATE_FRACTION_alt),
    late_rule_alt = paste0("Days_pi > ", signif(late_day_cutoff_alt, 4),
                           " (=", max_day_alt, "*", (1 - LATE_FRACTION_alt), ")")
  )

print(late_day_cutoff_tbl_alt)


urine_tp_alt <- urine_bc_levels %>%
  distinct(animal, Days_pi) %>%
  left_join(late_day_cutoff_tbl_alt, by = "animal") %>%
  mutate(
    is_late_alt  = Days_pi > late_day_cutoff_alt,
    is_early_alt = !is_late_alt
  ) %>%
  arrange(animal, Days_pi)

urine_flagged_alt <- urine_bc_levels %>%
  inner_join(urine_tp_alt %>% select(animal, Days_pi, is_early_alt, is_late_alt),
             by = c("animal","Days_pi"))

denom_tbl_alt <- urine_tp_alt %>%
  group_by(animal) %>%
  summarise(
    n_early_tp_alt = sum(is_early_alt, na.rm = TRUE),
    n_late_tp_alt  = sum(is_late_alt,  na.rm = TRUE),
    .groups = "drop"
  )

bc_day_alt <- urine_flagged_alt %>%
  group_by(animal, barcode, Days_pi, is_early_alt, is_late_alt) %>%
  summarise(bc_day_level_alt = sum(bc_level, na.rm = TRUE), .groups = "drop")

smolder_early_late_alt <- bc_day_alt %>%
  mutate(detected_alt = bc_day_level_alt >= DETECT_alt) %>%
  group_by(animal, barcode) %>%
  summarise(
    n_early_detect_alt  = sum(detected_alt & is_early_alt, na.rm = TRUE),
    n_late_detect_alt   = sum(detected_alt & is_late_alt,  na.rm = TRUE),
    late_total_load_alt = sum(if_else(is_late_alt, bc_day_level_alt, 0), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(denom_tbl_alt, by = "animal") %>%
  mutate(
    pct_detect_early_alt = if_else(n_early_tp_alt > 0,
                                   100 * n_early_detect_alt / n_early_tp_alt, NA_real_),
    pct_detect_late_alt  = if_else(n_late_tp_alt  > 0,
                                   100 * n_late_detect_alt  / n_late_tp_alt,  NA_real_)
  ) %>%
  left_join(late_day_cutoff_tbl_alt %>% select(animal, max_day_alt, late_day_cutoff_alt, late_rule_alt),
            by = "animal")


LABEL_RANDOM_alt   <- "Random control"
LABEL_HIGH_alt     <- "Top late high-shedders"
LABEL_SMOLDER_alt  <- "Top late smolderers"

topLate_smolder_alt <- smolder_early_late_alt %>%
  group_by(animal) %>%
  arrange(desc(pct_detect_late_alt), desc(late_total_load_alt), .by_group = TRUE) %>%
  slice_head(n = TOP_LATE_N_alt) %>%
  ungroup() %>%
  transmute(animal, barcode, group_alt = LABEL_SMOLDER_alt)

topLate_highshed_alt <- smolder_early_late_alt %>%
  group_by(animal) %>%
  arrange(desc(late_total_load_alt), .by_group = TRUE) %>%
  slice_head(n = TOP_LATE_N_alt) %>%
  ungroup() %>%
  transmute(animal, barcode, group_alt = LABEL_HIGH_alt)

urine_pool_alt <- urine_bc_levels %>% distinct(animal, barcode)

set.seed(2100)
random_late_alt <- urine_pool_alt %>%
  anti_join(bind_rows(topLate_smolder_alt, topLate_highshed_alt),
            by = c("animal","barcode")) %>%
  group_by(animal) %>%
  mutate(.rand_alt = runif(n())) %>%
  arrange(.rand_alt, .by_group = TRUE) %>%
  slice_head(n = TOP_LATE_N_alt) %>%
  ungroup() %>%
  transmute(animal, barcode, group_alt = LABEL_RANDOM_alt)

all_sets_alt <- bind_rows(topLate_smolder_alt, topLate_highshed_alt, random_late_alt) %>%
  left_join(smolder_early_late_alt, by = c("animal","barcode")) %>%
  mutate(
    group_alt = factor(group_alt, levels = c(LABEL_RANDOM_alt, LABEL_HIGH_alt, LABEL_SMOLDER_alt)),
    early_smolder_alt = !is.na(pct_detect_early_alt) &
      pct_detect_early_alt >= (EARLY_SMOLDER_FRAC_alt * 100)
  )


group_fill_alt <- c(
  "Random control"         = "#7F7F7F",
  "Top late high-shedders" = "#E69F00",
  "Top late smolderers"    = "#0072B2"
)

p_early_pct_alt <- all_sets_alt %>%
  filter(!is.na(pct_detect_early_alt)) %>%
  ggplot(aes(x = group_alt, y = pct_detect_early_alt, fill = group_alt)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.7) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(color = "black", width = 0.08, height = 0, size = 0.9, alpha = 0.35) +
  facet_wrap(vars(animal)) +
  labs(
    title = paste0(
      "Early detection frequency of late-defined barcode groups (Top ", TOP_LATE_N_alt, ")"
    ),
    x = "",
    y = "Early urine timepoints with detection (%)"
  ) +
  scale_fill_manual(values = group_fill_alt) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.major.x = element_blank()
  )

print(p_early_pct_alt)

ggsave(
  filename = paste0("fig7_update_alt_dayCutoff_early_detection_by_lateDefinedGroups_Top", TOP_LATE_N_alt, ".pdf"),
  plot     = p_early_pct_alt,
  path     = SAVE_DIR_alt,
  width    = 7.24, height = 7
)


supp_tbl_alt <- all_sets_alt %>%
  left_join(late_day_cutoff_tbl_alt, by = "animal") %>%
  select(
    animal,
    max_day_alt,
    late_day_cutoff_alt,
    late_rule_alt,
    barcode,
    group_alt,
    n_early_tp_alt,
    n_late_tp_alt,
    n_early_detect_alt,
    n_late_detect_alt,
    pct_detect_early_alt,
    pct_detect_late_alt,
    late_total_load_alt,
    early_smolder_alt
  ) %>%
  arrange(animal, group_alt, desc(pct_detect_late_alt), desc(late_total_load_alt))

write_tsv(
  supp_tbl_alt,
  file = file.path(SAVE_DIR_alt,
                   paste0("supp_table_fig7_alt_dayCutoff_Top", TOP_LATE_N_alt, ".tsv"))
)

write_tsv(
  late_day_cutoff_tbl_alt,
  file = file.path(SAVE_DIR_alt, "supp_table_late_day_cutoffs.tsv")
)










####### Fig 9: Spearman correlation between organs (kidney split into late-shed vs late non-shed)


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(corrplot)
})



SAVE_DIR <- "../plots/urine_tissue_manuscript5"
dir.create(SAVE_DIR, showWarnings = FALSE, recursive = TRUE)


kidney_late_shed_all <- tissue_bc_levels %>%
  filter(is_kidney(organ)) %>%
  inner_join(present_lastK, by = c("animal","barcode")) %>%
  distinct(animal, barcode)

kidney_late_non_shed_all <- tissue_bc_levels %>%
  filter(is_kidney(organ)) %>%
  inner_join(absent_lastK, by = c("animal","barcode")) %>%
  distinct(animal, barcode)


compute_spearman_correlation_splitKidney <- function(data,
                                                     animal_id,
                                                     kidney_late_shed_tbl,
                                                     kidney_late_non_shed_tbl) {
  shed_bc <- kidney_late_shed_tbl %>%
    filter(animal == animal_id) %>%
    distinct(barcode)
  
  nonshed_bc <- kidney_late_non_shed_tbl %>%
    filter(animal == animal_id) %>%
    distinct(barcode)
  
  bc_union <- bind_rows(shed_bc, nonshed_bc) %>%
    distinct(barcode)
  
  if (nrow(bc_union) == 0) {
    stop(paste("No kidney late-shed or late non-shed barcodes for animal", animal_id))
  }
  
  animal_data <- data %>%
    filter(animal == animal_id) %>%
    semi_join(bc_union, by = "barcode") %>%
    mutate(
      organ = case_when(
        is_kidney(organ) & barcode %in% shed_bc$barcode    ~ "Kidney LS",
        is_kidney(organ) & barcode %in% nonshed_bc$barcode ~ "Kidney LNS",
        TRUE                                               ~ organ
      )
    ) %>%
    select(organ, barcode, frac_weighted_count) %>%
    pivot_wider(
      names_from  = organ,
      values_from = frac_weighted_count,
      values_fill = 0
    ) %>%
    column_to_rownames("barcode")
  
  if (ncol(animal_data) < 2) {
    stop(paste("Not enough organs after splitting kidney for animal", animal_id))
  }
  
  cor_matrix <- cor(
    as.matrix(animal_data),
    method = "spearman",
    use    = "pairwise.complete.obs"
  )
  
  
  colnames(cor_matrix) <- gsub(" ", "\n", colnames(cor_matrix))
  rownames(cor_matrix) <- gsub(" ", "\n", rownames(cor_matrix))
  
  cor_matrix
}

generate_correlation_plot <- function(cor_matrix, animal_color) {
  corrplot.mixed(
    cor_matrix,
    upper      = "circle",
    lower      = "number",
    tl.col     = "black",
    tl.font    = 2,
    tl.cex     = 0.6,
    tl.srt     = 45,
    number.cex = 0.5,
    lower.col  = "black",
    is.corr    = TRUE,
    upper.col  = colorRampPalette(c("gray", "white", animal_color))(400),
    col.lim    = c(-1, 1)
  )
}




out_pdf <- file.path(SAVE_DIR, "Fig8_kidneySplit_spearmanCorrelation_between_organs.pdf")

pdf(out_pdf, onefile = TRUE, height = 11, width = 11)


par(
  mfrow = c(2, 2),
  mar  = c(0, 0, 1.8, 0),   # tighter per-panel margins
  oma  = c(0, 0, 4, 0.0)    # outer margin, extra on top for title
)

for (animal_id in names(animal_colors)) {
  cor_matrix <- compute_spearman_correlation_splitKidney(
    data                     = tissue_bc_levels,
    animal_id                = animal_id,
    kidney_late_shed_tbl     = kidney_late_shed_all,
    kidney_late_non_shed_tbl = kidney_late_non_shed_all
  )
  generate_correlation_plot(cor_matrix, animal_colors[[animal_id]])
  title(main = animal_id, line = 1.0, cex.main = 1.15)
}


mtext(
  "Spearman correlation between organs\nKidney split: late-shed (Kidney LS) and late non-shed (Kidney LNS)",
  outer = TRUE, side = 3, line = 1.0, cex = 1.2, font = 2
)

par(mfrow = c(1, 1))
dev.off()



























###### Fig 9B: Spearman correlation between mice (kidney-only)

suppressPackageStartupMessages({
  library(corrplot)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

SAVE_DIR <- "../plots/urine_tissue_manuscript5"
dir.create(SAVE_DIR, showWarnings = FALSE, recursive = TRUE)


COL_LS  <- "#F8766D"  
COL_LNS <- "#00BFC4"  

compute_mouse_corr_for_group <- function(tissue_tbl, group_tbl) {
  kidney_dat <- tissue_tbl %>%
    filter(is_kidney(organ)) %>%
    inner_join(group_tbl, by = c("animal", "barcode")) %>%
    select(animal, barcode, frac_weighted_count)
  
  if (nrow(kidney_dat) == 0) stop("No data for this group.")
  
  wide <- kidney_dat %>%
    pivot_wider(
      names_from  = animal,
      values_from = frac_weighted_count,
      values_fill = 0
    ) %>%
    column_to_rownames("barcode")
  
  if (ncol(wide) < 2) stop("Not enough mice for correlation.")
  
  cor(as.matrix(wide), method = "spearman", use = "pairwise.complete.obs")
}

out_pdf <- file.path(SAVE_DIR, "FigS1_kidney_between_mice_corr_LS_vs_LNS_samePage.pdf")
pdf(out_pdf, width = 9, height = 5, onefile = TRUE)

par(
  mfrow = c(1, 2),
  mar  = c(1.5, 1.5, 3.2, .8),
  oma  = c(0.8, 0.8, 2.2, 0.4)
)


cor_kidney_ls_mice <- compute_mouse_corr_for_group(
  tissue_tbl = tissue_bc_levels,
  group_tbl  = kidney_present_all
)

corrplot.mixed(
  cor_kidney_ls_mice,
  upper       = "circle",
  lower       = "number",
  tl.col      = "black",
  tl.cex      = 0.95,
  number.cex  = 0.85,
  lower.col   = "black",
  is.corr     = TRUE,
  upper.col   = colorRampPalette(c("gray", "white", COL_LS))(400),
  col.lim     = c(-1, 1),
  addgrid.col = "gray90",
  addgrid.lwd = 0.3
)
title(main = "Kidney late-shed", line = 1.4, cex.main = .8)

cor_kidney_lns_mice <- compute_mouse_corr_for_group(
  tissue_tbl = tissue_bc_levels,
  group_tbl  = kidney_notLastK_all
)

corrplot.mixed(
  cor_kidney_lns_mice,
  upper       = "circle",
  lower       = "number",
  tl.col      = "black",
  tl.cex      = 0.95,
  number.cex  = 0.85,
  lower.col   = "black",
  is.corr     = TRUE,
  upper.col   = colorRampPalette(c("gray", "white", COL_LNS))(400),
  col.lim     = c(-1, 1),
  addgrid.col = "gray90",
  addgrid.lwd = 0.3
)
title(main = "Kidney late non-shed", line = 1.4, cex.main = .8)


mtext(
  "Spearman correlation between mice (kidney barcodes)",
  outer = TRUE, side = 3, line = 0.6, cex = 1.2, font = 2
)

dev.off()
par(mfrow = c(1, 1))














# Fig S2: GC content of kidney barcodes


suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(patchwork)
})

SAVE_DIR <- "../plots/urine_tissue_manuscript5"
dir.create(SAVE_DIR, showWarnings = FALSE, recursive = TRUE)


gc_pct_all_barcodes <- tissue_bc_levels %>%
  distinct(barcode) %>%
  mutate(
    GC_pct = str_count(barcode, "[GC]") / str_length(barcode),
    animal = "All barcodes"
  )


gc_kidney_late_shed_top50 <- kidney_rank %>%
  group_by(animal) %>%
  arrange(kidney_rank, .by_group = TRUE) %>%
  slice_head(n = 50) %>%
  ungroup() %>%
  distinct(animal, barcode) %>%
  mutate(GC_pct = str_count(barcode, "[GC]") / str_length(barcode))

gc_late_shed_df <- gc_kidney_late_shed_top50 %>%
  bind_rows(gc_pct_all_barcodes) %>%
  mutate(animal = fct_relevel(animal, c("FL", "FR", "ML", "MR", "All barcodes")))

p_gc_late_shed <- ggplot(gc_late_shed_df, aes(animal, GC_pct)) +
  geom_boxplot(aes(color = animal), linewidth = 0.8, show.legend = FALSE) +
  labs(
    title = "Late-shed kidney barcodes (Top 50) vs all tissue barcodes",
    x = "",
    y = "GC content"
  ) +
  scale_color_manual(values = c(animal_colors, `All barcodes` = "gray70")) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  theme_minimal(base_size = 10)


gc_kidney_late_non_shed_top50 <- kidney_rank_notLastK %>%
  group_by(animal) %>%
  arrange(kidney_rank, .by_group = TRUE) %>%
  slice_head(n = 50) %>%
  ungroup() %>%
  distinct(animal, barcode) %>%
  mutate(GC_pct = str_count(barcode, "[GC]") / str_length(barcode))

gc_late_non_shed_df <- gc_kidney_late_non_shed_top50 %>%
  bind_rows(gc_pct_all_barcodes) %>%
  mutate(animal = fct_relevel(animal, c("FL", "FR", "ML", "MR", "All barcodes")))

p_gc_late_non_shed <- ggplot(gc_late_non_shed_df, aes(animal, GC_pct)) +
  geom_boxplot(aes(color = animal), linewidth = 0.8, show.legend = FALSE) +
  labs(
    title = "Late non-shed kidney barcodes (Top 50) vs all tissue barcodes",
    x = "",
    y = "GC content"
  ) +
  scale_color_manual(values = c(animal_colors, `All barcodes` = "gray70")) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  theme_minimal(base_size = 10)


figS2_panel <- p_gc_late_shed | p_gc_late_non_shed +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 12),
    plot.tag.position = c(0.0, 1.02)
  )

print(figS2_panel)

ggsave(
  filename = "FigS2_GC_content_kidney_late_shed_vs_non_shed.pdf",
  plot     = figS2_panel,
  path     = SAVE_DIR,
  width    = 11,
  height   = 6
)








# Fig S3: miRNA-hit enrichment in kidney late-shed vs late non-shed

'''
########miRNA-mRNA match

RNAhybrid -s 3utr_human -t kidney_barcodes_unique.fa -q mmu.mirna.cleaned.fa -b 1 -m 10000 > rnahybrid.kidney_mirna.txt


'''
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

SAVE_DIR <- "../plots/urine_tissue_manuscript5"
dir.create(SAVE_DIR, showWarnings = FALSE, recursive = TRUE)

LAST_K <- 2


rnahybrid_tsv_mouse <- "/stor/work/Sullivan/anik/barcode_project/data/mirna_target_kidney/analysis.rnahybrid/rev_rnahybrid/true_rev_rnahybrid.kidney_mirna.parsed.tsv"
rnahybrid_tsv_muPyV <- "/stor/work/Sullivan/anik/barcode_project/data/mirna_target_kidney/analysis.rnahybrid/rev_rnahybrid/true.muPyV_rev_rnahybrid.kidney_mirna.parsed.tsv"

fasta_with_rev <- "/stor/work/Sullivan/anik/barcode_project/data/mirna_target_kidney/analysis.rnahybrid/rev_rnahybrid/kidney_barcodes_unique_with_rev.fa"

is_kidney <- function(x) grepl("kidney|renal", x, ignore.case = TRUE)


fill_colors <- c(
  "miRNA-hit"      = "#3C8DBC",   # blue-ish
  "Random control" = "#7F7F7F"    # gray
)


read_fasta_named <- function(path) {
  lines <- readLines(path)
  hdr_idx <- grep("^>", lines)
  ends <- c(hdr_idx[-1] - 1, length(lines))
  tibble(
    name = sub("^>", "", lines[hdr_idx]),
    seq  = purrr::map2_chr(hdr_idx + 1, ends, ~ gsub("\\s+", "", paste(lines[.x:.y], collapse = "")))
  )
}

panel_title_theme <- theme(
  plot.title = element_text(face = "plain", size = 10, hjust = 0.5),
  plot.title.position = "plot"
)


last_k_days <- urine_bc_levels %>%
  distinct(animal, Days_pi) %>%
  group_by(animal) %>%
  arrange(Days_pi, .by_group = TRUE) %>%
  slice_tail(n = LAST_K) %>%
  mutate(in_lastK = TRUE) %>%
  ungroup()

present_lastK <- urine_bc_levels %>%
  inner_join(last_k_days, by = c("animal","Days_pi")) %>%
  group_by(animal, barcode) %>%
  summarise(present_lastK = any(bc_level > 0, na.rm = TRUE), .groups = "drop") %>%
  filter(present_lastK) %>%
  select(animal, barcode)

kidney_barcodes_all <- tissue_bc_levels %>%
  filter(is_kidney(organ)) %>%
  distinct(animal, barcode)

kidney_late_shed_set <- kidney_barcodes_all %>%
  semi_join(present_lastK, by = c("animal","barcode")) %>%
  mutate(kidney_status = "Late-shed")

kidney_late_non_shed_set <- kidney_barcodes_all %>%
  anti_join(present_lastK, by = c("animal","barcode")) %>%
  mutate(kidney_status = "Late non-shed")

kidney_status_tbl <- bind_rows(kidney_late_shed_set, kidney_late_non_shed_set) %>%
  mutate(kidney_status = factor(kidney_status, levels = c("Late-shed", "Late non-shed")))


fasta_tbl <- read_fasta_named(fasta_with_rev)

forward_seq_tbl <- fasta_tbl %>%
  mutate(forward_name = sub("-rev$", "", name)) %>%
  group_by(forward_name) %>%
  summarise(forward_seq = seq[which.min(grepl("-rev$", name))], .groups = "drop")


make_mirna_enrichment_plot <- function(rnahybrid_tsv, panel_title, seed = 42) {
  
  miRNA_interaction <- readr::read_tsv(rnahybrid_tsv, show_col_types = FALSE) %>%
    mutate(mirna_number = as.integer(stringr::str_extract(mirna, "\\d+"))) %>%
    filter(!is.na(mirna_number), mirna_number >= 1, mirna_number <= 800) %>%
    mutate(
      barcode_id_raw  = target,
      barcode_id_norm = sub("-rev$", "", barcode_id_raw)
    ) %>%
    left_join(forward_seq_tbl, by = c("barcode_id_norm" = "forward_name"))
  
  
  mirna_hit_barcodes <- miRNA_interaction %>%
    distinct(forward_seq) %>%
    rename(barcode = forward_seq) %>%
    filter(!is.na(barcode), barcode != "")
  
  mirna_hit_in_kidney <- kidney_status_tbl %>%
    semi_join(mirna_hit_barcodes, by = "barcode")
  
  
  set.seed(seed)
  hits_per_mouse <- mirna_hit_in_kidney %>%
    count(animal, name = "n_hits")
  
  rand_control_by_mouse <- kidney_barcodes_all %>%
    group_split(animal) %>%
    purrr::map_dfr(function(df_mouse) {
      this_animal <- df_mouse$animal[1]
      n_hits <- hits_per_mouse %>% filter(animal == this_animal) %>% pull(n_hits)
      n_hits <- ifelse(length(n_hits) == 0 || is.na(n_hits) || n_hits <= 0, 0, n_hits)
      
      if (n_hits == 0) return(tibble(animal = this_animal, barcode = character(0)))
      
      df_mouse %>%
        slice_sample(n = min(n_hits, nrow(df_mouse))) %>%
        select(animal, barcode)
    })
  
  rand_in_kidney_status <- kidney_status_tbl %>%
    semi_join(rand_control_by_mouse, by = c("animal","barcode"))
  
  summarize_status <- function(tbl, label) {
    tbl %>%
      count(animal, kidney_status, name = "n") %>%
      group_by(animal) %>%
      mutate(
        n_total = sum(n),
        prop    = ifelse(n_total > 0, n / n_total, NA_real_)
      ) %>%
      ungroup() %>%
      mutate(group = label)
  }
  
  hit_summary  <- summarize_status(mirna_hit_in_kidney, "miRNA-hit")
  rand_summary <- summarize_status(rand_in_kidney_status, "Random control")
  
  plot_df <- bind_rows(hit_summary, rand_summary) %>%
    mutate(
      group = factor(group, levels = c("Random control", "miRNA-hit"))
    )
  
  ggplot(plot_df, aes(x = kidney_status, y = prop, fill = group)) +
    geom_col(position = position_dodge(width = 0.6), width = 0.55) +
    facet_wrap(vars(animal), nrow = 1) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    scale_fill_manual(values = fill_colors, drop = FALSE) +
    labs(
      title = panel_title,
      x = "",
      y = "Proportion of kidney barcodes",
      fill = ""
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      panel.grid.major.x = element_blank(),
      legend.position = "top"
    ) +
    panel_title_theme
}


p_A <- make_mirna_enrichment_plot(
  rnahybrid_tsv = rnahybrid_tsv_mouse,
  panel_title   = "Murine miRNAs: kidney status of miRNA-hit barcodes vs random control"
)

p_B <- make_mirna_enrichment_plot(
  rnahybrid_tsv = rnahybrid_tsv_muPyV,
  panel_title   = "muPyV miRNA: kidney status of miRNA-hit barcodes vs random control"
)

figS3_panel <- p_A / p_B +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 12),
    plot.tag.position = c(0.0, 1.02)
  )

print(figS3_panel)

ggsave(
  filename = "FigS3_miRNA_hit_enrichment_twoPanels.pdf",
  plot     = figS3_panel,
  path     = SAVE_DIR,
  width    = 11,
  height   = 8
)









# Fig S4: Barcode length distribution (Top-50 kidney barcodes)


suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
})


TOP_N_KIDNEY <- if (exists("TOP_N_KIDNEY")) TOP_N_KIDNEY else if (exists("TOP_N")) TOP_N else 50
SAVE_DIR     <- "../plots/urine_tissue_manuscript5"
dir.create(SAVE_DIR, showWarnings = FALSE, recursive = TRUE)

LABEL_SHED    <- "Late-shed"
LABEL_NONSHED <- "Late non-shed"

fill_cols <- c(
  "Late-shed"     = "#F8766D",  # ggplot default coral
  "Late non-shed" = "#00BFC4"   # ggplot default teal
)

topN_length_dist <- bind_rows(
  kidney_present_topN %>%
    mutate(kidney_status = LABEL_SHED),
  kidney_absent_topN %>%
    mutate(kidney_status = LABEL_NONSHED)
) %>%
  mutate(
    kidney_status  = factor(kidney_status, levels = c(LABEL_SHED, LABEL_NONSHED)),
    barcode_length = str_length(barcode)
  ) %>%
  count(animal, kidney_status, barcode_length, name = "n")



topN_length_dist <- bind_rows(
  kidney_present_topN %>% mutate(kidney_status = LABEL_SHED),
  kidney_absent_topN  %>% mutate(kidney_status = LABEL_NONSHED)
) %>%
  mutate(
    kidney_status  = factor(kidney_status, levels = c(LABEL_SHED, LABEL_NONSHED)),
    barcode_length = str_length(barcode)
  ) %>%
  count(animal, kidney_status, barcode_length, name = "n") %>%
  tidyr::complete(
    animal,
    barcode_length,
    kidney_status,
    fill = list(n = 0)
  )


p_len_dist <- ggplot(
  topN_length_dist,
  aes(x = barcode_length, y = n, fill = kidney_status)
) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  facet_wrap(vars(animal), nrow = 1) +
  scale_x_continuous(breaks = sort(unique(topN_length_dist$barcode_length))) +
  scale_fill_manual(values = fill_cols) +
  labs(
    title = sprintf("Barcode length distribution of Top-%d kidney barcodes by group", TOP_N_KIDNEY),
    x = "Barcode length (nt)",
    y = "Number of barcodes",
    fill = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )

print(p_len_dist)


ggsave(
  filename = paste0("FigS4_barcode_length_distribution_Top", TOP_N_KIDNEY, ".pdf"),
  plot     = p_len_dist,
  path     = SAVE_DIR,
  width    = 7.24,
  height   = 3.75
)






if (!exists("SAVE_DIR_LOCAL")) SAVE_DIR_LOCAL <- if (exists("SAVE_DIR")) SAVE_DIR else "."
dir.create(SAVE_DIR_LOCAL, showWarnings = FALSE, recursive = TRUE)

.write_xlsx <- function(df, out_path, sheet = "Sheet1") {
  if (!is.data.frame(df)) df <- as.data.frame(df)
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeDataTable(wb, sheet, df)
    openxlsx::saveWorkbook(wb, out_path, overwrite = TRUE)
  } else if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(list(!!sheet := df), out_path)
  } else {
    stop("Need either {openxlsx} or {writexl} installed to write .xlsx")
  }
}

.save_tbl <- function(obj, filename, sheet = "data") {
  out_path <- file.path(SAVE_DIR_LOCAL, filename)
  .write_xlsx(obj, out_path, sheet = sheet)
  message("✅ Saved: ", out_path)
}


if (exists("urine_bc_levels")) .save_tbl(urine_bc_levels, "TableS1_urine_bc_levels.xlsx")
if (exists("tissue_bc_levels")) .save_tbl(tissue_bc_levels, "TableS2_tissue_bc_levels.xlsx")
if (exists("cutoff_99pct_stock_barcodes")) .save_tbl(cutoff_99pct_stock_barcodes, "TableS3_cutoff_99pct_stock_barcodes.xlsx")


if (exists("kidney_counts")) .save_tbl(kidney_counts, "TableS4_Fig1A_kidney_barcode_counts_by_status.xlsx")
if (exists("kidney_total_by_status")) .save_tbl(kidney_total_by_status, "TableS5_Fig1B_total_kidney_barcode_level_by_status.xlsx")
if (exists("kidney_bc_levels_labeled")) .save_tbl(kidney_bc_levels_labeled, "TableS6_Fig1C_kidney_barcode_level_distributions.xlsx")


.make_donut_plot_df <- function(bc_table, highlights) {
  valid_highlights <- highlights[highlights %in% bc_table$barcode]
  if (length(valid_highlights) == 0) return(NULL)
  
  bc_table %>%
    group_by(sample) %>%
    arrange(desc(frac_weighted_count), .by_group = TRUE) %>%
    mutate(
      ymax = cumsum(frac_weighted_count),
      ymin = c(0, head(ymax, n = -1)),
      bc_color = ifelse(barcode %in% valid_highlights,
                        barcode,
                        ifelse(row_number() %% 2 == 0, "background1", "background2"))
    ) %>%
    ungroup()
}

if (exists("winner_barcodes_by_animal") && exists("urine_bc_levels")) {
  donut_data_ls <- purrr::map_dfr(animals, function(a) {
    bc_table <- urine_bc_levels %>% dplyr::filter(animal == a)
    hl <- winner_barcodes_by_animal[[a]]
    dd <- .make_donut_plot_df(bc_table, hl)
    if (is.null(dd)) return(NULL)
    dd %>% mutate(panel = "Fig2A_LateShed", animal = a)
  })
  if (nrow(donut_data_ls) > 0) .save_tbl(donut_data_ls, "TableS7_Fig2A_donut_plot_data_late_shed.xlsx")
}

if (exists("winner_barcodes_by_animal_notLastK") && exists("urine_bc_levels")) {
  donut_data_lns <- purrr::map_dfr(animals, function(a) {
    bc_table <- urine_bc_levels %>% dplyr::filter(animal == a)
    hl <- winner_barcodes_by_animal_notLastK[[a]]
    hl2 <- intersect(hl, bc_table$barcode)
    dd <- .make_donut_plot_df(bc_table, hl2)
    if (is.null(dd)) return(NULL)
    dd %>% mutate(panel = "Fig2B_LateNonShed", animal = a)
  })
  if (nrow(donut_data_lns) > 0) .save_tbl(donut_data_lns, "TableS8_Fig2B_donut_plot_data_late_non_shed.xlsx")
}


if (exists("sum_selected")) .save_tbl(sum_selected, "TableS9_Fig2C_sum_overlay_late_shed.xlsx")
if (exists("sum_selected_notLastK")) .save_tbl(sum_selected_notLastK, "TableS10_Fig2D_sum_overlay_late_non_shed.xlsx")


if (exists("selected_late_shed_tbl") && exists("selected_late_nonshed_tbl") && exists("rand_tbl")) {
  sel_sets <- bind_rows(
    selected_late_shed_tbl %>% mutate(cohort = "Late-shed"),
    selected_late_nonshed_tbl %>% mutate(cohort = "Late non-shed"),
    rand_tbl %>% mutate(cohort = "Random control")
  )
  .save_tbl(sel_sets, "TableS11_Fig3_selected_barcode_sets.xlsx")
}

if (exists("urine_bc_levels") && exists("selected_late_shed_tbl") &&
    exists("selected_late_nonshed_tbl") && exists("rand_tbl") && exists("CAP_VAL")) {
  
  .make_traj_df <- function(selected_tbl, cohort_label, cap) {
    urine_bc_levels %>%
      semi_join(selected_tbl %>% distinct(animal, barcode), by = c("animal","barcode")) %>%
      group_by(animal, barcode, Days_pi) %>%
      summarise(bc_day_level = sum(bc_level, na.rm = TRUE), .groups = "drop") %>%
      mutate(
        bc_capped = pmin(bc_day_level, cap),
        cohort = cohort_label,
        cap_val = cap
      )
  }
  
  traj_df <- bind_rows(
    .make_traj_df(selected_late_shed_tbl,    "Late-shed",     CAP_VAL),
    .make_traj_df(selected_late_nonshed_tbl, "Late non-shed", CAP_VAL),
    .make_traj_df(rand_tbl,                 "Random control",CAP_VAL)
  )
  
  .save_tbl(traj_df, paste0("TableS12_Fig3_urine_trajectories_capped_cap", CAP_VAL, ".xlsx"))
}


if (exists("kidney_present_with_stock")) .save_tbl(kidney_present_with_stock, "TableS13_Fig4_late_shed_stock_rank_scatter_data.xlsx")
if (exists("kidney_notLastK_with_stock")) .save_tbl(kidney_notLastK_with_stock, "TableS14_Fig4_late_non_shed_stock_rank_scatter_data.xlsx")


if (exists("cohorts_with_max")) .save_tbl(cohorts_with_max, "TableS15_Fig5_peak_urine_level_per_barcode_by_cohort.xlsx")


if (exists("enrich_complete")) .save_tbl(enrich_complete, "TableS16_Fig6_smolderer_enrichment_vs_random.xlsx")


if (exists("all_sets")) {
  fig7_tbl <- all_sets %>% filter(!is.na(pct_detect_early))
  .save_tbl(fig7_tbl, "TableS17_Fig7_early_detection_frequency_by_group.xlsx")
}


if (exists("gc_late_shed_df")) .save_tbl(gc_late_shed_df, "TableS18_FigS2A_GC_late_shed_vs_all.xlsx")
if (exists("gc_late_non_shed_df")) .save_tbl(gc_late_non_shed_df, "TableS19_FigS2B_GC_late_non_shed_vs_all.xlsx")


.make_mirna_plot_df <- function(rnahybrid_tsv, seed = 42) {
  
  
  miRNA_interaction <- readr::read_tsv(rnahybrid_tsv, show_col_types = FALSE) %>%
    mutate(mirna_number = as.integer(stringr::str_extract(mirna, "\\d+"))) %>%
    filter(!is.na(mirna_number), mirna_number >= 1, mirna_number <= 800) %>%
    mutate(barcode_id_norm = sub("-rev$", "", target)) %>%
    left_join(forward_seq_tbl, by = c("barcode_id_norm" = "forward_name"))
  
  mirna_hit_barcodes <- miRNA_interaction %>%
    distinct(forward_seq) %>%
    rename(barcode = forward_seq) %>%
    filter(!is.na(barcode), barcode != "")
  
  mirna_hit_in_kidney <- kidney_status_tbl %>%
    semi_join(mirna_hit_barcodes, by = "barcode")
  
  set.seed(seed)
  hits_per_mouse <- mirna_hit_in_kidney %>% count(animal, name = "n_hits")
  
  rand_control_by_mouse <- kidney_barcodes_all %>%
    group_split(animal) %>%
    purrr::map_dfr(function(df_mouse) {
      this_animal <- df_mouse$animal[1]
      n_hits <- hits_per_mouse %>% filter(animal == this_animal) %>% pull(n_hits)
      n_hits <- ifelse(length(n_hits) == 0 || is.na(n_hits) || n_hits <= 0, 0, n_hits)
      if (n_hits == 0) return(tibble(animal = this_animal, barcode = character(0)))
      df_mouse %>% slice_sample(n = min(n_hits, nrow(df_mouse))) %>% select(animal, barcode)
    })
  
  rand_in_kidney_status <- kidney_status_tbl %>%
    semi_join(rand_control_by_mouse, by = c("animal","barcode"))
  
  summarize_status <- function(tbl, label) {
    tbl %>%
      count(animal, kidney_status, name = "n") %>%
      group_by(animal) %>%
      mutate(n_total = sum(n), prop = ifelse(n_total > 0, n / n_total, NA_real_)) %>%
      ungroup() %>%
      mutate(group = label)
  }
  
  hit_summary  <- summarize_status(mirna_hit_in_kidney, "miRNA-hit")
  rand_summary <- summarize_status(rand_in_kidney_status, "Random control")
  
  bind_rows(hit_summary, rand_summary)
}

if (exists("kidney_status_tbl") && exists("kidney_barcodes_all") && exists("forward_seq_tbl") &&
    exists("rnahybrid_tsv_mouse") && exists("rnahybrid_tsv_muPyV")) {
  
  figS3A_tbl <- .make_mirna_plot_df(rnahybrid_tsv_mouse, seed = 42) %>% mutate(panel = "FigS3A_Murine")
  figS3B_tbl <- .make_mirna_plot_df(rnahybrid_tsv_muPyV, seed = 42) %>% mutate(panel = "FigS3B_muPyV")
  
  .save_tbl(figS3A_tbl, "TableS20_FigS3A_murine_miRNA_hit_enrichment.xlsx")
  .save_tbl(figS3B_tbl, "TableS21_FigS3B_muPyV_miRNA_hit_enrichment.xlsx")
}


if (exists("topN_length_dist")) .save_tbl(topN_length_dist, "TableS22_FigS4_barcode_length_distribution.xlsx")










# Fig 3 - UNCAPPED trajectory table (for Supplementary Table)


fig3_traj_uncapped <- urine_bc_levels %>%
  semi_join(
    bind_rows(
      selected_late_shed_tbl,
      selected_late_nonshed_tbl,
      rand_tbl
    ) %>% distinct(animal, barcode),
    by = c("animal", "barcode")
  ) %>%
  group_by(animal, barcode, Days_pi) %>%
  summarise(
    bc_day_level = sum(bc_level, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    group = case_when(
      barcode %in% selected_late_shed_tbl$barcode    ~ "Late-shed",
      barcode %in% selected_late_nonshed_tbl$barcode ~ "Late non-shed",
      barcode %in% rand_tbl$barcode                  ~ "Random control",
      TRUE                                           ~ NA_character_
    )
  ) %>%
  arrange(group, animal, barcode, Days_pi)









# ONE EXCEL FILE, MANY SHEETS (all tables)
# Creates: SAVE_DIR_LOCAL/Supplementary_Tables_All.xlsx


if (!exists("SAVE_DIR_LOCAL")) SAVE_DIR_LOCAL <- if (exists("SAVE_DIR")) SAVE_DIR else "."
dir.create(SAVE_DIR_LOCAL, showWarnings = FALSE, recursive = TRUE)

OUT_XLSX <- file.path(SAVE_DIR_LOCAL, "Supplementary_Tables_All.xlsx")

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("For a single multi-sheet .xlsx, please install {openxlsx}: install.packages('openxlsx')")
}

wb <- openxlsx::createWorkbook()

.add_sheet <- function(wb, sheet_name, df) {
  if (is.null(df)) return(invisible(NULL))
  if (!is.data.frame(df)) df <- as.data.frame(df)
  
  
  sheet_name <- substr(sheet_name, 1, 31)
  if (sheet_name %in% names(wb)) {
    k <- 2
    base <- substr(sheet_name, 1, 28)
    while (paste0(base, "_", k) %in% names(wb)) k <- k + 1
    sheet_name <- paste0(base, "_", k)
  }
  
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeDataTable(wb, sheet = sheet_name, x = df)
  invisible(NULL)
}


if (exists("urine_bc_levels")) .add_sheet(wb, "S1_urine_bc_levels", urine_bc_levels)
if (exists("tissue_bc_levels")) .add_sheet(wb, "S2_tissue_bc_levels", tissue_bc_levels)
if (exists("cutoff_99pct_stock_barcodes")) .add_sheet(wb, "S3_stock_cutoff_99pct", cutoff_99pct_stock_barcodes)


if (exists("kidney_counts")) .add_sheet(wb, "S4_Fig1A_counts", kidney_counts)
if (exists("kidney_total_by_status")) .add_sheet(wb, "S5_Fig1B_totals", kidney_total_by_status)
if (exists("kidney_bc_levels_labeled")) .add_sheet(wb, "S6_Fig1C_levels", kidney_bc_levels_labeled)


if (exists("donut_data_ls")) .add_sheet(wb, "S7_Fig2A_donut_LS", donut_data_ls)
if (exists("donut_data_lns")) .add_sheet(wb, "S8_Fig2B_donut_LNS", donut_data_lns)
if (exists("sum_selected")) .add_sheet(wb, "S9_Fig2C_overlay_LS", sum_selected)
if (exists("sum_selected_notLastK")) .add_sheet(wb, "S10_Fig2D_overlay_LNS", sum_selected_notLastK)


if (exists("sel_sets")) .add_sheet(wb, "S11_Fig3_selected_sets", sel_sets)
if (exists("fig3_traj_uncapped")) .add_sheet(wb, "S12_Fig3_trajectories", fig3_traj_uncapped)


if (exists("kidney_present_with_stock")) .add_sheet(wb, "S13_Fig4_LS_stockrank", kidney_present_with_stock)
if (exists("kidney_notLastK_with_stock")) .add_sheet(wb, "S14_Fig4_LNS_stockrank", kidney_notLastK_with_stock)


if (exists("kidney_present_median_stock_rank")) .add_sheet(wb, "S15_Fig4_LS_medianRank", kidney_present_median_stock_rank)
if (exists("kidney_notLastK_median_stock_rank")) .add_sheet(wb, "S16_Fig4_LNS_medianRank", kidney_notLastK_median_stock_rank)


if (exists("cohorts_with_max")) .add_sheet(wb, "S17_Fig5_peak_urine", cohorts_with_max)


if (exists("enrich_complete")) .add_sheet(wb, "S18_Fig6_enrichment", enrich_complete)


# if (exists("all_sets")) {
#   fig7_tbl <- all_sets %>% dplyr::filter(!is.na(pct_detect_early))
#   .add_sheet(wb, "S19_Fig7_early_detect", fig7_tbl)
# }
if (exists("all_sets_alt")) {
  fig7_alt_tbl <- all_sets_alt %>%
    dplyr::filter(!is.na(pct_detect_early_alt))
  
  .add_sheet(
    wb,
    "S19_Fig7_early_detect_dayCutoff",
    fig7_alt_tbl
  )
}


if (exists("gc_late_shed_df")) .add_sheet(wb, "S20_FigS2A_GC_LS", gc_late_shed_df)
if (exists("gc_late_non_shed_df")) .add_sheet(wb, "S21_FigS2B_GC_LNS", gc_late_non_shed_df)


if (exists("figS3A_tbl")) .add_sheet(wb, "S22_FigS3A_miRNA_mur", figS3A_tbl)
if (exists("figS3B_tbl")) .add_sheet(wb, "S23_FigS3B_miRNA_pyv", figS3B_tbl)


if (exists("topN_length_dist")) .add_sheet(wb, "S24_FigS4_lengthdist", topN_length_dist)


openxlsx::saveWorkbook(wb, OUT_XLSX, overwrite = TRUE)





























# =========================================================
# Fig S3: miRNA-hit enrichment in kidney late-shed vs late non-shed
# Panel A: Murine miRNAs (true_rev_rnahybrid...)
# Panel B: muPyV miRNA  (true.muPyV_rev_rnahybrid...)
# Per-mouse (no pooled panel)
# Saves to: ../plots/urine_tissue_manuscript5/FigS3_miRNA_hit_enrichment_twoPanels.pdf
# =========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

SAVE_DIR <- "../plots/urine_tissue_manuscript5"
dir.create(SAVE_DIR, showWarnings = FALSE, recursive = TRUE)

LAST_K <- 2

# ---- Your two RNAhybrid inputs ----
rnahybrid_tsv_mouse <- "/stor/work/Sullivan/anik/barcode_project/data/mirna_target_kidney/analysis.rnahybrid/rev_rnahybrid/true_rev_rnahybrid.kidney_mirna.parsed.tsv"
rnahybrid_tsv_muPyV <- "/stor/work/Sullivan/anik/barcode_project/data/mirna_target_kidney/analysis.rnahybrid/rev_rnahybrid/true.muPyV_rev_rnahybrid.kidney_mirna.parsed.tsv"

fasta_with_rev <- "/stor/work/Sullivan/anik/barcode_project/data/mirna_target_kidney/analysis.rnahybrid/rev_rnahybrid/kidney_barcodes_unique_with_rev.fa"

is_kidney <- function(x) grepl("kidney|renal", x, ignore.case = TRUE)

# ---- fixed colors for miRNA-hit vs random ----
fill_colors <- c(
  "miRNA-hit"      = "#3C8DBC",   # blue-ish
  "Random control" = "#7F7F7F"    # gray
)

# =========================================================
# Helpers
# =========================================================
read_fasta_named <- function(path) {
  lines <- readLines(path)
  hdr_idx <- grep("^>", lines)
  ends <- c(hdr_idx[-1] - 1, length(lines))
  tibble(
    name = sub("^>", "", lines[hdr_idx]),
    seq  = purrr::map2_chr(hdr_idx + 1, ends, ~ gsub("\\s+", "", paste(lines[.x:.y], collapse = "")))
  )
}

panel_title_theme <- theme(
  plot.title = element_text(face = "plain", size = 10, hjust = 0.5),
  plot.title.position = "plot"
)

# =========================================================
# Build kidney status table ONCE (Late-shed vs Late non-shed)
# =========================================================
last_k_days <- urine_bc_levels %>%
  distinct(animal, Days_pi) %>%
  group_by(animal) %>%
  arrange(Days_pi, .by_group = TRUE) %>%
  slice_tail(n = LAST_K) %>%
  mutate(in_lastK = TRUE) %>%
  ungroup()

present_lastK <- urine_bc_levels %>%
  inner_join(last_k_days, by = c("animal","Days_pi")) %>%
  group_by(animal, barcode) %>%
  summarise(present_lastK = any(bc_level > 0, na.rm = TRUE), .groups = "drop") %>%
  filter(present_lastK) %>%
  select(animal, barcode)

kidney_barcodes_all <- tissue_bc_levels %>%
  filter(is_kidney(organ)) %>%
  distinct(animal, barcode)

kidney_late_shed_set <- kidney_barcodes_all %>%
  semi_join(present_lastK, by = c("animal","barcode")) %>%
  mutate(kidney_status = "Late-shed")

kidney_late_non_shed_set <- kidney_barcodes_all %>%
  anti_join(present_lastK, by = c("animal","barcode")) %>%
  mutate(kidney_status = "Late non-shed")

kidney_status_tbl <- bind_rows(kidney_late_shed_set, kidney_late_non_shed_set) %>%
  mutate(kidney_status = factor(kidney_status, levels = c("Late-shed", "Late non-shed")))

# =========================================================
# FASTA mapping table ONCE
# =========================================================
fasta_tbl <- read_fasta_named(fasta_with_rev)

forward_seq_tbl <- fasta_tbl %>%
  mutate(forward_name = sub("-rev$", "", name)) %>%
  group_by(forward_name) %>%
  summarise(forward_seq = seq[which.min(grepl("-rev$", name))], .groups = "drop")

# =========================================================
# Core function: given RNAhybrid TSV -> make per-mouse barplot
# =========================================================
make_mirna_enrichment_plot <- function(rnahybrid_tsv, panel_title, seed = 42) {
  
  miRNA_interaction <- readr::read_tsv(rnahybrid_tsv, show_col_types = FALSE) %>%
    mutate(mirna_number = as.integer(stringr::str_extract(mirna, "\\d+"))) %>%
    filter(!is.na(mirna_number), mirna_number >= 1, mirna_number <= 800) %>%
    mutate(
      barcode_id_raw  = target,
      barcode_id_norm = sub("-rev$", "", barcode_id_raw)
    ) %>%
    left_join(forward_seq_tbl, by = c("barcode_id_norm" = "forward_name"))
  
  # miRNA-hit barcode set = unique forward_seq observed in RNAhybrid hits
  mirna_hit_barcodes <- miRNA_interaction %>%
    distinct(forward_seq) %>%
    rename(barcode = forward_seq) %>%
    filter(!is.na(barcode), barcode != "")
  
  mirna_hit_in_kidney <- kidney_status_tbl %>%
    semi_join(mirna_hit_barcodes, by = "barcode")
  
  # matched random control per mouse (same n as miRNA hits in that mouse)
  set.seed(seed)
  hits_per_mouse <- mirna_hit_in_kidney %>%
    count(animal, name = "n_hits")
  
  rand_control_by_mouse <- kidney_barcodes_all %>%
    group_split(animal) %>%
    purrr::map_dfr(function(df_mouse) {
      this_animal <- df_mouse$animal[1]
      n_hits <- hits_per_mouse %>% filter(animal == this_animal) %>% pull(n_hits)
      n_hits <- ifelse(length(n_hits) == 0 || is.na(n_hits) || n_hits <= 0, 0, n_hits)
      
      if (n_hits == 0) return(tibble(animal = this_animal, barcode = character(0)))
      
      df_mouse %>%
        slice_sample(n = min(n_hits, nrow(df_mouse))) %>%
        select(animal, barcode)
    })
  
  rand_in_kidney_status <- kidney_status_tbl %>%
    semi_join(rand_control_by_mouse, by = c("animal","barcode"))
  
  summarize_status <- function(tbl, label) {
    tbl %>%
      count(animal, kidney_status, name = "n") %>%
      group_by(animal) %>%
      mutate(
        n_total = sum(n),
        prop    = ifelse(n_total > 0, n / n_total, NA_real_)
      ) %>%
      ungroup() %>%
      mutate(group = label)
  }
  
  hit_summary  <- summarize_status(mirna_hit_in_kidney, "miRNA-hit")
  rand_summary <- summarize_status(rand_in_kidney_status, "Random control")
  
  plot_df <- bind_rows(hit_summary, rand_summary) %>%
    mutate(
      group = factor(group, levels = c("Random control", "miRNA-hit"))
    )
  
  ggplot(plot_df, aes(x = kidney_status, y = prop, fill = group)) +
    geom_col(position = position_dodge(width = 0.6), width = 0.55) +
    facet_wrap(vars(animal), nrow = 1) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    scale_fill_manual(values = fill_colors, drop = FALSE) +
    labs(
      title = panel_title,
      x = "",
      y = "Proportion of kidney barcodes",
      fill = ""
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      panel.grid.major.x = element_blank(),
      legend.position = "top"
    ) +
    panel_title_theme
}

# =========================================================
# Build panels
# =========================================================
p_A <- make_mirna_enrichment_plot(
  rnahybrid_tsv = rnahybrid_tsv_mouse,
  panel_title   = "Murine miRNAs: kidney status of miRNA-hit barcodes vs random control"
)

p_B <- make_mirna_enrichment_plot(
  rnahybrid_tsv = rnahybrid_tsv_muPyV,
  panel_title   = "muPyV miRNA: kidney status of miRNA-hit barcodes vs random control"
)

figS3_panel <- p_A / p_B +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 12),
    plot.tag.position = c(0.0, 1.02)
  )

print(figS3_panel)

ggsave(
  filename = "FigS3_miRNA_hit_enrichment_twoPanels.pdf",
  plot     = figS3_panel,
  path     = SAVE_DIR,
  width    = 11,
  height   = 8
)


# ---- INPUT ----
# Your tibble is assumed to be called `cutoff_99pct_stock_barcodes`
# with a column named `barcode`

out_dir  <- "/stor/work/Sullivan/anik/barcode_project/data/mirna_target_kidney/analysis.rnahybrid/analysis_7feb26"
out_fa   <- file.path(out_dir, "cutoff_99pct_stock_barcodes.fa")

# ---- MAKE FASTA ----
con <- file(out_fa, open = "w")

for (i in seq_len(nrow(cutoff_99pct_stock_barcodes))) {
  bc_seq <- cutoff_99pct_stock_barcodes$barcode[i]
  bc_name <- paste0("barcode", i)
  
  writeLines(paste0(">", bc_name), con)
  writeLines(bc_seq, con)
}

close(con)

# ---- CONFIRM ----
cat("FASTA written to:\n", out_fa, "\n")







#######miRNA

library(readr)

dir <- "/stor/work/Sullivan/anik/barcode_project/data/mirna_target_kidney/analysis.rnahybrid/analysis_7feb26"

# mmu
seed6_mmu <- read_tsv(
  file.path(dir, "seedmatch_mmu_seed6_cutoffFile.tsv"),
  show_col_types = FALSE
)

seed7_mmu <- read_tsv(
  file.path(dir, "seedmatch_mmu_seed7_cutoffFile.tsv"),
  show_col_types = FALSE
)

seed8_mmu <- read_tsv(
  file.path(dir, "seedmatch_mmu_seed8_cutoffFile.tsv"),
  show_col_types = FALSE
)

# muPyV
seed6_muPyV <- read_tsv(
  file.path(dir, "seedmatch_muPyV_seed6_cutoffFile.tsv"),
  show_col_types = FALSE
)

seed7_muPyV <- read_tsv(
  file.path(dir, "seedmatch_muPyV_seed7_cutoffFile.tsv"),
  show_col_types = FALSE
)

seed8_muPyV <- read_tsv(
  file.path(dir, "seedmatch_muPyV_seed8_cutoffFile.tsv"),
  show_col_types = FALSE
)




#Get the set of barcodes from the seed table


seed_barcodes <- seed7_mmu %>%
  distinct(original_barcode_seq)

#2⃣ Find matching rows in urine

urine_hits <- urine_bc_levels %>%
  filter(bc_level > 0) %>%
  inner_join(
    seed_barcodes,
    by = c("barcode" = "original_barcode_seq")
  )



#3 Find matching rows in tissue
tissue_hits_kidney <- tissue_bc_levels %>%
  filter(
    organ == "Kidney",
    bc_level > 0
  ) %>%
  inner_join(
    seed_barcodes,
    by = c("barcode" = "original_barcode_seq")
  )



View(tissue_hits_kidney)

View(urine_hits)



library(dplyr)


# ============================================================
# Seed6 miRNA hit analysis (mouse-by-mouse)
# - Late-shed vs Late non-shed kidney barcodes
# - Unweighted (presence/absence per barcode)
# - Abundance-weighted (kidney bc_level)
# ============================================================

library(dplyr)

# ------------------------------------------------------------
# 1) Define seed6 hit sets (>=1 hit per barcode)
# ------------------------------------------------------------

hits_muPyV_seed6 <- seed6_muPyV %>%
  distinct(original_barcode_seq) %>%
  rename(barcode = original_barcode_seq) %>%
  mutate(hit_muPyV_seed6 = TRUE)

hits_mmu_seed6 <- seed6_mmu %>%
  distinct(original_barcode_seq) %>%
  rename(barcode = original_barcode_seq) %>%
  mutate(hit_mmu_seed6 = TRUE)

# ------------------------------------------------------------
# 2) Unweighted: fraction of barcodes with >=1 seed6 hit
#    (per mouse, per Late-shed / Late non-shed group)
# ------------------------------------------------------------

presence_by_mouse_muPyV_seed6 <- kidney_classified %>%
  left_join(hits_muPyV_seed6, by = "barcode") %>%
  mutate(hit_muPyV_seed6 = if_else(is.na(hit_muPyV_seed6), FALSE, hit_muPyV_seed6)) %>%
  group_by(animal, status) %>%
  summarise(
    n_barcodes = n(),
    n_hit = sum(hit_muPyV_seed6),
    pct_hit = n_hit / n_barcodes,
    .groups = "drop"
  )

presence_by_mouse_mmu_seed6 <- kidney_classified %>%
  left_join(hits_mmu_seed6, by = "barcode") %>%
  mutate(hit_mmu_seed6 = if_else(is.na(hit_mmu_seed6), FALSE, hit_mmu_seed6)) %>%
  group_by(animal, status) %>%
  summarise(
    n_barcodes = n(),
    n_hit = sum(hit_mmu_seed6),
    pct_hit = n_hit / n_barcodes,
    .groups = "drop"
  )

# ------------------------------------------------------------
# 3) Prepare kidney abundance (bc_level) per mouse x barcode
# ------------------------------------------------------------

kidney_abundance <- tissue_bc_levels %>%
  filter(is_kidney(organ), bc_level > 0) %>%
  group_by(animal, barcode) %>%
  summarise(kidney_bc_level = sum(bc_level, na.rm = TRUE), .groups = "drop")

kidney_classified_abund <- kidney_classified %>%
  left_join(kidney_abundance, by = c("animal", "barcode"))

# ------------------------------------------------------------
# 4) Abundance-weighted: fraction of kidney bc_level
#    carried by seed6-hit barcodes (per mouse, per group)
# ------------------------------------------------------------

weighted_by_mouse_muPyV_seed6 <- kidney_classified_abund %>%
  left_join(hits_muPyV_seed6, by = "barcode") %>%
  mutate(
    hit_muPyV_seed6 = if_else(is.na(hit_muPyV_seed6), FALSE, hit_muPyV_seed6),
    kidney_bc_level = if_else(is.na(kidney_bc_level), 0, kidney_bc_level)
  ) %>%
  group_by(animal, status) %>%
  summarise(
    total_kidney_bc_level = sum(kidney_bc_level),
    hit_kidney_bc_level   = sum(kidney_bc_level[hit_muPyV_seed6]),
    pct_kidney_bc_level_hit = hit_kidney_bc_level / total_kidney_bc_level,
    .groups = "drop"
  )

weighted_by_mouse_mmu_seed6 <- kidney_classified_abund %>%
  left_join(hits_mmu_seed6, by = "barcode") %>%
  mutate(
    hit_mmu_seed6 = if_else(is.na(hit_mmu_seed6), FALSE, hit_mmu_seed6),
    kidney_bc_level = if_else(is.na(kidney_bc_level), 0, kidney_bc_level)
  ) %>%
  group_by(animal, status) %>%
  summarise(
    total_kidney_bc_level = sum(kidney_bc_level),
    hit_kidney_bc_level   = sum(kidney_bc_level[hit_mmu_seed6]),
    pct_kidney_bc_level_hit = hit_kidney_bc_level / total_kidney_bc_level,
    .groups = "drop"
  )

# ------------------------------------------------------------
# 5) Outputs now available:
# ------------------------------------------------------------
presence_by_mouse_muPyV_seed6
presence_by_mouse_mmu_seed6
weighted_by_mouse_muPyV_seed6
weighted_by_mouse_mmu_seed6








# ============================================================
# Seed7 miRNA hit analysis (mouse-by-mouse)
# - Late-shed vs Late non-shed kidney barcodes
# - Unweighted (presence/absence per barcode)
# - Abundance-weighted (kidney bc_level)
# ============================================================

library(dplyr)

# ------------------------------------------------------------
# 1) Define seed7 hit sets (>=1 hit per barcode)
# ------------------------------------------------------------

hits_muPyV_seed7 <- seed7_muPyV %>%
  distinct(original_barcode_seq) %>%
  rename(barcode = original_barcode_seq) %>%
  mutate(hit_muPyV_seed7 = TRUE)

hits_mmu_seed7 <- seed7_mmu %>%
  distinct(original_barcode_seq) %>%
  rename(barcode = original_barcode_seq) %>%
  mutate(hit_mmu_seed7 = TRUE)

# ------------------------------------------------------------
# 2) Unweighted: fraction of barcodes with >=1 seed7 hit
#    (per mouse, per Late-shed / Late non-shed group)
# ------------------------------------------------------------

presence_by_mouse_muPyV_seed7 <- kidney_classified %>%
  left_join(hits_muPyV_seed7, by = "barcode") %>%
  mutate(hit_muPyV_seed7 = if_else(is.na(hit_muPyV_seed7), FALSE, hit_muPyV_seed7)) %>%
  group_by(animal, status) %>%
  summarise(
    n_barcodes = n(),
    n_hit = sum(hit_muPyV_seed7),
    pct_hit = n_hit / n_barcodes,
    .groups = "drop"
  )

presence_by_mouse_mmu_seed7 <- kidney_classified %>%
  left_join(hits_mmu_seed7, by = "barcode") %>%
  mutate(hit_mmu_seed7 = if_else(is.na(hit_mmu_seed7), FALSE, hit_mmu_seed7)) %>%
  group_by(animal, status) %>%
  summarise(
    n_barcodes = n(),
    n_hit = sum(hit_mmu_seed7),
    pct_hit = n_hit / n_barcodes,
    .groups = "drop"
  )

# ------------------------------------------------------------
# 3) Prepare kidney abundance (bc_level) per mouse x barcode
#    (reuses kidney_abundance + kidney_classified_abund if already made)
# ------------------------------------------------------------

if (!exists("kidney_abundance")) {
  kidney_abundance <- tissue_bc_levels %>%
    filter(is_kidney(organ), bc_level > 0) %>%
    group_by(animal, barcode) %>%
    summarise(kidney_bc_level = sum(bc_level, na.rm = TRUE), .groups = "drop")
}

if (!exists("kidney_classified_abund")) {
  kidney_classified_abund <- kidney_classified %>%
    left_join(kidney_abundance, by = c("animal", "barcode"))
}

# ------------------------------------------------------------
# 4) Abundance-weighted: fraction of kidney bc_level
#    carried by seed7-hit barcodes (per mouse, per group)
# ------------------------------------------------------------

weighted_by_mouse_muPyV_seed7 <- kidney_classified_abund %>%
  left_join(hits_muPyV_seed7, by = "barcode") %>%
  mutate(
    hit_muPyV_seed7 = if_else(is.na(hit_muPyV_seed7), FALSE, hit_muPyV_seed7),
    kidney_bc_level = if_else(is.na(kidney_bc_level), 0, kidney_bc_level)
  ) %>%
  group_by(animal, status) %>%
  summarise(
    total_kidney_bc_level = sum(kidney_bc_level),
    hit_kidney_bc_level   = sum(kidney_bc_level[hit_muPyV_seed7]),
    pct_kidney_bc_level_hit = hit_kidney_bc_level / total_kidney_bc_level,
    .groups = "drop"
  )

weighted_by_mouse_mmu_seed7 <- kidney_classified_abund %>%
  left_join(hits_mmu_seed7, by = "barcode") %>%
  mutate(
    hit_mmu_seed7 = if_else(is.na(hit_mmu_seed7), FALSE, hit_mmu_seed7),
    kidney_bc_level = if_else(is.na(kidney_bc_level), 0, kidney_bc_level)
  ) %>%
  group_by(animal, status) %>%
  summarise(
    total_kidney_bc_level = sum(kidney_bc_level),
    hit_kidney_bc_level   = sum(kidney_bc_level[hit_mmu_seed7]),
    pct_kidney_bc_level_hit = hit_kidney_bc_level / total_kidney_bc_level,
    .groups = "drop"
  )

# ------------------------------------------------------------
# 5) Outputs now available:
# ------------------------------------------------------------
presence_by_mouse_muPyV_seed7
presence_by_mouse_mmu_seed7
weighted_by_mouse_muPyV_seed7
weighted_by_mouse_mmu_seed7















# ============================================================
# Seed7 host miRNA analysis restricted to kidney-enriched miRNAs
# - Unweighted and abundance-weighted analyses
# - Mouse-by-mouse, Late-shed vs Late non-shed
# ============================================================

library(dplyr)

# ------------------------------------------------------------
# 1) Define kidney-enriched / kidney-relevant murine miRNAs
# ------------------------------------------------------------

kidney_mmu_mirnas <- c(
  "mmu-miR-192-5p",
  "mmu-miR-194-5p",
  "mmu-miR-204-5p",
  "mmu-miR-215-5p",
  "mmu-miR-30a-5p",
  "mmu-miR-30b-5p",
  "mmu-miR-30c-5p",
  "mmu-miR-30d-5p",
  "mmu-miR-30e-5p",
  "mmu-miR-21-5p",
  "mmu-miR-146a-5p",
  "mmu-miR-155-5p"
)


# df <- read.delim("/stor/work/Sullivan/anik/barcode_project/data/mirna_target_kidney/analysis.rnahybrid/analysis_7feb26/PNAS2020_miRNA_atlas/pnas.2002277117.sd02.tsv",
#                  check.names = FALSE,
#                  stringsAsFactors = FALSE)
# 
# names(df)[ncol(df)] <- "biotype"
# 
# mir <- df[df$biotype == "miRNA", , drop = FALSE]
# 
# kidney_cols <- grep("^Kidney_", names(mir), value = TRUE)
# 
# mir$kidney_mean <- rowMeans(mir[, kidney_cols, drop = FALSE], na.rm = TRUE)
# 
# kidney_abundant <- mir[order(mir$kidney_mean, decreasing = TRUE),
#                        c("kidney_mean", kidney_cols)]
# 
# head(kidney_abundant, 50)
# 
# # ============================================================
# # Filter seed7_mmu to only miRNAs in top-50 kidney-abundant atlas list
# # (handles naming differences between Mir10b vs mmu-miR-10b-5p)
# # ============================================================
# 
# library(dplyr)
# library(stringr)
# 
# # ---- 1) Get top 50 kidney-abundant miRNAs from your atlas result ----
# top50_atlas <- rownames(head(kidney_abundant, 50))   # e.g., "Mir10b", "Mir26a.2", ...
# 
# # ---- 2) Make a comparable "key" for atlas names ----
# atlas_to_key <- function(x) {
#   x %>%
#     tolower() %>%
#     str_replace("^mir", "mir") %>%      # keep mir prefix stable
#     str_replace_all("\\.", "-") %>%     # Mir26a.2 -> Mir26a-2
#     str_replace_all("[^a-z0-9]", "")    # remove punctuation: mir26a-2 -> mir26a2
# }
# 
# top50_keys <- atlas_to_key(top50_atlas)
# 
# # ---- 3) Make a comparable "key" for seed7_mmu matched_miRNA ----
# mirbase_to_key <- function(x) {
#   x %>%
#     tolower() %>%
#     str_replace("^mmu-", "") %>%                         # drop species prefix
#     str_replace("-(3p|5p)$", "") %>%                     # drop arm if present at end
#     str_replace_all("[^a-z0-9]", "")                     # remove punctuation
# }
# 
# seed7_mmu_tagged <- seed7_mmu %>%
#   mutate(
#     matched_miRNA_key = mirbase_to_key(matched_miRNA)
#   )
# 
# # ---- 4) Filter seed7_mmu to matches ----
# seed7_mmu_top50kidney <- seed7_mmu_tagged %>%
#   filter(matched_miRNA_key %in% top50_keys)
# 
# # ---- 5) Sanity checks: which atlas miRNAs didn’t match anything? ----
# matched_keys_in_seed <- unique(seed7_mmu_tagged$matched_miRNA_key)
# unmatched_atlas <- top50_atlas[!(top50_keys %in% matched_keys_in_seed)]
# 
# cat("Top50 atlas miRNAs:", length(top50_atlas), "\n")
# cat("Rows retained from seed7_mmu:", nrow(seed7_mmu_top50kidney), "\n")
# cat("Atlas miRNAs with no direct match in seed7_mmu (check naming):\n")
# print(unmatched_atlas)
# 
# # ---- Output object ----
# seed7_mmu_top50kidney
library(readxl)
library(dplyr)

# Load the scaled miRNA file
df <- read_xlsx("/stor/work/Sullivan/anik/barcode_project/data/mirna_target_kidney/analysis.rnahybrid/analysis_7feb26/PNAS2020_miRNA_atlas/pnas.2002277117.sd04.xlsx")

# Filter: only Kidney rows
kidney_miR <- df %>%
  filter(Tissue == "Kidney")

# Select miRNA and mean scaled expression
kidney_ranked <- kidney_miR %>%
  select(miRNA, `Mean scaled expression *`) %>%
  arrange(desc(`Mean scaled expression *`))

# Print top 50 (or remove head() to show all)
head(kidney_ranked, 50)
top50_miRs <- kidney_ranked$miRNA[1:100]

# 5) Filter your seed7_mmu table by these miRNAs
# Adjust "miRNA" column name below if your table uses something like "miRNA_name"
seed7_mmu_kidney <- seed7_mmu %>%
  filter(matched_miRNA %in% top50_miRs)
seed7_mmu_kidney

# ------------------------------------------------------------
# 2) Subset seed7 mmu hits to kidney-relevant miRNAs only
# ------------------------------------------------------------

#seed7_mmu_kidney <- seed7_mmu %>%
#  filter(matched_miRNA %in% kidney_mmu_mirnas)
#seed7_mmu_kidney <- seed7_mmu_top50kidney
# ------------------------------------------------------------
# 3) Define barcode hit set (>=1 kidney-miRNA seed match)
# ------------------------------------------------------------

hits_mmu_seed7_kidney <- seed7_mmu_kidney %>%
  distinct(original_barcode_seq) %>%
  rename(barcode = original_barcode_seq) %>%
  mutate(hit_mmu_seed7_kidney = TRUE)

# ------------------------------------------------------------
# 4) Unweighted analysis:
#    Fraction of kidney barcodes with kidney-miRNA seed matches
# ------------------------------------------------------------

presence_by_mouse_mmu_seed7_kidney <- kidney_classified %>%
  left_join(hits_mmu_seed7_kidney, by = "barcode") %>%
  mutate(
    hit_mmu_seed7_kidney = if_else(
      is.na(hit_mmu_seed7_kidney),
      FALSE,
      hit_mmu_seed7_kidney
    )
  ) %>%
  group_by(animal, status) %>%
  summarise(
    n_barcodes = n(),
    n_hit = sum(hit_mmu_seed7_kidney),
    pct_hit = n_hit / n_barcodes,
    .groups = "drop"
  )

# ------------------------------------------------------------
# 5) Abundance-weighted analysis:
#    Fraction of kidney viral abundance in kidney-miRNA-hit barcodes
# ------------------------------------------------------------

weighted_by_mouse_mmu_seed7_kidney <- kidney_classified_abund %>%
  left_join(hits_mmu_seed7_kidney, by = "barcode") %>%
  mutate(
    hit_mmu_seed7_kidney = if_else(
      is.na(hit_mmu_seed7_kidney),
      FALSE,
      hit_mmu_seed7_kidney
    ),
    kidney_bc_level = if_else(is.na(kidney_bc_level), 0, kidney_bc_level)
  ) %>%
  group_by(animal, status) %>%
  summarise(
    total_kidney_bc_level = sum(kidney_bc_level),
    hit_kidney_bc_level   = sum(kidney_bc_level[hit_mmu_seed7_kidney]),
    pct_kidney_bc_level_hit =
      hit_kidney_bc_level / total_kidney_bc_level,
    .groups = "drop"
  )

# ------------------------------------------------------------
# 6) Outputs
# ------------------------------------------------------------

presence_by_mouse_mmu_seed7_kidney
weighted_by_mouse_mmu_seed7_kidney







library(openxlsx)

# Create a new Excel workbook
wb <- createWorkbook()

# -----------------------------
# Sheet 1: Unweighted muPyV seed7
# -----------------------------
addWorksheet(wb, "muPyV_seed7_unweighted")
writeData(
  wb,
  sheet = "muPyV_seed7_unweighted",
  x = presence_by_mouse_muPyV_seed7,
  startRow = 1,
  startCol = 1
)

# -----------------------------
# Sheet 2: Unweighted mmu seed7
# -----------------------------
addWorksheet(wb, "mmu_seed7_unweighted")
writeData(
  wb,
  sheet = "mmu_seed7_unweighted",
  x = presence_by_mouse_mmu_seed7,
  startRow = 1,
  startCol = 1
)

# -----------------------------
# Sheet 3: Weighted muPyV seed7
# -----------------------------
addWorksheet(wb, "muPyV_seed7_weighted")
writeData(
  wb,
  sheet = "muPyV_seed7_weighted",
  x = weighted_by_mouse_muPyV_seed7,
  startRow = 1,
  startCol = 1
)

# -----------------------------
# Sheet 4: Weighted mmu seed7
# -----------------------------
addWorksheet(wb, "mmu_seed7_weighted")
writeData(
  wb,
  sheet = "mmu_seed7_weighted",
  x = weighted_by_mouse_mmu_seed7,
  startRow = 1,
  startCol = 1
)

saveWorkbook(
  wb,
  file = "Seed7_miRNA_hit_summary_tables.xlsx",
  overwrite = TRUE
)














# ============================================================
# Seed7 plots + unified controls
# Panels:
# A) Unweighted pct_hit per mouse (mmu vs muPyV)
# B) Weighted pct_kidney_bc_level_hit per mouse (mmu vs muPyV)
# C) mmu control: length+GC-matched resampling null (violin behind points)
# D) muPyV control: permutation null for "# hit barcodes in Late-shed"
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)

set.seed(1)

# -----------------------------
# 0) Hit sets (seed7)
# -----------------------------
hits_muPyV_seed7 <- seed7_muPyV %>%
  distinct(original_barcode_seq) %>%
  rename(barcode = original_barcode_seq) %>%
  mutate(hit_muPyV = TRUE)

hits_mmu_seed7 <- seed7_mmu %>%
  distinct(original_barcode_seq) %>%
  rename(barcode = original_barcode_seq) %>%
  mutate(hit_mmu = TRUE)

# -----------------------------
# 1) Kidney abundance per mouse x barcode (bc_level)
# -----------------------------
kidney_abundance <- tissue_bc_levels %>%
  filter(is_kidney(organ), bc_level > 0) %>%
  group_by(animal, barcode) %>%
  summarise(kidney_bc_level = sum(bc_level, na.rm = TRUE), .groups = "drop")

# -----------------------------
# 2) Master table: one row per mouse x barcode in kidney_classified
#    Includes status, kidney_bc_level, hit flags, length + GC
# -----------------------------
barcode_master_seed7 <- kidney_classified %>%
  left_join(kidney_abundance, by = c("animal", "barcode")) %>%
  left_join(hits_mmu_seed7, by = "barcode") %>%
  left_join(hits_muPyV_seed7, by = "barcode") %>%
  mutate(
    kidney_bc_level = if_else(is.na(kidney_bc_level), 0, kidney_bc_level),
    hit_mmu  = if_else(is.na(hit_mmu), FALSE, hit_mmu),
    hit_muPyV = if_else(is.na(hit_muPyV), FALSE, hit_muPyV),
    len = nchar(barcode),
    gc  = (str_count(barcode, "[Gg]") + str_count(barcode, "[Cc]")) / pmax(len, 1)
  )

# -----------------------------
# 3) Panel A + B observed summaries (per mouse, per status)
# -----------------------------
obs_unweighted <- barcode_master_seed7 %>%
  group_by(animal, status) %>%
  summarise(
    pct_hit_mmu   = mean(hit_mmu),
    pct_hit_muPyV = mean(hit_muPyV),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = starts_with("pct_hit_"),
               names_to = "dataset",
               values_to = "pct_hit") %>%
  mutate(dataset = recode(dataset,
                          pct_hit_mmu = "mmu",
                          pct_hit_muPyV = "muPyV"))

obs_weighted <- barcode_master_seed7 %>%
  group_by(animal, status) %>%
  summarise(
    pct_weighted_mmu = if_else(sum(kidney_bc_level) > 0,
                               sum(kidney_bc_level[hit_mmu]) / sum(kidney_bc_level),
                               NA_real_),
    pct_weighted_muPyV = if_else(sum(kidney_bc_level) > 0,
                                 sum(kidney_bc_level[hit_muPyV]) / sum(kidney_bc_level),
                                 NA_real_),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = starts_with("pct_weighted_"),
               names_to = "dataset",
               values_to = "pct_weighted") %>%
  mutate(dataset = recode(dataset,
                          pct_weighted_mmu = "mmu",
                          pct_weighted_muPyV = "muPyV"))

# -----------------------------
# 4) Helper: length+GC binning for matched resampling (mmu control)
#    - length bins: exact length
#    - GC bins: 5 bins within each mouse (quantile-based)
# -----------------------------
master_for_null <- barcode_master_seed7 %>%
  group_by(animal) %>%
  mutate(
    gc_bin = ntile(gc, 5),
    len_bin = len
  ) %>%
  ungroup()

# Stratified sample: for a given mouse & status target,
# sample from same mouse pool matching counts per (len_bin, gc_bin).
stratified_sample <- function(pool_df, target_df, strata = c("len_bin", "gc_bin")) {
  target_counts <- target_df %>%
    count(across(all_of(strata)), name = "n_target")
  
  sampled <- target_counts %>%
    group_split(across(all_of(strata))) %>%
    map_dfr(function(one_stratum) {
      n_take <- one_stratum$n_target[1]
      # values for this stratum
      st_vals <- one_stratum %>% select(all_of(strata)) %>% slice(1)
      
      candidates <- pool_df
      for (s in strata) {
        candidates <- candidates %>% filter(.data[[s]] == st_vals[[s]][1])
      }
      
      # If not enough candidates in this stratum, sample with replacement as fallback
      if (nrow(candidates) == 0) {
        return(pool_df %>% slice_sample(n = n_take, replace = TRUE))
      } else if (nrow(candidates) < n_take) {
        return(candidates %>% slice_sample(n = n_take, replace = TRUE))
      } else {
        return(candidates %>% slice_sample(n = n_take, replace = FALSE))
      }
    })
  
  sampled
}

# -----------------------------
# 5) Panel C: mmu null distributions (length+GC-matched)
#    Build null for BOTH unweighted and weighted metrics
# -----------------------------
N_REP <- 2000  # increase if you want smoother violins

make_mmu_null <- function(animal_id, status_level) {
  pool <- master_for_null %>% filter(animal == animal_id)
  target <- pool %>% filter(status == status_level)
  
  # if target is empty, skip
  if (nrow(target) == 0) return(NULL)
  
  reps <- map_dfr(seq_len(N_REP), function(i) {
    samp <- stratified_sample(pool, target, strata = c("len_bin", "gc_bin"))
    
    tibble(
      animal = animal_id,
      status = status_level,
      rep = i,
      pct_hit = mean(samp$hit_mmu),
      pct_weighted = if_else(sum(samp$kidney_bc_level) > 0,
                             sum(samp$kidney_bc_level[samp$hit_mmu]) / sum(samp$kidney_bc_level),
                             NA_real_)
    )
  })
  reps
}

animals <- sort(unique(master_for_null$animal))
statuses <- levels(master_for_null$status)

mmu_null <- map_dfr(animals, function(a) {
  map_dfr(statuses, function(s) make_mmu_null(a, s))
})

# -----------------------------
# 6) Panel D: muPyV permutation null
#    For each mouse: null distribution of "# muPyV-hit barcodes in Late-shed"
# -----------------------------
N_PERM <- 20000

# muPyV_perm <- map_dfr(animals, function(a) {
#   pool <- barcode_master_seed7 %>% filter(animal == a)
#   target <- pool %>% filter(status == "Late-shed")
#   
#   n_target <- nrow(target)
#   if (n_target == 0) return(NULL)
#   
#   obs_hits <- sum(target$hit_muPyV)
#   
#   null_hits <- replicate(N_PERM, {
#     samp_idx <- sample.int(nrow(pool), size = n_target, replace = FALSE)
#     sum(pool$hit_muPyV[samp_idx])
#   })
#   
#   tibble(
#     animal = a,
#     obs_hits = obs_hits,
#     null_hits = null_hits
#   )
# })
# 
# muPyV_perm_pvals <- muPyV_perm %>%
#   group_by(animal) %>%
#   summarise(
#     obs_hits = first(obs_hits),
#     p_null_le_obs = mean(null_hits <= obs_hits),   # for depletion (obs small)
#     p_null_eq_obs = mean(null_hits == obs_hits),   # probability of exact 0 etc.
#     .groups = "drop"
#   )
# 
# print(muPyV_perm_pvals)

# -----------------------------
# 7) Plotting
# -----------------------------

# Panel A: Unweighted pct_hit (observed)
pA <- ggplot(obs_unweighted, aes(x = status, y = pct_hit, group = animal)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 2) +
  facet_wrap(~ dataset, scales = "free_y") +
  labs(x = NULL, y = "Unweighted fraction of barcodes with seed7 hit") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

# Panel B: Weighted pct_hit (observed)
pB <- ggplot(obs_weighted, aes(x = status, y = pct_weighted, group = animal)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 2) +
  facet_wrap(~ dataset, scales = "free_y") +
  labs(x = NULL, y = "Weighted fraction of kidney bc_level in seed7-hit barcodes") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

# Panel C: mmu control violins behind observed points (unweighted + weighted)
# Unweighted mmu control
pC1 <- ggplot() +
  geom_violin(
    data = mmu_null,
    aes(x = status, y = pct_hit),
    trim = TRUE
  ) +
  geom_point(
    data = obs_unweighted %>% filter(dataset == "mmu"),
    aes(x = status, y = pct_hit, group = animal),
    size = 2
  ) +
  geom_line(
    data = obs_unweighted %>% filter(dataset == "mmu"),
    aes(x = status, y = pct_hit, group = animal),
    alpha = 0.6
  ) +
  facet_wrap(~ animal) +
  labs(x = NULL, y = "mmu: unweighted pct_hit (obs over null)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

# Weighted mmu control
pC2 <- ggplot() +
  geom_violin(
    data = mmu_null,
    aes(x = status, y = pct_weighted),
    trim = TRUE
  ) +
  geom_point(
    data = obs_weighted %>% filter(dataset == "mmu"),
    aes(x = status, y = pct_weighted, group = animal),
    size = 2
  ) +
  geom_line(
    data = obs_weighted %>% filter(dataset == "mmu"),
    aes(x = status, y = pct_weighted, group = animal),
    alpha = 0.6
  ) +
  facet_wrap(~ animal) +
  labs(x = NULL, y = "mmu: weighted pct_hit (obs over null)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

# Panel D: muPyV null distribution (# hits in Late-shed), per mouse
pD <- ggplot(muPyV_perm, aes(x = null_hits)) +
  geom_histogram(bins = 30) +
  geom_vline(
    data = muPyV_perm_pvals,
    aes(xintercept = obs_hits),
    linewidth = 1
  ) +
  facet_wrap(~ animal, scales = "free_y") +
  labs(x = "# muPyV seed7-hit barcodes in Late-shed (null)", y = "Count") +
  theme_classic()



# -----------------------------
# 5b) muPyV null distributions (length+GC-matched)
#     Same as mmu_null, but using hit_muPyV
# -----------------------------
make_muPyV_null <- function(animal_id, status_level) {
  pool <- master_for_null %>% filter(animal == animal_id)
  target <- pool %>% filter(status == status_level)
  
  if (nrow(target) == 0) return(NULL)
  
  reps <- map_dfr(seq_len(N_REP), function(i) {
    samp <- stratified_sample(pool, target, strata = c("len_bin", "gc_bin"))
    
    tibble(
      animal = animal_id,
      status = status_level,
      rep = i,
      pct_hit = mean(samp$hit_muPyV),
      pct_weighted = if_else(sum(samp$kidney_bc_level) > 0,
                             sum(samp$kidney_bc_level[samp$hit_muPyV]) / sum(samp$kidney_bc_level),
                             NA_real_)
    )
  })
  reps
}

muPyV_null <- map_dfr(animals, function(a) {
  map_dfr(statuses, function(s) make_muPyV_null(a, s))
})


# Panel C3: muPyV control violins behind observed points (unweighted)
pC3 <- ggplot() +
  geom_violin(
    data = muPyV_null,
    aes(x = status, y = pct_hit),
    trim = TRUE
  ) +
  geom_point(
    data = obs_unweighted %>% filter(dataset == "muPyV"),
    aes(x = status, y = pct_hit, group = animal),
    size = 2
  ) +
  geom_line(
    data = obs_unweighted %>% filter(dataset == "muPyV"),
    aes(x = status, y = pct_hit, group = animal),
    alpha = 0.6
  ) +
  facet_wrap(~ animal) +
  labs(x = NULL, y = "muPyV: unweighted pct_hit (obs over null)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

# Panel C4: muPyV control violins behind observed points (weighted)
pC4 <- ggplot() +
  geom_violin(
    data = muPyV_null,
    aes(x = status, y = pct_weighted),
    trim = TRUE
  ) +
  geom_point(
    data = obs_weighted %>% filter(dataset == "muPyV"),
    aes(x = status, y = pct_weighted, group = animal),
    size = 2
  ) +
  geom_line(
    data = obs_weighted %>% filter(dataset == "muPyV"),
    aes(x = status, y = pct_weighted, group = animal),
    alpha = 0.6
  ) +
  facet_wrap(~ animal) +
  labs(x = NULL, y = "muPyV: weighted pct_hit (obs over null)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))


# Print plots
print(pA)
print(pB)
print(pC1)
print(pC2)
print(pC3)
print(pC4)

print(pD)

pC1 <- ggplot() +
  geom_violin(
    data = mmu_null,
    aes(x = status, y = pct_hit),
    trim = TRUE
  ) +
  geom_point(
    data = obs_unweighted %>% filter(dataset == "mmu"),
    aes(x = status, y = pct_hit, group = animal),
    size = 2
  ) +
  geom_line(
    data = obs_unweighted %>% filter(dataset == "mmu"),
    aes(x = status, y = pct_hit, group = animal),
    alpha = 0.6
  ) +
  facet_wrap(~ animal) +
  labs(
    x = NULL,
    y = "Fraction of kidney barcodes\nwith host miRNA seed match",
    title = "A. Host miRNA seed targeting in kidney barcodes\n  (Observed vs Null distribution)"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))


pC2 <- ggplot() +
  geom_violin(
    data = mmu_null,
    aes(x = status, y = pct_weighted),
    trim = TRUE
  ) +
  geom_point(
    data = obs_weighted %>% filter(dataset == "mmu"),
    aes(x = status, y = pct_weighted, group = animal),
    size = 2
  ) +
  geom_line(
    data = obs_weighted %>% filter(dataset == "mmu"),
    aes(x = status, y = pct_weighted, group = animal),
    alpha = 0.6
  ) +
  facet_wrap(~ animal) +
  labs(
    x = NULL,
    y = "Fraction of kidney viral abundance in \n host miRNA seed targeted barcodes",
    title = "B. Host miRNA seed targeting weighted by kidney barcode levels\n  (Observed vs Null distribution)"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))


print(pC1)
print(pC2)
library(patchwork)

pC1/pC2


library(patchwork)

ggsave(
  filename = "../plots/urine_tissue_manuscript5/mmu_seed7_observed_vs_null.pdf",
  plot = pC1 / pC2,
  width = 6.3,
  height = 8,
  units = "in"
)





barcode_master_seed7 %>% count(animal, status)
barcode_master_seed7 %>% summarise(min_len=min(len), max_len=max(len))
barcode_master_seed7 %>% group_by(animal) %>% summarise(n_hit_muPyV=sum(hit_muPyV), n_hit_mmu=sum(hit_mmu))




scale_def_colors <- scale_fill_manual(values = c(
  "Late-shed" = "#F8766D",
  "Late non-shed" = "#00BFC4"
))

scale_def_colors_lines <- scale_color_manual(values = c(
  "Late-shed" = "#F8766D",
  "Late non-shed" = "#00BFC4"
))
pC1 <- ggplot() +
  geom_violin(
    data = mmu_null,
    aes(x = status, y = pct_hit, fill = status),
    trim = TRUE,
    alpha = 0.4
  ) +
  geom_point(
    data = obs_unweighted %>% filter(dataset == "mmu"),
    aes(x = status, y = pct_hit, group = animal, color = status),
    size = 2
  ) +
  geom_line(
    data = obs_unweighted %>% filter(dataset == "mmu"),
    aes(x = status, y = pct_hit, group = animal, color = status),
    alpha = 0.6
  ) +
  facet_wrap(~ animal) +
  labs(
    x = NULL,
    y = "Fraction of kidney barcodes\nwith host miRNA seed match",
    title = "A. Host miRNA seed targeting in kidney barcodes\n(Observed vs null distribution)"
  ) +
  scale_def_colors +
  scale_def_colors_lines +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "none"
  )

pC2 <- ggplot() +
  geom_violin(
    data = mmu_null,
    aes(x = status, y = pct_weighted, fill = status),
    trim = TRUE,
    alpha = 0.4
  ) +
  geom_point(
    data = obs_weighted %>% filter(dataset == "mmu"),
    aes(x = status, y = pct_weighted, group = animal, color = status),
    size = 2
  ) +
  geom_line(
    data = obs_weighted %>% filter(dataset == "mmu"),
    aes(x = status, y = pct_weighted, group = animal, color = status),
    alpha = 0.6
  ) +
  facet_wrap(~ animal) +
  labs(
    x = NULL,
    y = "Fraction of kidney viral abundance in\nhost miRNA seed–targeted barcodes",
    title = "B. Host miRNA seed targeting weighted by kidney barcode levels\n(Observed vs null distribution)"
  ) +
  scale_def_colors +
  scale_def_colors_lines +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "none"
  )







###
# ============================================================
# Seed7 (Top-50 kidney-abundant host miRNAs) observed vs null
# - Unweighted + abundance-weighted (mouse-by-mouse)
# - Length + GC–matched resampling null (N_REP simulations)
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(patchwork)

set.seed(1)

# -----------------------------
# 0) Hit set: seed7 hits restricted to top-50 kidney-abundant host miRNAs
#    (expects you already created: hits_mmu_seed7_kidney with column barcode + hit_mmu_seed7_kidney)
# -----------------------------
# If not already defined earlier in your session, uncomment and use:
# hits_mmu_seed7_kidney <- seed7_mmu_kidney %>%
#   distinct(original_barcode_seq) %>%
#   rename(barcode = original_barcode_seq) %>%
#   mutate(hit_mmu_seed7_kidney = TRUE)

# -----------------------------
# 1) Kidney abundance per mouse x barcode (bc_level)
# -----------------------------
kidney_abundance <- tissue_bc_levels %>%
  filter(is_kidney(organ), bc_level > 0) %>%
  group_by(animal, barcode) %>%
  summarise(kidney_bc_level = sum(bc_level, na.rm = TRUE), .groups = "drop")

# -----------------------------
# 2) Master table: one row per mouse x barcode in kidney_classified
#    Includes status, kidney_bc_level, hit flag, length + GC
# -----------------------------
barcode_master_seed7_kidneyMiR <- kidney_classified %>%
  left_join(kidney_abundance, by = c("animal", "barcode")) %>%
  left_join(hits_mmu_seed7_kidney, by = "barcode") %>%
  mutate(
    kidney_bc_level = if_else(is.na(kidney_bc_level), 0, kidney_bc_level),
    hit_mmu_kidneyMiR = if_else(is.na(hit_mmu_seed7_kidney), FALSE, hit_mmu_seed7_kidney),
    len = nchar(barcode),
    gc  = (str_count(barcode, "[Gg]") + str_count(barcode, "[Cc]")) / pmax(len, 1)
  )

# -----------------------------
# 3) Observed summaries (per mouse, per status)
# -----------------------------
obs_unweighted_kidneyMiR <- barcode_master_seed7_kidneyMiR %>%
  group_by(animal, status) %>%
  summarise(
    pct_hit = mean(hit_mmu_kidneyMiR),
    .groups = "drop"
  )

obs_weighted_kidneyMiR <- barcode_master_seed7_kidneyMiR %>%
  group_by(animal, status) %>%
  summarise(
    pct_weighted = if_else(sum(kidney_bc_level) > 0,
                           sum(kidney_bc_level[hit_mmu_kidneyMiR]) / sum(kidney_bc_level),
                           NA_real_),
    .groups = "drop"
  )

# -----------------------------
# 4) Helper: length+GC binning for matched resampling
#    - length bins: exact length
#    - GC bins: 5 bins within each mouse (quantile-based)
# -----------------------------
master_for_null_kidneyMiR <- barcode_master_seed7_kidneyMiR %>%
  group_by(animal) %>%
  mutate(
    gc_bin  = ntile(gc, 5),
    len_bin = len
  ) %>%
  ungroup()

# Stratified sample: for a given mouse & status target,
# sample from same mouse pool matching counts per (len_bin, gc_bin).
stratified_sample <- function(pool_df, target_df, strata = c("len_bin", "gc_bin")) {
  target_counts <- target_df %>%
    count(across(all_of(strata)), name = "n_target")
  
  sampled <- target_counts %>%
    group_split(across(all_of(strata))) %>%
    map_dfr(function(one_stratum) {
      n_take <- one_stratum$n_target[1]
      st_vals <- one_stratum %>% select(all_of(strata)) %>% slice(1)
      
      candidates <- pool_df
      for (s in strata) {
        candidates <- candidates %>% filter(.data[[s]] == st_vals[[s]][1])
      }
      
      if (nrow(candidates) == 0) {
        pool_df %>% slice_sample(n = n_take, replace = TRUE)
      } else if (nrow(candidates) < n_take) {
        candidates %>% slice_sample(n = n_take, replace = TRUE)
      } else {
        candidates %>% slice_sample(n = n_take, replace = FALSE)
      }
    })
  
  sampled
}

# -----------------------------
# 5) Null distributions (length+GC-matched)
#    Build null for BOTH unweighted and weighted metrics
# -----------------------------
N_REP <- 2000  # increase for smoother violins

make_mmu_null_kidneyMiR <- function(animal_id, status_level) {
  pool   <- master_for_null_kidneyMiR %>% filter(animal == animal_id)
  target <- pool %>% filter(status == status_level)
  if (nrow(target) == 0) return(NULL)
  
  map_dfr(seq_len(N_REP), function(i) {
    samp <- stratified_sample(pool, target, strata = c("len_bin", "gc_bin"))
    tibble(
      animal = animal_id,
      status = status_level,
      rep = i,
      pct_hit = mean(samp$hit_mmu_kidneyMiR),
      pct_weighted = if_else(sum(samp$kidney_bc_level) > 0,
                             sum(samp$kidney_bc_level[samp$hit_mmu_kidneyMiR]) / sum(samp$kidney_bc_level),
                             NA_real_)
    )
  })
}

animals  <- sort(unique(master_for_null_kidneyMiR$animal))
statuses <- levels(master_for_null_kidneyMiR$status)

mmu_null_kidneyMiR <- map_dfr(animals, function(a) {
  map_dfr(statuses, function(s) make_mmu_null_kidneyMiR(a, s))
})

# -----------------------------
# 6) Plots: observed points/lines over null violins
# -----------------------------
pKidneyMiR_C1 <- ggplot() +
  geom_violin(
    data = mmu_null_kidneyMiR,
    aes(x = status, y = pct_hit),
    trim = TRUE
  ) +
  geom_point(
    data = obs_unweighted_kidneyMiR,
    aes(x = status, y = pct_hit, group = animal),
    size = 2
  ) +
  geom_line(
    data = obs_unweighted_kidneyMiR,
    aes(x = status, y = pct_hit, group = animal),
    alpha = 0.6
  ) +
  facet_wrap(~ animal) +
  labs(
    x = NULL,
    y = "Fraction of kidney barcodes\nwith seed match to top kidney-abundant host miRNAs",
    title = "A. Host miRNA seed targeting (top kidney-abundant miRNAs)\n  (Observed vs null expectation)"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

pKidneyMiR_C2 <- ggplot() +
  geom_violin(
    data = mmu_null_kidneyMiR,
    aes(x = status, y = pct_weighted),
    trim = TRUE
  ) +
  geom_point(
    data = obs_weighted_kidneyMiR,
    aes(x = status, y = pct_weighted, group = animal),
    size = 2
  ) +
  geom_line(
    data = obs_weighted_kidneyMiR,
    aes(x = status, y = pct_weighted, group = animal),
    alpha = 0.6
  ) +
  facet_wrap(~ animal) +
  labs(
    x = NULL,
    y = "Fraction of kidney viral abundance in\nbarcodes seed-matched to top kidney-abundant host miRNAs",
    title = "B. Host miRNA seed targeting weighted by kidney barcode levels\n  (Observed vs null expectation)"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

# View
print(pKidneyMiR_C1)
print(pKidneyMiR_C2)

# Single page (2 rows)
pKidneyMiR_C1 / pKidneyMiR_C2

# Save PDF
ggsave(
  filename = "../plots/urine_tissue_manuscript5/mmu_seed7_top50KidneyMiRs_observed_vs_null.pdf",
  plot = pKidneyMiR_C1 / pKidneyMiR_C2,
  width = 6.3,
  height = 8,
  units = "in"
)





















#####top100


# ============================================================
# Seed7 – TOP 100 kidney-abundant host miRNAs
# Observed vs length+GC-matched null distributions
# Mouse-by-mouse analysis (Late-shed vs Late non-shed)
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(patchwork)

set.seed(1)

# -----------------------------
# 1) Load kidney miRNA atlas + define TOP_N
# -----------------------------

TOP_N <- 100

df <- read_xlsx("/stor/work/Sullivan/anik/barcode_project/data/mirna_target_kidney/analysis.rnahybrid/analysis_7feb26/PNAS2020_miRNA_atlas/pnas.2002277117.sd04.xlsx")

kidney_ranked <- df %>%
  filter(Tissue == "Kidney") %>%
  select(miRNA, `Mean scaled expression *`) %>%
  arrange(desc(`Mean scaled expression *`))

top100_miRs <- kidney_ranked$miRNA[1:TOP_N]

# -----------------------------
# 2) Filter seed7_mmu to TOP 100 kidney-abundant miRNAs
# -----------------------------

seed7_mmu_top100 <- seed7_mmu %>%
  filter(matched_miRNA %in% top100_miRs)

hits_mmu_seed7_top100 <- seed7_mmu_top100 %>%
  distinct(original_barcode_seq) %>%
  rename(barcode = original_barcode_seq) %>%
  mutate(hit_mmu_seed7_top100 = TRUE)

# -----------------------------
# 3) Kidney abundance per mouse x barcode
# -----------------------------

kidney_abundance <- tissue_bc_levels %>%
  filter(is_kidney(organ), bc_level > 0) %>%
  group_by(animal, barcode) %>%
  summarise(kidney_bc_level = sum(bc_level, na.rm = TRUE), .groups = "drop")

# -----------------------------
# 4) Master table
# -----------------------------

barcode_master_top100 <- kidney_classified %>%
  left_join(kidney_abundance, by = c("animal", "barcode")) %>%
  left_join(hits_mmu_seed7_top100, by = "barcode") %>%
  mutate(
    kidney_bc_level = if_else(is.na(kidney_bc_level), 0, kidney_bc_level),
    hit_top100 = if_else(is.na(hit_mmu_seed7_top100), FALSE, hit_mmu_seed7_top100),
    len = nchar(barcode),
    gc  = (str_count(barcode, "[Gg]") + str_count(barcode, "[Cc]")) / pmax(len, 1)
  )

# -----------------------------
# 5) Observed summaries
# -----------------------------

obs_unweighted_top100 <- barcode_master_top100 %>%
  group_by(animal, status) %>%
  summarise(
    pct_hit = mean(hit_top100),
    .groups = "drop"
  )

obs_weighted_top100 <- barcode_master_top100 %>%
  group_by(animal, status) %>%
  summarise(
    pct_weighted = if_else(sum(kidney_bc_level) > 0,
                           sum(kidney_bc_level[hit_top100]) / sum(kidney_bc_level),
                           NA_real_),
    .groups = "drop"
  )

# -----------------------------
# 6) Length + GC matched null distributions
# -----------------------------

master_for_null_top100 <- barcode_master_top100 %>%
  group_by(animal) %>%
  mutate(
    gc_bin = ntile(gc, 5),
    len_bin = len
  ) %>%
  ungroup()

stratified_sample <- function(pool_df, target_df, strata = c("len_bin", "gc_bin")) {
  target_counts <- target_df %>%
    count(across(all_of(strata)), name = "n_target")
  
  target_counts %>%
    group_split(across(all_of(strata))) %>%
    map_dfr(function(one_stratum) {
      n_take <- one_stratum$n_target[1]
      st_vals <- one_stratum %>% select(all_of(strata)) %>% slice(1)
      
      candidates <- pool_df
      for (s in strata) {
        candidates <- candidates %>% filter(.data[[s]] == st_vals[[s]][1])
      }
      
      if (nrow(candidates) == 0) {
        pool_df %>% slice_sample(n = n_take, replace = TRUE)
      } else if (nrow(candidates) < n_take) {
        candidates %>% slice_sample(n = n_take, replace = TRUE)
      } else {
        candidates %>% slice_sample(n = n_take, replace = FALSE)
      }
    })
}

N_REP <- 2000

make_null_top100 <- function(animal_id, status_level) {
  pool <- master_for_null_top100 %>% filter(animal == animal_id)
  target <- pool %>% filter(status == status_level)
  if (nrow(target) == 0) return(NULL)
  
  map_dfr(seq_len(N_REP), function(i) {
    samp <- stratified_sample(pool, target)
    tibble(
      animal = animal_id,
      status = status_level,
      rep = i,
      pct_hit = mean(samp$hit_top100),
      pct_weighted = if_else(sum(samp$kidney_bc_level) > 0,
                             sum(samp$kidney_bc_level[samp$hit_top100]) /
                               sum(samp$kidney_bc_level),
                             NA_real_)
    )
  })
}

animals <- sort(unique(master_for_null_top100$animal))
statuses <- levels(master_for_null_top100$status)

null_top100 <- map_dfr(animals, function(a) {
  map_dfr(statuses, function(s) make_null_top100(a, s))
})

# -----------------------------
# 7) Plots
# -----------------------------

p1_top100 <- ggplot() +
  geom_violin(data = null_top100,
              aes(x = status, y = pct_hit),
              trim = TRUE) +
  geom_point(data = obs_unweighted_top100,
             aes(x = status, y = pct_hit, group = animal),
             size = 2) +
  geom_line(data = obs_unweighted_top100,
            aes(x = status, y = pct_hit, group = animal),
            alpha = 0.6) +
  facet_wrap(~ animal) +
  labs(
    x = NULL,
    y = "Fraction of kidney barcodes with seed match\nto top kidney-abundant host miRNAs",
    title = "A. Host miRNA seed targeting (top 100 kidney-abundant miRNAs)\n  (Observed vs null expectation)"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

p2_top100 <- ggplot() +
  geom_violin(data = null_top100,
              aes(x = status, y = pct_weighted),
              trim = TRUE) +
  geom_point(data = obs_weighted_top100,
             aes(x = status, y = pct_weighted, group = animal),
             size = 2) +
  geom_line(data = obs_weighted_top100,
            aes(x = status, y = pct_weighted, group = animal),
            alpha = 0.6) +
  facet_wrap(~ animal) +
  labs(
    x = NULL,
    y = "Fraction of kidney viral abundance in barcodes\nseed-matched to top kidney-abundant host miRNAs",
    title = "B. Host miRNA seed targeting weighted by kidney barcode levels\n  (Observed vs null expectation)"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))




scale_def_colors <- scale_fill_manual(values = c(
  "Late-shed" = "#F8766D",
  "Late non-shed" = "#00BFC4"
))

scale_def_colors_line <- scale_color_manual(values = c(
  "Late-shed" = "#F8766D",
  "Late non-shed" = "#00BFC4"
))

# -----------------------------
# Panel A (Unweighted)
# -----------------------------
p1_top100 <- ggplot() +
  geom_violin(
    data = null_top100,
    aes(x = status, y = pct_hit, fill = status),
    trim = TRUE,
    alpha = 0.6,
    color = NA
  ) +
  geom_point(
    data = obs_unweighted_top100,
    aes(x = status, y = pct_hit, group = animal, color = status),
    size = 2
  ) +
  geom_line(
    data = obs_unweighted_top100,
    aes(x = status, y = pct_hit, group = animal, color = status),
    alpha = 0.6
  ) +
  facet_wrap(~ animal) +
  scale_def_colors +
  scale_def_colors_line +
  labs(
    x = NULL,
    y = "Fraction of kidney barcodes with seed match\nto top kidney-abundant host miRNAs",
    title = "A. Host miRNA seed targeting (top 100 kidney-abundant miRNAs)\n  (Observed vs null expectation)"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "none")


# -----------------------------
# Panel B (Weighted)
# -----------------------------
p2_top100 <- ggplot() +
  geom_violin(
    data = null_top100,
    aes(x = status, y = pct_weighted, fill = status),
    trim = TRUE,
    alpha = 0.6,
    color = NA
  ) +
  geom_point(
    data = obs_weighted_top100,
    aes(x = status, y = pct_weighted, group = animal, color = status),
    size = 2
  ) +
  geom_line(
    data = obs_weighted_top100,
    aes(x = status, y = pct_weighted, group = animal, color = status),
    alpha = 0.6
  ) +
  facet_wrap(~ animal) +
  scale_def_colors +
  scale_def_colors_line +
  labs(
    x = NULL,
    y = "Fraction of kidney viral abundance in barcodes\nseed-matched to top kidney-abundant host miRNAs",
    title = "B. Host miRNA seed targeting weighted by kidney barcode levels\n  (Observed vs null expectation)"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "none")


combined_plot_top100 <- p1_top100 / p2_top100

print(combined_plot_top100)

# -----------------------------
# 8) Save to separate folder
# -----------------------------

#dir.create("../plots/urine_tissue_manuscript5/top100_analysis", showWarnings = FALSE)

ggsave(
  filename = "../plots/urine_tissue_manuscript5/top100_analysis/mmu_seed7_top100_observed_vs_null2.pdf",
  plot = combined_plot_top100,
  width = 6.3,
  height = 8,
  units = "in"
)


# ============================================================
# Export TOP100 kidney-miRNA seed7 summary tables to 1 Excel file (2 sheets)
# ============================================================

library(openxlsx)

out_xlsx <- "../plots/urine_tissue_manuscript5/top100_analysis/mmu_seed7_top100_summary_tables.xlsx"

wb <- createWorkbook()

addWorksheet(wb, "Unweighted_presence")
writeDataTable(wb, "Unweighted_presence", presence_by_mouse_mmu_seed7_kidney)

addWorksheet(wb, "Weighted_abundance")
writeDataTable(wb, "Weighted_abundance", weighted_by_mouse_mmu_seed7_kidney)

saveWorkbook(wb, out_xlsx, overwrite = TRUE)

out_xlsx




