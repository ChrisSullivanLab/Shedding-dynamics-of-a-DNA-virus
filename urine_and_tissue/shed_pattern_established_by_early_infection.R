


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


animal_colors = c('FL' = 'palevioletred3', 'FR' = 'palevioletred1', 'ML' = 'skyblue3', 'MR' = 'skyblue1')
rgb_codes <- col2rgb(animal_colors)

# palevioletred3: RGB(205, 98, 145)
# palevioletred1: RGB(255, 130, 171)
# skyblue3: RGB(80, 170, 255)
# skyblue1: RGB(135, 206, 255)

##### bar plot of present and absent in kidney:
urine_bc_levels <- read_tsv("urine_bc_levels.tsv", show_col_types = FALSE)
tissue_bc_levels <- read_tsv("tissue_bc_levels.tsv", show_col_types = FALSE)
cutoff_99pct_stock_barcodes <- read_tsv("cutoff_99pct_stock_barcodes.tsv", show_col_types = FALSE)

SAVE_DIR <- "../plots/urine_tissue_manuscript4"

#Kidney barcodes that ARE vs ARE NOT in last-K urine


library(dplyr)
library(ggplot2)
library(tidyr)

LAST_K <- 2
#LAST_K <- if (exists("LAST_K")) LAST_K else 2

if (!exists("is_kidney")) {
  is_kidney <- function(x) grepl("kidney", x, ignore.case = TRUE)
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
    present_lastK = any(bc_level > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(present_lastK) %>%
  select(animal, barcode)

kidney_barcodes <- tissue_bc_levels %>%
  filter(is_kidney(organ)) %>%
  group_by(animal, barcode) %>%
  summarise(kidney_present = any(bc_level > 0, na.rm = TRUE), .groups = "drop") %>%
  filter(kidney_present) %>%
  select(animal, barcode)

kidney_classified <- kidney_barcodes %>%
  left_join(present_lastK %>% mutate(in_lastK = TRUE),
            by = c("animal","barcode")) %>%
  mutate(
    in_lastK = if_else(is.na(in_lastK), FALSE, in_lastK),
    status   = if_else(in_lastK,
                       "Present in last-K urine",
                       "Absent in last-K urine")
  )

kidney_counts <- kidney_classified %>%
  count(animal, status, name = "n_barcodes")

kidney_counts <- kidney_counts %>%
  mutate(status = factor(status,
                         levels = c("Present in last-K urine",
                                    "Absent in last-K urine")))

p_kidney_bar_lastK <- ggplot(kidney_counts,
                             aes(x = animal, y = n_barcodes, fill = status)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  labs(
    title = paste0("Kidney barcodes by last-", LAST_K, " urine presence"),
    x = "Mouse",
    y = "Number of kidney barcodes",
    fill = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )


print(p_kidney_bar_lastK)

ggsave(
  filename = paste0("fig1A_kidney_barcode_counts_present_vs_absent_lastK.pdf"),
  plot     = p_kidney_bar_lastK,
  path     = if (exists("SAVE_DIR")) SAVE_DIR else ".",
  width    = 6, height = 4.5
)




#########15oct___top10 stuff over time






TOP_N  <- 50   
LAST_K <- 2    
REQUIRE_ALL_LASTK <- FALSE  # TRUE = present in *all* of the last K; FALSE = present in *any* of the last K

is_kidney <- function(x) grepl("kidney", x, ignore.case = TRUE)


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



animals <- names(animal_colors)
plot_diversity_donut_w_colored_barcodes <- function(bc_table,
                                                    barcode_color_highlights,
                                                    faceting_variable) {
  
  valid_highlights <- barcode_color_highlights[barcode_color_highlights %in% bc_table$barcode]
  
  if (length(valid_highlights) == 0) {
    stop("None of the barcode highlights are present in the dataset.")
  }
  
  winner_and_bg_colors <- set_names(c(scales::hue_pal()(length(valid_highlights)), 'gray75', 'gray60'),
                                    nm = c(valid_highlights, 'background1', 'background2'))
  
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

donut_plot_list <- lapply(animals, function(a) {
  bc_table <- urine_bc_levels %>% filter(animal == a)  # full table so grey background includes all
  highlights <- winner_barcodes_by_animal[[a]]
  if (is.null(highlights) || length(highlights) == 0) {
    # fall back to empty panel if no qualifying barcodes
    return(ggplot() + theme_void() + labs(title = a))
  }
  plot_diversity_donut_w_colored_barcodes(
    bc_table = bc_table,
    barcode_color_highlights = highlights,
    faceting_variable = Days_pi
  ) + labs(title = a)
})

donut_combined <- gridExtra::grid.arrange(grobs = donut_plot_list, ncol = 2)
ggsave(
  filename = "fig2A_donut_topKidney_lastK.pdf",
  plot = donut_combined,
  path = "../plots/urine_tissue_manuscript4",
  width = 8, height = 8
)
donut_combined



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

pdf("../plots/urine_tissue_manuscript4/fig3A_sum_overlay_topKidney_lastK.pdf", onefile = TRUE)
grid.arrange(
  sum_overlay_plots[["FL"]],
  sum_overlay_plots[["FR"]],
  sum_overlay_plots[["ML"]],
  sum_overlay_plots[["MR"]],
  nrow = 4, ncol = 1
)
grid.text("Days Post-Injection", x = unit(.48, "npc"), y = unit(0.04, "npc"),
          just = c("left","top"), gp = gpar(fontsize = 11))
grid.text("Barcode level", x = unit(.005, "npc"), y = unit(0.5, "npc"),
          just = c("left","top"), gp = gpar(fontsize = 11), rot = 90)
dev.off()

grid.arrange(
  sum_overlay_plots[["FL"]],
  sum_overlay_plots[["FR"]],
  sum_overlay_plots[["ML"]],
  sum_overlay_plots[["MR"]],
  nrow = 4, ncol = 1
)
grid.text("Days Post-Injection", x = unit(.48, "npc"), y = unit(0.04, "npc"),
          just = c("left","top"), gp = gpar(fontsize = 11))
grid.text("Barcode level", x = unit(.005, "npc"), y = unit(0.5, "npc"),
          just = c("left","top"), gp = gpar(fontsize = 11), rot = 90)


##rank in stock

stock_rank_tbl <- cutoff_99pct_stock_barcodes %>%
  mutate(stock_rank = min_rank(desc(count))) %>%
  select(barcode, count, stock_rank)

kidney_present_with_stock <- kidney_rank %>%   # from the “present last-K” section
  left_join(stock_rank_tbl, by = "barcode") %>%
  mutate(kidney_rank_f = as_factor(kidney_rank))

p_kidney_present_in_stock <- kidney_present_with_stock %>%
  ggplot(aes(x = fct_relevel(kidney_rank_f, rev), y = stock_rank)) +
  geom_point(aes(color = as_factor(animal)), size = 3, show.legend = FALSE) +
  facet_wrap(vars(animal)) +
  labs(title = "",
       x = "Rank of Top Kidney (Present in last-K urine) Barcode",
       y = "Rank in Stock") +
  scale_y_continuous(limits = c(1, nrow(cutoff_99pct_stock_barcodes))) +
  scale_x_discrete(breaks = seq(1, TOP_N)) +
  scale_color_manual(values = animal_colors) +
  coord_flip() +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(linetype = "dotted", color = "gray70"),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(2, "lines"),
        strip.background = element_rect(fill = "gray90", color = "gray90"))

print(p_kidney_present_in_stock)

ggsave("fig5A_rank_kidney_present_in_stock_dot.pdf",
       plot = p_kidney_present_in_stock,
       path = "../plots/urine_tissue_manuscript4",
       width = 8, height = 8)

# median rank
kidney_present_median_stock_rank <- kidney_present_with_stock %>%
  group_by(animal) %>%
  summarise(median_stock_rank = median(stock_rank, na.rm = TRUE), .groups = "drop")

kidney_present_median_stock_rank







#### Now not present


TOP_N  <- 50   
LAST_K <- 2    
REQUIRE_ALL_LASTK <- FALSE  


is_kidney <- function(x) grepl("kidney|renal", x, ignore.case = TRUE)


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


animals <- names(animal_colors)

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

donut_combined_notLastK <- gridExtra::grid.arrange(grobs = donut_plot_list_notLastK, ncol = 2)
ggsave(
  filename = "fig2B_donut_topKidney_NOTlastK.pdf",
  plot = donut_combined_notLastK,
  path = "../plots/urine_tissue_manuscript4",
  width = 8, height = 8
)
donut_combined_notLastK




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

pdf("../plots/urine_tissue_manuscript4/fig3B_sum_overlay_topKidney_NOTlastK.pdf", onefile = TRUE)
grid.arrange(
  sum_overlay_plots_notLastK[["FL"]],
  sum_overlay_plots_notLastK[["FR"]],
  sum_overlay_plots_notLastK[["ML"]],
  sum_overlay_plots_notLastK[["MR"]],
  nrow = 4, ncol = 1
)
grid.text("Days Post-Injection", x = unit(.48, "npc"), y = unit(0.04, "npc"),
          just = c("left","top"), gp = gpar(fontsize = 11))
grid.text("Barcode level", x = unit(.005, "npc"), y = unit(0.5, "npc"),
          just = c("left","top"), gp = gpar(fontsize = 11), rot = 90)
dev.off()

grid.arrange(
  sum_overlay_plots_notLastK[["FL"]],
  sum_overlay_plots_notLastK[["FR"]],
  sum_overlay_plots_notLastK[["ML"]],
  sum_overlay_plots_notLastK[["MR"]],
  nrow = 4, ncol = 1
)
grid.text("Days Post-Injection", x = unit(.48, "npc"), y = unit(0.04, "npc"),
          just = c("left","top"), gp = gpar(fontsize = 11))
grid.text("Barcode level", x = unit(.005, "npc"), y = unit(0.5, "npc"),
          just = c("left","top"), gp = gpar(fontsize = 11), rot = 90)





#rank in stock - not present group

stock_rank_tbl <- cutoff_99pct_stock_barcodes %>%
  mutate(stock_rank = min_rank(desc(count))) %>%
  select(barcode, count, stock_rank)

kidney_notLastK_with_stock <- kidney_rank_notLastK %>%     # from the NOT-last-K selection
  left_join(stock_rank_tbl, by = "barcode") %>%
  mutate(kidney_rank_f = as_factor(kidney_rank))

p_kidney_notLastK_in_stock <- kidney_notLastK_with_stock %>%
  ggplot(aes(x = fct_relevel(kidney_rank_f, rev), y = stock_rank)) +
  geom_point(aes(color = as_factor(animal)), size = 3, show.legend = FALSE) +
  facet_wrap(vars(animal)) +
  labs(title = "",
       x = "Rank of Top Kidney (NOT last-K) Barcode",
       y = "Rank in Stock") +
  scale_y_continuous(limits = c(1, nrow(cutoff_99pct_stock_barcodes))) +
  scale_x_discrete(breaks = seq(1, TOP_N)) +
  scale_color_manual(values = animal_colors) +
  coord_flip() +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(linetype = "dotted", color = "gray70"),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(2, "lines"),
        strip.background = element_rect(fill = "gray90", color = "gray90"))

print(p_kidney_notLastK_in_stock)

ggsave("fig5B_rank_kidney_NOTlastK_in_stock_dot.pdf",
       plot = p_kidney_notLastK_in_stock,
       path = "../plots/urine_tissue_manuscript3",
       width = 8, height = 8)


kidney_notLastK_median_stock_rank <- kidney_notLastK_with_stock %>%
  group_by(animal) %>%
  summarise(median_stock_rank = median(stock_rank, na.rm = TRUE), .groups = "drop")

kidney_notLastK_median_stock_rank









###### customizable top how many


TOP_N         <- 50                     # or 10
TOP_LABEL     <- paste0("Top-", TOP_N)
CAP_LIST      <- c(10)   # order matters: highest -> lowest
SAVE_DIR      <- "../plots/urine_tissue_manuscript3"



is_brain <- function(x) grepl("brain|cerebr|cortex|hippo", x, ignore.case = TRUE)
is_lung  <- function(x) grepl("lung|pulmo",               x, ignore.case = TRUE)

build_topN_by_organ <- function(tissue_tbl, organ_predicate, load_name, rank_name, top_n = TOP_N) {
  ranked <- tissue_tbl %>%
    dplyr::filter(organ_predicate(organ)) %>%
    dplyr::group_by(animal, barcode) %>%
    dplyr::summarise(!!load_name := sum(bc_level, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(animal) %>%
    dplyr::arrange(dplyr::desc(.data[[load_name]]), .by_group = TRUE) %>%
    dplyr::mutate(!!rank_name := dplyr::row_number()) %>%
    dplyr::ungroup()
  
  ranked %>%
    dplyr::group_by(animal) %>%
    dplyr::arrange(.data[[rank_name]], .by_group = TRUE) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup()
}

random_topN_control <- function(urine_tbl, selected_tbl, n = TOP_N, seed = 1) {
  set.seed(seed)
  avail <- urine_tbl %>% dplyr::distinct(animal, barcode)
  pool  <- dplyr::anti_join(avail, selected_tbl %>% dplyr::distinct(animal, barcode),
                            by = c("animal","barcode"))
  pool %>%
    dplyr::group_by(animal) %>%
    dplyr::mutate(.rand = runif(dplyr::n())) %>%
    dplyr::arrange(.rand, .by_group = TRUE) %>%
    dplyr::mutate(.rank = dplyr::row_number()) %>%
    dplyr::filter(.rank <= n) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.rand, -.rank)
}

plot_barcode_lines_capped_for_mouse <- function(urine_tbl, selected_tbl, mouse_id,
                                                title_suffix = "", cap = 100) {
  sel_bcs <- selected_tbl %>%
    dplyr::filter(animal == mouse_id) %>%
    dplyr::pull(barcode) %>% unique()
  
  if (length(sel_bcs) == 0) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::labs(title = paste(mouse_id, "(no data)")))
  }
  
  bc_ts <- urine_tbl %>%
    dplyr::semi_join(tibble::tibble(animal = mouse_id, barcode = sel_bcs),
                     by = c("animal","barcode")) %>%
    dplyr::group_by(animal, barcode, Days_pi) %>%
    dplyr::summarise(bc_day_level = sum(bc_level, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(bc_capped = pmin(bc_day_level, cap))
  
  pal <- setNames(scales::hue_pal()(length(sel_bcs)), sel_bcs)
  
  ggplot2::ggplot(bc_ts, ggplot2::aes(x = Days_pi, y = bc_capped,
                                      group = barcode, color = barcode)) +
    ggplot2::geom_line(linewidth = 0.35, alpha = 0.85, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = pal) +
    ggplot2::scale_y_continuous(limits = c(0, cap),
                                name = paste0("Barcode level (capped @ ", cap, ")")) +
    ggplot2::labs(title = title_suffix, x = "Days post-injection") +
    ggplot2::theme_minimal()
}


safe_arrange <- function(plot_list, nrow = NULL, ncol = NULL) {
  ord <- c("FL","FR","ML","MR")
  have <- intersect(ord, names(plot_list))
  do.call(gridExtra::grid.arrange, c(plot_list[have], nrow = nrow, ncol = ncol))
}

make_overlays_across_caps <- function(selection_tbl, selection_label, file_stub) {
  animals <- names(animal_colors)
  for (cap in CAP_LIST) {
    grobs <- setNames(vector("list", length(animals)), animals)
    for (a in animals) {
      grobs[[a]] <- plot_barcode_lines_capped_for_mouse(
        urine_tbl    = urine_bc_levels,
        selected_tbl = selection_tbl,
        mouse_id     = a,
        title_suffix = paste0(a, " — ", selection_label, " (", TOP_LABEL, "), cap=", cap),
        cap          = cap
      )
    }
    out_file <- file.path(SAVE_DIR, paste0(file_stub, "_", TOP_LABEL, "_cap", cap, ".pdf"))
    pdf(out_file, width = 12, height = 10, onefile = TRUE)
    safe_arrange(grobs, nrow = 2, ncol = 2)
    dev.off()
  }
}



selected_kidney_present_tbl  <- if (exists("kidney_rank_topN"))           kidney_rank_topN           else kidney_rank
selected_kidney_notLastK_tbl <- if (exists("kidney_rank_notLastK_topN"))  kidney_rank_notLastK_topN  else kidney_rank_notLastK


if (!exists("brain_topN")) {
  brain_topN <- build_topN_by_organ(
    tissue_tbl = tissue_bc_levels,
    organ_predicate = is_brain,
    load_name = "brain_load",
    rank_name = "brain_rank",
    top_n = TOP_N
  )
}
if (!exists("lung_topN")) {
  lung_topN <- build_topN_by_organ(
    tissue_tbl = tissue_bc_levels,
    organ_predicate = is_lung,
    load_name = "lung_load",
    rank_name = "lung_rank",
    top_n = TOP_N
  )
}

# Random controls (per mouse) for Kidney PRESENT, Brain, Lung
rand_kidney_present_topN <- random_topN_control(urine_bc_levels, selected_kidney_present_tbl, n = TOP_N, seed = 401)
rand_brain_topN          <- random_topN_control(urine_bc_levels, brain_topN,                   n = TOP_N, seed = 402)
rand_lung_topN           <- random_topN_control(urine_bc_levels, lung_topN,                    n = TOP_N, seed = 403)



# Kidney Top-N PRESENT
make_overlays_across_caps(
  selection_tbl   = selected_kidney_present_tbl,
  selection_label = "Kidney PRESENT",
  file_stub       = "overlay_BARCODE_ONLY_Kidney_PRESENT"
)

# Kidney Top-N NOT-lastK
make_overlays_across_caps(
  selection_tbl   = selected_kidney_notLastK_tbl,
  selection_label = "Kidney NOT-lastK",
  file_stub       = "overlay_BARCODE_ONLY_Kidney_NOTlastK"
)

# Brain Top-N (no last-K)
make_overlays_across_caps(
  selection_tbl   = brain_topN,
  selection_label = "Brain",
  file_stub       = "overlay_BARCODE_ONLY_Brain"
)

# Lung Top-N (no last-K)
make_overlays_across_caps(
  selection_tbl   = lung_topN,
  selection_label = "Lung",
  file_stub       = "overlay_BARCODE_ONLY_Lung"
)

# Random controls (Kidney PRESENT / Brain / Lung)
make_overlays_across_caps(
  selection_tbl   = rand_kidney_present_topN,
  selection_label = paste0("Random (Kidney PRESENT baseline)"),
  file_stub       = "overlay_BARCODE_ONLY_Random_fromKidneyPresent"
)

make_overlays_across_caps(
  selection_tbl   = rand_brain_topN,
  selection_label = paste0("Random (Brain baseline)"),
  file_stub       = "overlay_BARCODE_ONLY_Random_fromBrain"
)

make_overlays_across_caps(
  selection_tbl   = rand_lung_topN,
  selection_label = paste0("Random (Lung baseline)"),
  file_stub       = "overlay_BARCODE_ONLY_Random_fromLung"
)


preview_tbl   <- selected_kidney_present_tbl   # 👈 change this to kidney_notLastK_tbl, brain_topN, lung_topN, etc.
preview_label <- "Kidney PRESENT"

for (cap in CAP_LIST) {
  grobs <- setNames(vector("list", length(names(animal_colors))), names(animal_colors))
  for (a in names(animal_colors)) {
    grobs[[a]] <- plot_barcode_lines_capped_for_mouse(
      urine_tbl    = urine_bc_levels,
      selected_tbl = preview_tbl,
      mouse_id     = a,
      title_suffix = paste0(a, " — ", preview_label, " (", TOP_LABEL, "), cap=", cap),
      cap          = cap
    )
  }
  cat("\n===== Previewing:", preview_label, " cap =", cap, "=====\n")
  print(safe_arrange(grobs, nrow = 2, ncol = 2))
}


combined_pdf <- file.path(SAVE_DIR, paste0("Fig4_ALL_overlay_BARCODE_ONLY_", TOP_LABEL, "_multiCap_combined.pdf"))

pdf(combined_pdf, width = 12, height = 10, onefile = TRUE)

render_overlay_set <- function(selection_tbl, selection_label) {
  for (cap in CAP_LIST) {
    grobs <- setNames(vector("list", length(names(animal_colors))), names(animal_colors))
    for (a in names(animal_colors)) {
      grobs[[a]] <- plot_barcode_lines_capped_for_mouse(
        urine_tbl    = urine_bc_levels,
        selected_tbl = selection_tbl,
        mouse_id     = a,
        title_suffix = paste0(a, " — ", selection_label, " (", TOP_LABEL, "), cap=", cap),
        cap          = cap
      )
    }
    grid::grid.newpage()
    grid::grid.text(
      paste0(selection_label, " (cap=", cap, ")"),
      x = 0.5, y = 0.97, gp = grid::gpar(fontsize = 14, fontface = "bold")
    )
    safe_arrange(grobs, nrow = 2, ncol = 2)
  }
}

render_overlay_set(selected_kidney_present_tbl,  "Kidney PRESENT")
render_overlay_set(selected_kidney_notLastK_tbl, "Kidney NOT-lastK")
render_overlay_set(brain_topN,                   "Brain Top-N")
render_overlay_set(lung_topN,                    "Lung Top-N")
render_overlay_set(rand_kidney_present_topN,     "Random")
#render_overlay_set(rand_brain_topN,              "Random (Brain baseline)")
#render_overlay_set(rand_lung_topN,               "Random (Lung baseline)")

dev.off()

cat("\n✅ Combined multi-page PDF saved to:\n", combined_pdf, "\n")



library(ggplot2)

bc_dist <- urine_bc_levels %>%
  mutate(bc_level = as.numeric(bc_level)) %>%
  filter(!is.na(bc_level) & bc_level > 0)   # remove zeros for log plotting

thr_values <- c(0.2, 1, 5)

ggplot(bc_dist, aes(x = bc_level)) +
  geom_histogram(bins = 200, fill = "gray70", color = "gray40") +
  geom_vline(xintercept = thr_values, linetype = "dashed",
             color = c("darkgreen", "orange", "red"), linewidth = 0.6) +
  scale_x_log10() +
  labs(
    title = "Distribution of bc_level values across all urine samples",
    x = "bc_level (log10 scale)",
    y = "Count of barcode–timepoints",
    subtitle = "Dashed lines: thresholds = 0.2 (green), 1 (orange), 5 (red)"
  ) +
  theme_minimal(base_size = 12)







stock_rank_tbl <- cutoff_99pct_stock_barcodes %>%
  mutate(stock_rank = min_rank(desc(count))) %>%
  select(barcode, count, stock_rank)

kidney_present_with_stock <- kidney_rank %>%   # from the “present last-K” section
  left_join(stock_rank_tbl, by = "barcode") %>%
  mutate(kidney_rank_f = as_factor(kidney_rank))

p_kidney_present_in_stock <- kidney_present_with_stock %>%
  ggplot(aes(x = fct_relevel(kidney_rank_f, rev), y = stock_rank)) +
  geom_point(aes(color = as_factor(animal)), size = 3, show.legend = FALSE) +
  facet_wrap(vars(animal)) +
  labs(title = "",
       x = "Rank of Top Kidney (Present in last-K urine) Barcode",
       y = "Rank in Stock") +
  scale_y_continuous(limits = c(1, nrow(cutoff_99pct_stock_barcodes))) +
  scale_x_discrete(breaks = seq(1, TOP_N)) +
  scale_color_manual(values = animal_colors) +
  coord_flip() +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(linetype = "dotted", color = "gray70"),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(2, "lines"),
        strip.background = element_rect(fill = "gray90", color = "gray90"))

print(p_kidney_present_in_stock)

ggsave("fig5A_rank_kidney_present_in_stock_dot.pdf",
       plot = p_kidney_present_in_stock,
       path = "../plots/urine_tissue_manuscript3",
       width = 8, height = 8)

# (Optional summary: median stock rank per animal)
kidney_present_median_stock_rank <- kidney_present_with_stock %>%
  group_by(animal) %>%
  summarise(median_stock_rank = median(stock_rank, na.rm = TRUE), .groups = "drop")

kidney_present_median_stock_rank




stock_rank_tbl <- cutoff_99pct_stock_barcodes %>%
  mutate(stock_rank = min_rank(desc(count))) %>%
  select(barcode, count, stock_rank)

kidney_notLastK_with_stock <- kidney_rank_notLastK %>%    
  left_join(stock_rank_tbl, by = "barcode") %>%
  mutate(kidney_rank_f = as_factor(kidney_rank))

p_kidney_notLastK_in_stock <- kidney_notLastK_with_stock %>%
  ggplot(aes(x = fct_relevel(kidney_rank_f, rev), y = stock_rank)) +
  geom_point(aes(color = as_factor(animal)), size = 3, show.legend = FALSE) +
  facet_wrap(vars(animal)) +
  labs(title = "",
       x = "Rank of Top Kidney (NOT last-K) Barcode",
       y = "Rank in Stock") +
  scale_y_continuous(limits = c(1, nrow(cutoff_99pct_stock_barcodes))) +
  scale_x_discrete(breaks = seq(1, TOP_N)) +
  scale_color_manual(values = animal_colors) +
  coord_flip() +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(linetype = "dotted", color = "gray70"),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(2, "lines"),
        strip.background = element_rect(fill = "gray90", color = "gray90"))

print(p_kidney_notLastK_in_stock)

ggsave("fig5B_rank_kidney_NOTlastK_in_stock_dot.pdf",
       plot = p_kidney_notLastK_in_stock,
       path = "../plots/urine_tissue_manuscript3",
       width = 8, height = 8)


kidney_notLastK_median_stock_rank <- kidney_notLastK_with_stock %>%
  group_by(animal) %>%
  summarise(median_stock_rank = median(stock_rank, na.rm = TRUE), .groups = "drop")

kidney_notLastK_median_stock_rank




# "History of a single high event" for kidney non-shed vs random/other organs


suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(purrr)
})

# -----------------
# Tunables
# -----------------
SMOLDER_TOP_N      <- if (exists("TOP_N")) TOP_N else 50  
HISTORY_THRESH     <- 1                                
LOG_EPS            <- 1e-6                           
SAVE_DIR_LOCAL     <- if (exists("SAVE_DIR")) SAVE_DIR else "."


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
  mutate(cohort = "Kidney_NOTlastK")

kidney_present_topN <- if (exists("kidney_rank_topN")) kidney_rank_topN else kidney_rank
cohort_kidney_present <- kidney_present_topN %>%
  distinct(animal, barcode) %>%
  mutate(cohort = "Kidney_PresentLastK")

cohort_brain <- if (exists("brain_topN")) {
  brain_topN %>% distinct(animal, barcode) %>% mutate(cohort = "Brain_TopN")
} else {
  tibble(animal = character(), barcode = character(), cohort = character())
}

cohort_lung <- if (exists("lung_topN")) {
  lung_topN %>% distinct(animal, barcode) %>% mutate(cohort = "Lung_TopN")
} else {
  tibble(animal = character(), barcode = character(), cohort = character())
}

# Random-from-urine baseline, matched to kidney_notLastK
urine_universe <- urine_bc_levels %>% distinct(animal, barcode)
random_for_kidney_notLastK <- random_topN_from_urine(
  urine_tbl   = urine_bc_levels,
  exclude_tbl = cohort_kidney_notLastK,
  n           = SMOLDER_TOP_N,
  seed        = 777
) %>%
  mutate(cohort = "Random_fromUrine")

cohorts_all <- bind_rows(
  cohort_kidney_notLastK,
  cohort_kidney_present,
  #cohort_brain,
  #cohort_lung,
  random_for_kidney_notLastK
) %>%
  distinct(animal, barcode, cohort)

#  "lowest single high amount" = max urine bc_level per barcode

max_urine_bc <- urine_bc_levels %>%
  group_by(animal, barcode) %>%
  summarise(max_bc = max(bc_level, na.rm = TRUE), .groups = "drop")

cohorts_with_max <- cohorts_all %>%
  left_join(max_urine_bc, by = c("animal","barcode")) %>%
  mutate(max_bc = replace_na(max_bc, 0))



p_max_bc_violin <- ggplot(
  cohorts_with_max,
  aes(x = cohort, y = pmax(max_bc, LOG_EPS))
) +
  geom_violin(trim = TRUE, scale = "width", fill = "grey80", color = "grey50") +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.4, size = 1) +
  scale_y_log10() +
  facet_wrap(vars(animal), nrow = 2) +
  labs(
    title = "Max urine bc_level per barcode (\"lowest single high amount\")",
    x     = "",
    y     = "Max bc_level across urine samples (log10 scale)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor.y = element_blank()
  )

print(p_max_bc_violin)

ggsave(
  filename = paste0("fig_max_urine_level_by_cohort_TopN", SMOLDER_TOP_N, ".pdf"),
  plot     = p_max_bc_violin,
  path     = SAVE_DIR_LOCAL,
  width    = 10, height = 6
)



history_summary <- cohorts_with_max %>%
  mutate(ever_high = max_bc >= HISTORY_THRESH) %>%
  group_by(animal, cohort) %>%
  summarise(
    n        = n(),
    n_ever   = sum(ever_high),
    pct_ever = 100 * n_ever / pmax(n, 1),
    .groups  = "drop"
  )

p_history_bar <- ggplot(history_summary,
                        aes(x = cohort, y = pct_ever, fill = cohort)) +
  geom_col(width = 0.65) +
  facet_wrap(vars(animal), nrow = 2) +
  labs(
    title = paste0("Fraction of barcodes ever ≥ ", HISTORY_THRESH, " in urine"),
    x     = "",
    y     = "% of barcodes with at least one high-shed event"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

print(p_history_bar)

ggsave(
  filename = paste0("fig_fraction_ever_high_by_cohort_TopN", SMOLDER_TOP_N,
                    "_thr", HISTORY_THRESH, ".pdf"),
  plot     = p_history_bar,
  path     = SAVE_DIR_LOCAL,
  width    = 10, height = 6
)


cat("\n--- Summary: max urine bc_level per cohort ---\n")
print(
  cohorts_with_max %>%
    group_by(cohort) %>%
    summarise(
      median_max_bc = median(max_bc, na.rm = TRUE),
      p75_max_bc    = quantile(max_bc, 0.75, na.rm = TRUE),
      p90_max_bc    = quantile(max_bc, 0.90, na.rm = TRUE),
      n             = n(),
      .groups = "drop"
    )
)

cat("\n--- Summary: fraction ever ≥ ", HISTORY_THRESH, " by cohort & mouse ---\n", sep = "")
print(history_summary)



















# =========================================================
# Top-N smolderers: enrichment in late-high-urine, kidney-present, kidney-notLastK vs random
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
})

# -----------------
# Tunables
# -----------------
SMOLDER_TOP_N <- 50
DETECT        <- .00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001#if (exists("DETECT_THRESH")) DETECT_THRESH else 0.2
LAST_K        <- 2
SAVE_DIR_LOCAL <- if (exists("SAVE_DIR")) SAVE_DIR else "."


if (!exists("compute_pct_detected_per_barcode")) {
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

if (!exists("is_kidney")) {
  is_kidney <- function(x) grepl("kidney|renal", x, ignore.case = TRUE)
}


urine_universe <- urine_bc_levels %>% distinct(animal, barcode)

kidney_present_topN <- if (exists("kidney_rank_topN"))          kidney_rank_topN          else kidney_rank
kidney_notLastK_topN <- if (exists("kidney_rank_notLastK_topN")) kidney_rank_notLastK_topN else kidney_rank_notLastK

cat_set_kidney_present <- kidney_present_topN %>%
  distinct(animal, barcode) %>%
  mutate(category = "Kidney_PresentLastK")

cat_set_kidney_notLastK <- kidney_notLastK_topN %>%
  distinct(animal, barcode) %>%
  mutate(category = "Kidney_NotLastK")

kidney_all_rank <- tissue_bc_levels %>%
  filter(is_kidney(organ)) %>%
  group_by(animal, barcode) %>%
  summarise(kidney_load = sum(bc_level, na.rm = TRUE), .groups = "drop") %>%
  group_by(animal) %>%
  arrange(desc(kidney_load), .by_group = TRUE) %>%
  mutate(kidney_rank_all = row_number()) %>%
  ungroup()

kidney_high_topN <- kidney_all_rank %>%
  group_by(animal) %>%
  arrange(kidney_rank_all, .by_group = TRUE) %>%
  slice_head(n = SMOLDER_TOP_N) %>%
  ungroup()

cat_set_kidney_high <- kidney_high_topN %>%
  distinct(animal, barcode) %>%
  mutate(category = "Kidney_HighLoad")

last_k_days <- urine_bc_levels %>%
  distinct(animal, Days_pi) %>%
  group_by(animal) %>%
  arrange(Days_pi, .by_group = TRUE) %>%
  slice_tail(n = LAST_K) %>%
  mutate(in_lastK = TRUE) %>%
  ungroup()

urine_lastK_sum <- urine_bc_levels %>%
  inner_join(last_k_days, by = c("animal","Days_pi")) %>%
  group_by(animal, barcode) %>%
  summarise(lastK_urine_load = sum(bc_level, na.rm = TRUE), .groups = "drop")

urineLate_highshed_topN <- urine_lastK_sum %>%
  group_by(animal) %>%
  arrange(desc(lastK_urine_load), .by_group = TRUE) %>%
  slice_head(n = SMOLDER_TOP_N) %>%
  ungroup()

cat_set_urine_lateHigh <- urineLate_highshed_topN %>%
  distinct(animal, barcode) %>%
  mutate(category = "Urine_LateHighShed")

categories_all <- bind_rows(
  cat_set_kidney_present,
  cat_set_kidney_notLastK,
  cat_set_kidney_high,
  cat_set_urine_lateHigh
) %>%
  distinct(animal, barcode, category)


smolder_all <- compute_pct_detected_per_barcode(
  urine_tbl    = urine_bc_levels,
  selected_tbl = urine_universe,
  detect_thresh = DETECT
) %>%
  mutate(pct_detected = replace_na(pct_detected, 0))

top_smolderers <- smolder_all %>%
  group_by(animal) %>%
  arrange(desc(pct_detected), .by_group = TRUE) %>%
  slice_head(n = SMOLDER_TOP_N) %>%
  ungroup() %>%
  select(animal, barcode, pct_detected) %>%
  mutate(group = "TopSmolder")

# Random control: same size per mouse, from urine universe
rand_smolder_ctrl <- random_topN_from_urine(
  urine_tbl   = urine_bc_levels,
  exclude_tbl = top_smolderers,
  n           = SMOLDER_TOP_N,
  seed        = 2025
) %>%
  mutate(group = "Random")

smolder_vs_random <- bind_rows(
  top_smolderers %>% select(animal, barcode, group),
  rand_smolder_ctrl %>% select(animal, barcode, group)
)


smolder_cat <- smolder_vs_random %>%
  left_join(categories_all, by = c("animal","barcode")) %>%
  mutate(category = replace_na(category, "None"))

# Summaries per animal, category, and group (TopSmolder vs Random)
enrich_summary <- smolder_cat %>%
  group_by(animal, category, group) %>%
  summarise(
    n_in_cat    = n(),           
    .groups     = "drop_last"
  ) %>%
  group_by(animal, group) %>%
  mutate(
    n_total_group = sum(n_in_cat),
    frac_in_cat   = n_in_cat / pmax(n_total_group, 1)
  ) %>%
  ungroup()

print(enrich_summary)



enrich_tests <- smolder_cat %>%
  # Mark whether barcode is in the given category vs not
  group_by(animal) %>%
  group_modify(~{
    df <- .x
    out_lst <- list()
    cats <- unique(df$category)
    
    for (cat_i in cats) {
      # in_cat vs not_in_cat for this category
      df_cat <- df %>%
        mutate(in_cat = (category == cat_i))
      
      tab <- table(df_cat$group, df_cat$in_cat)  
      
      if (all(dim(tab) == c(2,2))) {
        ft <- fisher.test(tab)
        out_lst[[length(out_lst) + 1]] <- tibble(
          category = cat_i,
          OR       = unname(ft$estimate),
          p_value  = ft$p.value,
          n_top    = sum(df_cat$group == "TopSmolder"),
          n_rand   = sum(df_cat$group == "Random")
        )
      }
    }
    bind_rows(out_lst)
  }) %>%
  ungroup()

print(enrich_tests)



plot_enrich <- enrich_summary %>%
  filter(category != "None") %>%
  mutate(category = factor(category,
                           levels = c("Urine_LateHighShed",
                                      "Kidney_PresentLastK",
                                      "Kidney_NotLastK",
                                      "Kidney_HighLoad")))

p_frac_cat <- ggplot(plot_enrich,
                     aes(x = category, y = frac_in_cat, fill = group)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  facet_wrap(vars(animal), nrow = 2) +
  labs(
    title = paste0("Top-", SMOLDER_TOP_N, " smolderers: category membership vs random"),
    x     = "Category",
    y     = "Fraction of barcodes in category",
    fill  = ""
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

print(p_frac_cat)

ggsave(
  filename = paste0("enrichment_TopSmolder_vs_Random_by_category_Top", SMOLDER_TOP_N, ".pdf"),
  plot     = p_frac_cat,
  path     = SAVE_DIR_LOCAL,
  width    = 10, height = 6
)




#####for grant-writing course



plot_enrich_ML_kidney <- enrich_summary %>%
  filter(
    animal == "ML",
    category %in% c("Kidney_PresentLastK", "Kidney_NotLastK")
  ) %>%
  mutate(
    category = recode(category,
                      "Kidney_PresentLastK" = "Late Shed (Kidney)",
                      "Kidney_NotLastK"    = "Non-Late Shed (Kidney)")
  ) %>%
  mutate(category = factor(
    category,
    levels = c("Late Shed (Kidney)", "Non-Late Shed (Kidney)")
  ))

p_frac_cat_ML_kidney <- ggplot(plot_enrich_ML_kidney,
                               aes(x = category, y = frac_in_cat, fill = group)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  labs(
    title = paste0("Top-", SMOLDER_TOP_N, " Smolderers: Kidney Categories"),
    x = "Kidney Category",
    y = "Fraction of Barcodes",
    fill = ""
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal(base_size = 14) +       # <-- bigger base font
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    legend.title = element_text(size = 14),
    plot.title   = element_text(size = 18, face = "bold"),
    panel.grid.major.x = element_blank()
  )

print(p_frac_cat_ML_kidney)

ggsave(
  filename = paste0("enrichment_TopSmolder_vs_Random_ML_KidneyOnly_Top", SMOLDER_TOP_N, "_styled.pdf"),
  plot     = p_frac_cat_ML_kidney,
  path     = "./",
  width    = 7, height = 7
)

ggsave(
  filename = paste0("enrichment_TopSmolder_vs_Random_ML_KidneyOnly_Top", SMOLDER_TOP_N, "_styled_colorblind.png"),
  plot     = p_frac_cat_ML_kidney,
  path     = "./",
  width    = 7,     # inches
  height   = 7,     # inches
  dpi      = 600    # high-quality PNG
)



# Early vs late smoldering and late winners


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
})


TOP_LATE_N        <- if (exists("TOP_N")) TOP_N else 50   
DETECT            <- if (exists("DETECT_THRESH")) DETECT_THRESH else 0.2
LATE_FRACTION     <- 0.25       # last quartile of timepoints
EARLY_SMOLDER_FRAC <- 0.5       # e.g. "smolderer" = detected in >= 50% of early samples
SAVE_DIR_LOCAL    <- if (exists("SAVE_DIR")) SAVE_DIR else "."


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
    # early
    n_early_detect = sum(detected & is_early, na.rm = TRUE),
    # late
    n_late_detect  = sum(detected & is_late,  na.rm = TRUE),
    late_total_load = sum(if_else(is_late, bc_day_level, 0), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(denom_tbl, by = "animal") %>%
  mutate(
    pct_detect_early = if_else(n_early_tp > 0, 100 * n_early_detect / n_early_tp, NA_real_),
    pct_detect_late  = if_else(n_late_tp  > 0, 100 * n_late_detect  / n_late_tp,  NA_real_)
  )


topLate_smolder <- smolder_early_late %>%
  group_by(animal) %>%
  arrange(desc(pct_detect_late), desc(late_total_load), .by_group = TRUE) %>%
  slice_head(n = TOP_LATE_N) %>%
  ungroup() %>%
  transmute(animal, barcode, group = "TopLateSmolder")

topLate_highshed <- smolder_early_late %>%
  group_by(animal) %>%
  arrange(desc(late_total_load), .by_group = TRUE) %>%
  slice_head(n = TOP_LATE_N) %>%
  ungroup() %>%
  transmute(animal, barcode, group = "TopLateHighShed")

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
  transmute(animal, barcode, group = "Random")


all_sets <- bind_rows(topLate_smolder, topLate_highshed, random_late) %>%
  left_join(smolder_early_late, by = c("animal","barcode"))

all_sets <- all_sets %>%
  mutate(
    early_smolder = !is.na(pct_detect_early) &
      pct_detect_early >= (EARLY_SMOLDER_FRAC * 100)
  )


# How much did late smolderers shed early?


p_early_pct <- all_sets %>%
  filter(!is.na(pct_detect_early)) %>%
  mutate(group = factor(group, levels = c("Random", "TopLateHighShed", "TopLateSmolder"))) %>%
  ggplot(aes(x = group, y = pct_detect_early, fill = group)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.6) +
  geom_boxplot(width = 0.12, outlier.shape = NA) +
  geom_jitter(width = 0.08, height = 0, size = 1, alpha = 0.4) +
  facet_wrap(vars(animal)) +
  labs(
    title = paste0("Early smoldering of Top-", TOP_LATE_N,
                   " late smolderers and late high-shedders"),
    x = "",
    y = "Pct of EARLY urine timepoints with detection (≥ threshold)"
  ) +
  scale_fill_manual(values = c(Random = "salmon", TopLateHighShed = "goldenrod", TopLateSmolder = "deepskyblue3")) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

print(p_early_pct)

ggsave(
  filename = paste0("fig8A_early_pct_TopLateSmolder_vs_LateHigh_vs_Random_Top", TOP_LATE_N, ".pdf"),
  plot     = p_early_pct,
  path     = SAVE_DIR_LOCAL,
  width    = 10, height = 6
)



early_flag_summary <- all_sets %>%
  group_by(animal, group) %>%
  summarise(
    n_barcodes = n(),
    n_early_smolder = sum(early_smolder),
    frac_early_smolder = n_early_smolder / pmax(n_barcodes, 1),
    .groups = "drop"
  )

p_early_flag <- early_flag_summary %>%
  mutate(group = factor(group, levels = c("Random", "TopLateHighShed", "TopLateSmolder"))) %>%
  ggplot(aes(x = group, y = frac_early_smolder, fill = group)) +
  geom_col(width = 0.6) +
  facet_wrap(vars(animal)) +
  labs(
    title = paste0("Fraction of barcodes that are EARLY smolderers\n(threshold = ",
                   EARLY_SMOLDER_FRAC * 100, "% of early timepoints)"),
    x = "",
    y = "Fraction of barcodes classified as early smolderers"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c(Random = "salmon", TopLateHighShed = "goldenrod", TopLateSmolder = "deepskyblue3")) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid.major.x = element_blank()
  )

print(p_early_flag)

ggsave(
  filename = paste0("fig8B_frac_early_smolder_TopLateSmolder_vs_LateHigh_vs_Random_Top", TOP_LATE_N, ".pdf"),
  plot     = p_early_flag,
  path     = SAVE_DIR_LOCAL,
  width    = 10, height = 6
)



print(early_flag_summary)






library(dplyr)
library(tidyr)
library(tibble)
library(corrplot)



kidney_present_all <- tissue_bc_levels %>%
  filter(is_kidney(organ)) %>%                 # only kidney tissues
  inner_join(present_lastK, by = c("animal","barcode")) %>%
  distinct(animal, barcode)

kidney_notLastK_all <- tissue_bc_levels %>%
  filter(is_kidney(organ)) %>%
  inner_join(absent_lastK, by = c("animal","barcode")) %>%
  distinct(animal, barcode)


#Spearman correlation with kidney split
compute_spearman_correlation_splitKidney <- function(data,
                                                     animal_id,
                                                     kidney_present_tbl,
                                                     kidney_notLastK_tbl) {
  #kidney barcode lists for this mouse (ALL, not top-N)
  present_bc <- kidney_present_tbl %>%
    filter(animal == animal_id) %>%
    distinct(barcode)
  
  notLast_bc <- kidney_notLastK_tbl %>%
    filter(animal == animal_id) %>%
    distinct(barcode)
  
  bc_union <- bind_rows(present_bc, notLast_bc) %>%
    distinct(barcode)
  
  if (nrow(bc_union) == 0) {
    stop(paste("No kidney-present or kidney-notLastK barcodes for animal", animal_id))
  }
  
  animal_data <- data %>%
    filter(animal == animal_id) %>%
    semi_join(bc_union, by = "barcode") %>%
    mutate(
      organ = case_when(
        organ == "Kidney" & barcode %in% present_bc$barcode   ~ "Kidney_presentLastK",
        organ == "Kidney" & barcode %in% notLast_bc$barcode  ~ "Kidney_notLastK",
        TRUE                                                  ~ organ
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
  
  cor_matrix <- cor(as.matrix(animal_data),
                    method = "spearman",
                    use    = "pairwise.complete.obs")
  
  colnames(cor_matrix) <- gsub(" ", "\n", colnames(cor_matrix))
  rownames(cor_matrix) <- gsub(" ", "\n", rownames(cor_matrix))
  
  return(cor_matrix)
}

generate_correlation_plot <- function(cor_matrix, animal_id, animal_color) {
  corrplot.mixed(
    cor_matrix,
    upper      = "circle",
    lower      = "number",
    tl.col     = "black",
    tl.font    = 2,
    tl.cex     = .6,
    tl.srt     = 45,
    number.cex = 0.5,
    lower.col  = "black",
    is.corr    = TRUE,
    upper.col  = colorRampPalette(c("gray", "white", animal_color))(400),
    col.lim    = c(-1, 1)
  )
}


pdf("../plots/urine_tissue_manuscript4/Fig8A_kidneySplit_correlation_by_animal_ALL.pdf",
    onefile = TRUE, height = 10, width = 10)

par(mfrow = c(2, 2), mar = c(2, 2, 4, 2))

for (animal_id in names(animal_colors)) {
  cor_matrix <- compute_spearman_correlation_splitKidney(
    data               = tissue_bc_levels,
    animal_id          = animal_id,
    kidney_present_tbl = kidney_present_all,
    kidney_notLastK_tbl= kidney_notLastK_all
  )
  generate_correlation_plot(cor_matrix, animal_id, animal_colors[animal_id])
  title(main = animal_id, line = 3, cex.main = 1.2)
}

par(mfrow = c(1, 1))
dev.off()




library(corrplot)
library(dplyr)
library(tidyr)
library(tibble)

compute_mouse_corr_for_group <- function(tissue_tbl, group_tbl) {
  kidney_dat <- tissue_tbl %>%
    filter(is_kidney(organ)) %>%  # only kidney data
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

pdf("../plots/urine_tissue_manuscript4/Fig8b_kidney_between_mice_corr_ALL.pdf",
    width = 7, height = 7, onefile = TRUE)


cor_kidney_absent_mice <- compute_mouse_corr_for_group(
  tissue_tbl = tissue_bc_levels,
  group_tbl  = kidney_notLastK_all
)

corrplot.mixed(
  cor_kidney_absent_mice,
  upper      = "circle",
  lower      = "number",
  tl.col     = "black",
  tl.cex     = 0.9,
  number.cex = 0.8,
  lower.col  = "black",
  is.corr    = TRUE,
  upper.col  = colorRampPalette(c("gray", "white", "steelblue"))(400),
  col.lim    = c(-1, 1)
)
title(main = "Kidney_absent (NOT last-K)\nSpearman correlation between mice",
      line = 2, cex.main = 1.3)


cor_kidney_present_mice <- compute_mouse_corr_for_group(
  tissue_tbl = tissue_bc_levels,
  group_tbl  = kidney_present_all
)

corrplot.mixed(
  cor_kidney_present_mice,
  upper      = "circle",
  lower      = "number",
  tl.col     = "black",
  tl.cex     = 0.9,
  number.cex = 0.8,
  lower.col  = "black",
  is.corr    = TRUE,
  upper.col  = colorRampPalette(c("gray", "white", "darkred"))(400),
  col.lim    = c(-1, 1)
)
title(main = "Kidney_present (last-K)\nSpearman correlation between mice",
      line = 2, cex.main = 1.3)

dev.off()






library(dplyr)
library(stringr)
library(forcats)
library(ggplot2)


#GC% for ALL tissue barcodes (background)

gc_pct_all_barcodes <- tissue_bc_levels %>%
  distinct(barcode) %>%
  mutate(
    GC_pct = str_count(barcode, "[GC]") / str_length(barcode),
    animal = "All Barcodes"
  )


# 2) GC% for Top 50 kidney-present-lastK per mouse

gc_kidney_present_top50 <- kidney_rank %>%
  group_by(animal) %>%
  arrange(kidney_rank, .by_group = TRUE) %>%
  slice_head(n = 50) %>%
  ungroup() %>%
  distinct(animal, barcode) %>%
  mutate(GC_pct = str_count(barcode, "[GC]") / str_length(barcode))

gc_present_plot_df <- gc_kidney_present_top50 %>%
  bind_rows(gc_pct_all_barcodes) %>%
  mutate(
    animal = fct_relevel(animal, c("FL", "FR", "ML", "MR", "All Barcodes"))
  )

p_gc_kidney_present <- gc_present_plot_df %>%
  ggplot(aes(animal, GC_pct)) +
  geom_boxplot(aes(color = animal), linewidth = 0.8, show.legend = FALSE) +
  labs(
    title = "GC content — Top 50 kidney-present (last-K) vs all tissue barcodes",
    x     = "Mouse / All Barcodes",
    y     = "GC Content"
  ) +
  scale_color_manual(values = c(animal_colors, `All Barcodes` = "gray70")) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  theme_minimal()

p_gc_kidney_present


# 3) GC% for Top 50 kidney-absent-lastK per mouse

gc_kidney_absent_top50 <- kidney_rank_notLastK %>%
  group_by(animal) %>%
  arrange(kidney_rank, .by_group = TRUE) %>%
  slice_head(n = 50) %>%
  ungroup() %>%
  distinct(animal, barcode) %>%
  mutate(GC_pct = str_count(barcode, "[GC]") / str_length(barcode))

gc_absent_plot_df <- gc_kidney_absent_top50 %>%
  bind_rows(gc_pct_all_barcodes) %>%
  mutate(
    animal = fct_relevel(animal, c("FL", "FR", "ML", "MR", "All Barcodes"))
  )

p_gc_kidney_absent <- gc_absent_plot_df %>%
  ggplot(aes(animal, GC_pct)) +
  geom_boxplot(aes(color = animal), linewidth = 0.8, show.legend = FALSE) +
  labs(
    title = "GC content — Top 50 kidney-absent (NOT last-K) vs all tissue barcodes",
    x     = "Mouse / All Barcodes",
    y     = "GC Content"
  ) +
  scale_color_manual(values = c(animal_colors, `All Barcodes` = "gray70")) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  theme_minimal()

p_gc_kidney_absent


ggsave(
  "../plots/urine_tissue_manuscript2/GC_top50_kidney_present_vs_all.pdf",
  p_gc_kidney_present, width = 6, height = 4
)

ggsave(
  "../plots/urine_tissue_manuscript2/GC_top50_kidney_absent_vs_all.pdf",
  p_gc_kidney_absent, width = 6, height = 4
)




###make fasta for miRNA targets


library(dplyr)
library(purrr)

# 1) Get unique kidney barcodes with bc_level > 0
kidney_barcodes_unique <- tissue_bc_levels %>%
  filter(is_kidney(organ), bc_level > 0) %>%      # keep only kidney rows with >0
  distinct(barcode) %>%                            # unique barcode sequences
  arrange(barcode) %>%                             # optional: sort for reproducibility
  mutate(id = paste0("bc", row_number()))          # bc1, bc2, ...

# 2) Build FASTA lines
fasta_lines <- map2_chr(
  kidney_barcodes_unique$id,
  kidney_barcodes_unique$barcode,
  ~ paste0(">", .x, "\n", .y)
)

# 3) Write to a single FASTA file (change path/name if you want)
output_fasta <- "../plots/urine_tissue_manuscript3/kidney_barcodes_unique.fa"
output_fasta <- "./kidney_barcodes_unique.fa"

writeLines(fasta_lines, con = output_fasta)

output_fasta




########miRNA-mRNA match

RNAhybrid -s 3utr_human -t kidney_barcodes_unique.fa -q mmu.mirna.cleaned.fa -b 1 -m 10000 > rnahybrid.kidney_mirna.txt




library(tidyverse)






library(tidyverse)

rnahybrid_tsv <- "/stor/work/Sullivan/anik/barcode_project/data/mirna_target_kidney/analysis.rnahybrid/rev_rnahybrid/true.muPyV_rev_rnahybrid.kidney_mirna.parsed.tsv" #true_rev_rnahybrid.kidney_mirna.parsed.tsv" #"        
fasta_with_rev <- "/stor/work/Sullivan/anik/barcode_project/data/mirna_target_kidney/analysis.rnahybrid/rev_rnahybrid/kidney_barcodes_unique_with_rev.fa"
LAST_K <- 2

is_kidney <- function(x) grepl("kidney|renal", x, ignore.case = TRUE)


# 1) Load RNAhybrid results → miRNA_interaction

miRNA_interaction <- readr::read_tsv(rnahybrid_tsv, show_col_types = FALSE)
miRNA_interaction

library(dplyr)
library(stringr)

miRNA_interaction <- miRNA_interaction %>%
  mutate(
    mirna_number = as.integer(str_extract(mirna, "\\d+"))  # grab first run of digits
  ) %>%
  filter(
    !is.na(mirna_number),
    mirna_number >= 1,
    mirna_number <= 800
  ) 

# 
# #%>%
#   filter(
#     mfe <= -10,
#     p_value <=0.05
#   )



miRNA_interaction
# Normalize target IDs and keep forward barcode seq (original orientation)
miRNA_interaction <- miRNA_interaction %>%
  mutate(barcode_id_raw = target,
         barcode_id_norm = sub("-rev$", "", barcode_id_raw))


read_fasta_named <- function(path) {
  lines <- readLines(path)
  hdr_idx <- grep("^>", lines)
  ends <- c(hdr_idx[-1] - 1, length(lines))
  tibble(
    name = sub("^>", "", lines[hdr_idx]),
    seq  = map2_chr(hdr_idx + 1, ends, ~ gsub("\\s+", "", paste(lines[.x:.y], collapse = "")))
  )
}
fasta_tbl <- read_fasta_named(fasta_with_rev)


forward_seq_tbl <- fasta_tbl %>%
  mutate(forward_name = sub("-rev$", "", name)) %>%
  group_by(forward_name) %>%
  summarise(forward_seq = seq[which.min(grepl("-rev$", name))], .groups = "drop")

miRNA_interaction <- miRNA_interaction %>%
  left_join(forward_seq_tbl, by = c("barcode_id_norm" = "forward_name"))


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

kidney_present_set <- kidney_barcodes_all %>%
  semi_join(present_lastK, by = c("animal","barcode")) %>%
  mutate(kidney_status = "present_lastK")

kidney_absent_set <- kidney_barcodes_all %>%
  anti_join(present_lastK, by = c("animal","barcode")) %>%
  mutate(kidney_status = "absent_lastK")

kidney_status_tbl <- bind_rows(kidney_present_set, kidney_absent_set)


mirna_hit_barcodes <- miRNA_interaction %>%
  distinct(forward_seq) %>%
  rename(barcode = forward_seq)

mirna_hit_in_kidney <- kidney_status_tbl %>%
  semi_join(mirna_hit_barcodes, by = "barcode")

#rqandom control

set.seed(42)

# how many miRNA-hit barcodes in kidney for each mouse?
hits_per_mouse <- mirna_hit_in_kidney %>%
  count(animal, name = "n_hits")

# build random control: for each mouse, sample n_hits barcodes from all kidney barcodes
rand_control_by_mouse <- kidney_barcodes_all %>%
  group_split(animal) %>%
  purrr::map_dfr(function(df_mouse) {
    this_animal <- df_mouse$animal[1]
    n_hits <- hits_per_mouse %>%
      filter(animal == this_animal) %>%
      pull(n_hits)
    
    n_hits <- ifelse(length(n_hits) == 0 || is.na(n_hits) || n_hits <= 0, 0, n_hits)
    
    if (n_hits == 0) {
      return(tibble(animal = this_animal, barcode = character(0)))
    }
    
    df_mouse %>%
      slice_sample(n = min(n_hits, nrow(df_mouse))) %>%
      select(animal, barcode)
  })

# map those random barcodes into present/absent status table
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

hit_summary   <- summarize_status(mirna_hit_in_kidney, "miRNA_hit")
rand_summary  <- summarize_status(rand_in_kidney_status, "random_match")

enrichment_summary_per_mouse <- bind_rows(hit_summary, rand_summary) %>%
  mutate(kidney_status = factor(kidney_status,
                                levels = c("present_lastK","absent_lastK"))) %>%
  arrange(animal, group, kidney_status)

# Overall roll-up
overall_summary <- bind_rows(
  mirna_hit_in_kidney %>% mutate(group = "miRNA_hit"),
  rand_in_kidney_status %>% mutate(group = "random_match")
) %>%
  count(group, kidney_status, name = "n") %>%
  group_by(group) %>%
  mutate(
    n_total = sum(n),
    prop    = ifelse(n_total > 0, n / n_total, NA_real_)
  ) %>%
  ungroup() %>%
  mutate(kidney_status = factor(kidney_status,
                                levels = c("present_lastK","absent_lastK")))


enrichment_summary_per_mouse
overall_summary


ggplot(enrichment_summary_per_mouse,
       aes(x = kidney_status, y = prop, fill = group)) +
  geom_col(position = position_dodge(width = 0.6)) +
  facet_wrap(vars(animal), nrow = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Kidney last-K status of miRNA-hit barcodes vs random (per mouse)",
       x = "", y = "Proportion of barcodes") +
  theme_minimal()

#fisher
contingency <- overall_summary %>%
  select(group, kidney_status, n) %>%
  tidyr::pivot_wider(names_from = kidney_status, values_from = n, values_fill = 0) %>%
  tibble::column_to_rownames("group") %>%
  as.matrix()

if (all(c("present_lastK","absent_lastK") %in% colnames(contingency)) &&
    all(c("miRNA_hit","random_match") %in% rownames(contingency))) {
  fisher_res <- fisher.test(contingency)
  print(contingency)
  print(fisher_res)
} else {
  message("Contingency table missing categories; skipping Fisher test.")
}

###length distribution




library(dplyr)
library(stringr)
library(ggplot2)

topN_length_dist <- bind_rows(
  kidney_present_topN %>%
    mutate(kidney_status = "present_lastK"),
  kidney_absent_topN %>%
    mutate(kidney_status = "absent_lastK")
) %>%
  mutate(barcode_length = str_length(barcode)) %>%   # barcode is the sequence
  count(animal, kidney_status, barcode_length, name = "n")

topN_length_dist

ggplot(topN_length_dist,
       aes(x = barcode_length, y = n, fill = kidney_status)) +
  geom_col(position = position_dodge(width = 0.7)) +
  facet_wrap(vars(animal), nrow = 1) +
  scale_x_continuous(breaks = sort(unique(topN_length_dist$barcode_length))) +
  labs(
    title = sprintf("Length distribution of Top-%d kidney barcodes\npresent vs absent in last-K urine", TOP_N_KIDNEY),
    x = "Barcode length (nt)",
    y = "Count of barcodes"
  ) +
  scale_fill_manual(
    values = c(present_lastK = "#1b9e77", absent_lastK = "#d95f02"),
    name   = "Kidney status"
  ) +
  theme_minimal()





