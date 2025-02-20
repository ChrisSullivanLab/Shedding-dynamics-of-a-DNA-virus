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

animal_colors = c('FL' = 'palevioletred3', 'FR' = 'palevioletred1', 'ML' = 'skyblue3', 'MR' = 'skyblue1')
rgb_codes <- col2rgb(animal_colors)
rgb_codes
# palevioletred3: RGB(205, 98, 145)
# palevioletred1: RGB(255, 130, 171)
# skyblue3: RGB(80, 170, 255)
# skyblue1: RGB(135, 206, 255)
# List all relevant files in the current directory
stock_raw_files = list.files(pattern = 'BC7[_b].*cutadapt_counts.txt', full.names = TRUE)

stock_raw_counts = stock_raw_files %>%
  map(~read_table(., col_names = c(str_extract(., 'BC7[^_]*'), 'barcode'))) %>%
  map(~mutate(., barcode = replace_na(barcode, ''))) %>%
  reduce(full_join) %>%
  replace(is.na(.), 0) %>%
  pivot_longer(-barcode, names_to = 'BC7_run', values_to = 'count')

stock_raw_counts %<>%
  group_by(barcode) %>%
  summarize(count = sum(count))

stock_raw_counts %>%
  filter(barcode != '') %>%
  write_tsv('BC7_combined_stock_counts.tsv', col_names = FALSE)

stock_L3_counts = read_tsv('BC7_combined_stock_mp_L3_clusters.tsv',
                           col_names = c('barcode', 'count', 'elements'),
                           progress = FALSE) %>%
  select(-elements)

urine_raw_counts = list.files('/stor/work/Sullivan/anik/barcode_project/data/sra/counts3', # Update path to the directory
                              pattern='.*Urine.*\\.fastq\\.counts\\.txt$', # Match all relevant files
                              full.names=TRUE) %>%
  set_names(., gsub('^[^_]*_[^_]*_(.*)\\.fastq\\.counts\\.txt$', '\\1', ., perl=TRUE)) %>%
  imap(~ read_table(.x, col_names=c(paste0('sample_', .y), 'barcode'), col_types='ic')) %>%
  reduce(full_join, by='barcode') %>%
  pivot_longer(-barcode, names_to='run_sample', values_to='count') %>%
  replace_na(list(count = 0)) %>%
  mutate(sample = gsub('^sample_(.*)$', '\\1', run_sample)) %>%
  group_by(sample, barcode) %>%
  summarize(count = sum(count), .groups = 'drop') %>%
  ungroup()




#####

tissue_raw_counts = list.files('/stor/work/Sullivan/anik/barcode_project/data/sra/counts3', # Directory containing all files
                               pattern='.*(Brain|Bladder|Lung|Liver|Kidney|Heart|Gut|Gonads|Whole|Testicles|Spleen|Salivary|Muscle).*\\.fastq\\.counts\\.txt$', 
                               full.names=TRUE) %>%
  set_names(., gsub('^[^_]*_[^_]*_(.*)\\.fastq\\.counts\\.txt$', '\\1', ., perl=TRUE)) %>%
  imap(~ read_table(.x, col_names=c(paste0('sample_', .y), 'barcode'), col_types='ic')) %>%
  reduce(full_join, by='barcode') %>%
  select(-contains("Plasma")) %>%
  pivot_longer(-barcode, names_to='run_sample', values_to='count') %>%
  replace_na(list(count = 0)) %>%
  mutate(sample = gsub('^sample_(.*)$', '\\1', run_sample)) %>%
  group_by(sample, barcode) %>%
  summarize(count = sum(count), .groups = 'drop') %>%
  ungroup()

urine_raw_counts %>%
  distinct(sample) %>%
  pull(sample)

tissue_raw_counts %>%
  distinct(sample) %>%
  pull(sample)


# --- Read L3 Clustered Stock, Apply Cutoff -------------------------------

clustered_combined_stock_barcodes = read_tsv('BC7_combined_stock_mp_L3_clusters.tsv',
                                             col_names = c('barcode', 'count', 'elements'),
                                             progress = FALSE) %>%
  select(-elements)

cutoff_99pct_stock_barcodes = clustered_combined_stock_barcodes %>%
  mutate(total = sum(count),
         pct = count / total,
         cum_pct = cumsum(pct)) %>%
  filter(cum_pct <= 0.99)  %>%
  filter(!(barcode %in% c('A', 'TTT')))

# Specify the output file path
output_file <- "cutoff_99pct_stock_barcodes.tsv"

# Save the variable as a TSV file
write.table(
  cutoff_99pct_stock_barcodes,  
  file = "cutoff_99pct_stock_barcodes.tsv",          
  sep = "\t",                   
  row.names = FALSE,            
  col.names = FALSE,            
  quote = FALSE               
)

cat("Data has been saved to", output_file)

# --- Compute Distances Between Sample And Stock Barcodes -----------------
## Cut off distances past threshold, choose smallest distances to stock barcode
if (!dir.exists("split_sample_barcodes")) {
  dir.create("split_sample_barcodes")
}
cutoff_99pct_stock_barcodes %>%
  pull(barcode) %>%
  saveRDS('split_sample_barcodes/stock_cum99pct_cutoff.rds')

## Split unique urine barcodes into smaller pieces to parallelize adist computation
## This splitting factor splits into list of 10,000 barcode vectors
## Split unique urine barcodes into smaller pieces to parallelize adist computation
## This splitting factor splits into list of 10,000 barcode vectors

unique_sample_barcodes = c(urine_raw_counts$barcode, tissue_raw_counts$barcode) %>%
  unique()

sample_splitting_factor = ceiling(seq(1, length(unique_sample_barcodes)) / 10000)

split(unique_sample_barcodes, sample_splitting_factor) %>%
  imap(~ saveRDS(.x, paste0('split_sample_barcodes/sample_barcodes_', .y, '.rds')))

# Distance to stock barcode threshold
threshold = 3

getwd()
readRDS('./split_sample_barcodes/sample_barcodes_1.rds') %>% head(20)
readRDS('./split_sample_barcodes/stock_cum99pct_cutoff.rds') %>% head(20)


# Define the path to stock barcodes file
stock_barcodes <- readRDS("./split_sample_barcodes/stock_cum99pct_cutoff.rds")

# Get the list of all sample barcode files in the split_sample_barcodes folder
sample_files <- list.files("./split_sample_barcodes", pattern = "sample_barcodes_.*\\.rds", full.names = TRUE)



# Function to calculate distances and filter based on the threshold
process_sample_file <- function(sample_file) {
  sample_barcodes <- readRDS(sample_file)
  distances <- expand.grid(sample = sample_barcodes, stock = stock_barcodes) %>%
    mutate(dist = stringdist::stringdist(sample, stock, method = "lv")) %>%
    filter(dist <= threshold)
    output_file <- gsub("sample_barcodes_", "dist_", sample_file)
  saveRDS(distances, output_file)
}

walk(sample_files, process_sample_file)

cat("Distance files have been generated successfully in the ./split_sample_barcodes/ directory!")

sample_barcodes_within_threshold_distance = list.files('split_sample_barcodes',
                                                       pattern = '.*dist_[[:digit:]]+\\.rds',
                                                       full.names = TRUE) %>%
  map(~ readRDS(.) %>% filter(., dist <= threshold)) %>%
  map(~ group_by(., sample) %>% filter(., dist == min(dist))) %>%
  reduce(bind_rows) %>%
  mutate(n = n()) %>%
  rename(sample_bc = sample, stock_bc = stock)


## Associate stock barcode to every individual sample barcode, then divide sample count by the number of
## stock barcodes equidistance from the sample barcode. Suppose a sample barcode has count 3, and there are
## three stock barcodes equidistance (with minimum distance) from sample barcode; then 1 is assigned to each
## stock barcode from the original sample count of 3.

urine_to_stock_barcodes = inner_join(urine_raw_counts,
                                     sample_barcodes_within_threshold_distance,
                                     by = c('barcode' = 'sample_bc')) %>%
  mutate(count_div_n = count / n) %>%
  select(sample, stock_bc, count_div_n) %>%
  group_by(sample, stock_bc) %>%
  summarize(weighted_count = sum(count_div_n)) # Sum the weighted values of the same stock barcodes

#urine_to_stock_barcodes <- urine_to_stock_barcodes %>% filter(weighted_count >= 10)
urine_to_stock_barcodes <- urine_to_stock_barcodes %>%
  mutate(weighted_count = if_else(weighted_count < 10, 0, weighted_count))

urine_to_stock_barcodes %>%
  distinct(sample) %>%
  pull(sample)
## ^ In above code block, many distinct sample barcodes were converted to the same stock barcodes,
## so we need to sum those up, which the group_by and mutate do.

tissue_to_stock_barcodes = inner_join(tissue_raw_counts,
                                      sample_barcodes_within_threshold_distance,
                                      by = c('barcode' = 'sample_bc')) %>%
  mutate(count_div_n = count / n) %>%
  select(sample, stock_bc, count_div_n) %>%
  group_by(sample, stock_bc) %>%
  summarize(weighted_count = sum(count_div_n)) # Sum the weighted values of the same stock barcodes

tissue_to_stock_barcodes <- tissue_to_stock_barcodes %>%
  mutate(weighted_count = if_else(weighted_count < 10, 0, weighted_count))

#tissue_to_stock_barcodes <- tissue_to_stock_barcodes %>% filter(weighted_count >= 10)


tissue_to_stock_barcodes %>%
  distinct(sample) %>%
  pull(sample)

# --- Urine Analysis ------------------------------------------------------

## Read in data on concentrations of barcodes in urine,
## calculate

urine_pcr <- read_tsv('urine_pcr.tsv') %>%
  select(-starts_with('Copies'), -Comments) %>%
  rename(
    Days_pi = Days_pi,
    animal = animal,
    sample = sample,
    mean_ul_urine = `Quantity_Mean/ul_urine`,
    genomes_reaction = `Genomes/reaction`,
    quantity_mean = Quantity_Mean,
    vol_urine_pcr = `Total_Vol_(ul)_used_in_PCR`
  ) %>%
  select(sample, animal, Days_pi, mean_ul_urine, everything())

urine_pcr2 %>%
  distinct(sample) %>%
  pull(sample)


urine_pcr2 <- urine_pcr %>%
  mutate(sample2 = paste0(animal, "_Urine_Day_", Days_pi))

head(urine_pcr2 %>% select(sample, sample2))

# Rename columns in urine_pcr2
urine_pcr2 <- urine_pcr2 %>%
  rename(
    sample_backup = sample,  
    sample = sample2       
  )

colnames(urine_pcr2)

# Process Barcodes with weighted count data
urine_bc_levels <- urine_to_stock_barcodes %>%
  group_by(sample) %>%
  mutate(frac_weighted_count = weighted_count / sum(weighted_count)) %>%
  ungroup() %>%
  left_join(urine_pcr2, by = 'sample') %>%
  mutate(bc_level = frac_weighted_count * mean_ul_urine) %>%
  select(sample, animal, Days_pi, stock_bc, bc_level, frac_weighted_count, mean_ul_urine) %>%
  rename(barcode = stock_bc)

urine_bc_levels %>%
  distinct(sample) %>%
  pull(sample)

dim(urine_bc_levels)

print(head(urine_bc_levels))



## -- Common Functions ----------------------------------------------------

top_barcodes_by_max <- function(barcode_table, top_n) {
  top_barcode_levels = barcode_table %>%
    group_by(animal, barcode) %>%      # Get highest level for each barcode across days
    mutate(max_level = max(bc_level)) %>%
    ungroup() %>%
    select(-frac_weighted_count) %>%
    nest_by(animal, max_level) %>%     # Nest all the longitudinal data to focus on max levels
    group_by(animal) %>%               # Get rid of rowwise grouping
    slice_max(max_level, n = top_n) %>%   # Get max 10 bc levels per animal
    ungroup() %>% 
    unnest(cols = data)                # Expand all the data that corresponds to just the max levels
  return(top_barcode_levels)
}

top_barcodes_by_frac_max <- function(barcode_table, top_n) {
  top_barcode_levels = barcode_table %>%
    group_by(animal, barcode) %>%      # Get highest level for each barcode across days
    mutate(max_level = max(frac_weighted_count)) %>%
    ungroup() %>%
    select(-frac_weighted_count) %>%
    nest_by(animal, max_level) %>%     # Nest all the longitudinal data to focus on max levels
    group_by(animal) %>%               # Get rid of rowwise grouping
    slice_max(max_level, n = top_n) %>%   # Get max 10 bc levels per animal
    ungroup() %>% 
    unnest(cols = data)                # Expand all the data that corresponds to just the max levels
  return(top_barcode_levels)
}

slice_barcodes_by_max <- function(barcode_table, top_rank, bottom_rank) {
  top_barcode_levels = barcode_table %>%
    group_by(animal, barcode) %>%      # Get highest level for each barcode across days
    mutate(max_level = max(bc_level)) %>%
    ungroup() %>%
    select(-frac_weighted_count) %>%
    nest_by(animal, barcode, max_level) %>%     # Nest all the longitudinal data to focus on max levels
    group_by(animal) %>%               # Get rid of rowwise grouping
    arrange(desc(max_level)) %>%
    slice(top_rank:bottom_rank) %>%   # Get max 10 bc levels per animal
    ungroup() %>% 
    unnest(cols = data)                # Expand all the data that corresponds to just the max levels
  return(top_barcode_levels)
}



## -- Urine Area Plots Of Top 10 Winners ----------------------------------

individual_winner_area_plots = urine_bc_levels %>%
  top_barcodes_by_max(10) %>%
  split(.$animal) %>%
  imap(~ mutate(., barcode = fct_reorder(as.factor(barcode), max_level, .desc=TRUE))) %>%
  imap(~ ggplot(.x, aes(x = Days_pi, fill = barcode)) +
         geom_area(aes(y = mean_ul_urine), fill = 'grey85') +
         geom_area(aes(y = bc_level)) +
         facet_wrap(vars(barcode), ncol = 1) +
         labs(#title = 'Top 10 Winners',
           subtitle = paste(.y),
           x = 'Days post-injection',
           y = 'Barcode level') +
         theme_minimal() +
         theme(legend.position = 'none'))

summed_winner_area_plots = urine_bc_levels %>%
  top_barcodes_by_max(10) %>%
  group_by(animal, Days_pi, mean_ul_urine) %>%
  summarize(sum_winner_barcodes = sum(bc_level))  %>%
  split(.$animal) %>%
  imap(~ ggplot(.x, aes(x = Days_pi)) +
         geom_area(aes(y = mean_ul_urine), fill = 'grey85') +
         geom_area(aes(y = sum_winner_barcodes), fill = 'navajowhite3') +
         labs(title = .y, x='', y=''
              #x = 'Days post-injection', y = 'Barcode level'
         ) +
         theme_minimal() +
         theme(legend.position = 'none'))


pdf("../plots/urine_tissue_manuscript2/Fig7_sum_pf_top_10_winners.pdf", onefile = TRUE)

grid.arrange(
  summed_winner_area_plots[["FL"]], 
  summed_winner_area_plots[["FR"]],
  summed_winner_area_plots[["ML"]], 
  summed_winner_area_plots[["MR"]],
  nrow = 4, ncol = 1
  #top = "Sum of Top 10 Winner Barcodes for Each Animal" 
)

#grid.text("Days Post-Infection", x = 0.5, y = -0.1, gp = gpar(fontsize = 12, col = "red"))
grid.text("Days Post-Injection", x = unit(.48, "npc"), y = unit(0.04, "npc"), just = c("left", "top"), gp = gpar(fontsize = 11)) #fontface = "bold"))
grid.text("Barcode level", x = unit(.005, "npc"), y = unit(0.5, "npc"), just = c("left", "top"), gp = gpar(fontsize = 11), rot = 90) #fontface = "bold"))

dev.off()

layout_matrix <- rbind(
  c(1, 2, 3, 4),  
  c(5, 6, 7, 8)  
)
pdf(file = "../plots/urine_tissue_manuscript2/Fig2s_ridge_plots_of_the_10_most_abundant_barcodes.pdf", width = 10, height = 12)
figs2<- grid.arrange(
  individual_winner_area_plots[["FL"]], 
  individual_winner_area_plots[["FR"]],
  individual_winner_area_plots[["ML"]], 
  individual_winner_area_plots[["MR"]],
  summed_winner_area_plots[["FL"]], 
  summed_winner_area_plots[["FR"]],
  summed_winner_area_plots[["ML"]], 
  summed_winner_area_plots[["MR"]],
  layout_matrix = layout_matrix,
  heights = c(4, 1)#,  # First row 4x height, second row 1x height
  #top = "Urine Barcode Plots: Individual (Row 1) and Summed (Row 2)"  # Add a title
)

grid.text("A", x = unit(.005, "npc"), y = unit(.99, "npc"), just = c("left", "top"), gp = gpar(fontsize = 18, fontface = "bold"))
grid.text("B", x = unit(.005, "npc"), y = unit(0.2, "npc"), just = c("left", "top"), gp = gpar(fontsize = 18, fontface = "bold"))

dev.off()



#-- Number Of Unique Barcodes -------------------------------------------

cleaned_urine_data <- urine_to_stock_barcodes %>%
  filter(weighted_count >= 10) %>% 
  distinct()


processed_urine_data <- cleaned_urine_data %>%
  mutate(
    animal = sub("_.*", "", sample), # Extract animal (e.g., FL, FR, etc.)
    day = as.numeric(sub(".*Day_", "", sample)) # Extract day as a number
  )

barcode_counts <- processed_urine_data %>%
  group_by(animal, day) %>%
  summarise(unique_barcodes = n_distinct(stock_bc), .groups = "drop") # Count unique barcodes


library(writexl)

write_xlsx(
  barcode_counts,
  path = "barcode_counts.xlsx")
library(ggplot2)
library(dplyr)

# Create bins with 100-unit ranges
binned_counts <- barcode_counts %>%
  mutate(barcode_range = cut(unique_barcodes,
                             breaks = seq(0, 3500, by = 100),
                             labels = paste0(seq(0, 3400, by = 100), "-", seq(100, 3500, by = 100)),
                             include.lowest = TRUE)) %>%
  count(barcode_range)

# Create bar plot
ggplot(binned_counts, aes(x = barcode_range, y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Distribution of Unique Barcodes by Range",
       x = "Number of Unique Barcodes (Range)",
       y = "Count of Occurrences") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))  # Prevent label overlap





random_folder <- "split_sample_barcodes/random/"

if (!dir.exists(random_folder)) {
  dir.create(random_folder, recursive = TRUE)
}

# Function to shuffle nucleotides in a barcode
shuffle_barcode <- function(barcode) {
  paste0(sample(unlist(strsplit(barcode, ""))), collapse = "")
}
sample_files <- list.files("split_sample_barcodes", pattern = "sample_barcodes_.*\\.rds", full.names = TRUE)

walk(sample_files, function(file) {
  sample_barcodes <- readRDS(file)
  randomized_barcodes <- map_chr(sample_barcodes, shuffle_barcode)
  saveRDS(randomized_barcodes, gsub("split_sample_barcodes", random_folder, file))
})


random_folder <- "split_sample_barcodes/random"
stock_barcodes <- readRDS("./split_sample_barcodes/stock_cum99pct_cutoff.rds")

randomized_files <- list.files(random_folder, pattern = "sample_barcodes_.*\\.rds", full.names = TRUE)

randomized_data <- randomized_files %>%
  map(~ readRDS(.) %>%
        tibble(barcode = .)) %>%
  set_names(basename(randomized_files)) %>%
  bind_rows(.id = "source_file")

process_random_file <- function(file) {
  sample_barcodes <- readRDS(file)
  distances <- expand.grid(sample = sample_barcodes, stock = stock_barcodes) %>%
    mutate(dist = stringdist::stringdist(sample, stock, method = "lv")) %>%
    filter(dist <= threshold)
  output_file <- gsub("sample_barcodes_", "dist_", file)
  saveRDS(distances, output_file)
}

walk(randomized_files , process_random_file)


randomized_distances <- list.files(random_folder, pattern = "dist_.*\\.rds", full.names = TRUE) %>%
  map(~ readRDS(.) %>% filter(dist <= threshold)) %>%
  map(~ group_by(., sample) %>% filter(dist == min(dist))) %>%
  reduce(bind_rows) %>%
  mutate(n = n()) %>%
  rename(sample_bc = sample, stock_bc = stock)

randomized_to_stock_barcodes <- inner_join(urine_raw_counts,
                                           randomized_distances,
                                           by = c('barcode' = 'sample_bc')) %>%
  mutate(count_div_n = count / n) %>%
  select(sample, stock_bc, count_div_n) %>%
  group_by(sample, stock_bc) %>%
  summarize(weighted_count = sum(count_div_n))

cleaned_randomized_data <- randomized_to_stock_barcodes %>%
  filter(weighted_count >= 10) %>%
  distinct()

processed_randomized_data <- cleaned_randomized_data %>%
  mutate(
    animal = sub("_.*", "", sample),
    day = as.numeric(sub(".*Day_", "", sample))
  )

barcode_counts_random <- processed_randomized_data %>%
  group_by(animal, day) %>%
  summarise(unique_barcodes = n_distinct(stock_bc), .groups = "drop")

# Merge original and randomized data for comparison
comparison_data <- full_join(barcode_counts, barcode_counts_random, by = c("animal", "day"), suffix = c("_original", "_random")) %>%
  mutate(difference = unique_barcodes_original - unique_barcodes_random)

create_animal_plots <- function(animal_name) {
  animal_data <- comparison_data %>% filter(animal == animal_name)
  
  p1 <- ggplot(animal_data, aes(x = day)) +
    geom_line(aes(y = unique_barcodes_original, color = animal_name), size = 1) +  # Map animal to color
    geom_line(aes(y = unique_barcodes_random, color = "Randomized"), size = 1, linetype = "dashed") +
    labs(
      title = paste(animal_name),
      #x = "Day",
      #y = "Unique Barcodes",
      x = "",
      y = ""
    ) +
    scale_color_manual(
      values = c(animal_colors, "Randomized" = "gray") 
    ) +
    theme_minimal() +
    theme(legend.title = element_blank(),
          plot.title = element_text(size = 10))
  
  p4 <- animal_data %>%
    mutate(percentage_change = (difference / unique_barcodes_original) * 100) %>%
    ggplot(aes(x = day, y = percentage_change, color = animal_name)) +
    geom_line(size = 1) +
    #geom_point(size = 2) +
    labs(
      title = paste(animal_name),
      # x = "Day",
      # y = "Percentage Change (%)",
      x = "",
      y = ""
    ) +
    scale_color_manual(values = animal_colors) +  # Use animal_colors for the lines
    theme_minimal() +
    theme(legend.position = "none",  plot.title = element_text(size = 10) )  # Hide legend since animal name is in the title
  
  list(p1 = p1, p4 = p4)
}

# Generate plots for each animal
animal_plots <- comparison_data %>%
  pull(animal) %>%
  unique() %>%
  map(~ create_animal_plots(.x))

# Extract Plot 1 and Plot 4 for each animal
plot_list <- map(animal_plots, ~ list(.x$p1, .x$p4)) 
plot_list <- map(animal_plots, ~ list(.x$p1)) 

# Flatten the list of plots
plot_list_flat <- flatten(plot_list)

combined_plot <- grid.arrange(
  grobs = plot_list_flat,  
  nrow = 4,               
  ncol = 1#2                 
)



pdf('../plots/urine_tissue_manuscript2/fig2_unique_barcodes_detected_at_each_time_point.pdf',
    width = 8, height = 8) 

grid.draw(combined_plot)
# grid.text("A", x = 0.20, y = .99, gp = gpar(fontsize = 16, fontface = "bold"), just = "top")  
# grid.text("B", x = 0.75, y = .99, gp = gpar(fontsize = 16, fontface = "bold"), just = "top")  
#grid.text("A", x = 0.01, y = .99, gp = gpar(fontsize = 16, fontface = "bold"), just = "top")  
#grid.text("B", x = 0.5, y = .99, gp = gpar(fontsize = 16, fontface = "bold"), just = "top")  
grid.text("Day", x = unit(.45, "npc"), y = unit(0.04, "npc"), just = c("left", "top"), gp = gpar(fontsize = 11)) #fontface = "bold"))
#grid.text("Day", x = unit(.75, "npc"), y = unit(0.04, "npc"), just = c("left", "top"), gp = gpar(fontsize = 11)) #fontface = "bold"))

grid.text("Unique Barcodes", x = unit(.01, "npc"), y = unit(0.45, "npc"), just = c("left", "top"), gp = gpar(fontsize = 11), rot = 90) #fontface = "bold"))
#grid.text("Enrichment Over Background (%)", x = unit(.5, "npc"), y = unit(0.4, "npc"), just = c("left", "top"), gp = gpar(fontsize = 11), rot = 90) #fontface = "bold"))

dev.off()


## -- Compare Stock And Day 1 Repertoires ---------------------------------

# Taken from https://stats.stackexchange.com/questions/31565/compute-a-cosine-dissimilarity-matrix-in-r
# But modified to handle the incoming table has samples as columns and observations as rows

cosine_similarity <- function(table) {
  matrix = as.matrix(table)
  sim = t(t(matrix) / sqrt(colSums(matrix * matrix)))
  return(t(sim) %*% sim)
}

first_day_and_stock_frac_counts = urine_bc_levels %>%
  filter(Days_pi == 1) %>%
  select(animal, barcode, frac_weighted_count) %>%
  rename(repertoire = animal, pct = frac_weighted_count) %>%
  bind_rows(
    cutoff_99pct_stock_barcodes %>%
      mutate(total = sum(count), pct = count / total) %>%
      select(barcode, pct) %>%
      mutate(repertoire = 'stock')
  ) %>%
  pivot_wider(names_from = repertoire, values_from = pct, values_fill = 0) %>%
  column_to_rownames(var = 'barcode')

pdf('../plots/urine_tissue_manuscript2/day1_and_stock_cosine_corrplot.pdf', onefile = FALSE)

cosine_similarity(first_day_and_stock_frac_counts) %>%
  corrplot.mixed(.,
                 upper = 'circle',
                 lower.col = 'grey',
                 tl.col = 'black',
                 is.corr = FALSE,
                 col.lim = c(0, max(.)),  # For some reason, explicitly setting upper limit to 1 gave error
                 upper.col = COL1('Purples', 200))
dev.off()

## -- Cosine Similarities Per Animal Across Time --------------------------

scale_cos_sim_to_mean_ul_urine <- function(table, animal) {
  table = table %>% filter(animal == animal)
  scale_factor = ceiling(max(table$mean_ul_urine) / 1000) * 1000  # Multiples of 1000 seem appropriate for this particular data
  
  return(scale_factor)
}

cos_sim_across_days = urine_bc_levels %>%
  select(animal, Days_pi, barcode, frac_weighted_count) %>%
  rename(pct = frac_weighted_count) %>%
  split(.$animal) %>%
  imap(~ select(.x, -animal)) %>%
  imap(~ arrange(.x, Days_pi)) %>%
  imap(~ pivot_wider(.x, names_from = Days_pi, values_from = pct, values_fill = 0)) %>%
  imap(~ column_to_rownames(.x, var = 'barcode')) %>%
  imap(~ cosine_similarity(.x))

cos_sim_subsequent_timepoints = urine_bc_levels %>%
  select(animal, Days_pi, mean_ul_urine) %>%
  split(.$animal) %>%
  imap(~ distinct(.x)) %>%
  imap(~ arrange(.x, Days_pi)) %>%
  imap(~ mutate(.x, second_timepoint = lead(Days_pi))) %>%
  imap(~ slice(.x, 1:n() -1)) %>%
  imap(~ rowwise(.x)) %>%
  imap(~ mutate(.x, cos_sim = cos_sim_across_days[[.y]][as.character(Days_pi), as.character(second_timepoint)]))

cos_sim_subsequent_timepoints %>%
  bind_rows() %>%
  ggplot(aes(x = second_timepoint)) +
  geom_area(aes(y = mean_ul_urine), stat = 'identity', fill = 'gray85') +
  geom_line(aes(y = cos_sim * scale_cos_sim_to_mean_ul_urine(urine_bc_levels),
                color = as_factor(animal)),
            linewidth = 0.8,
            show.legend = FALSE) +
  scale_x_continuous(name = 'Days Post-Injection') +
  scale_y_continuous(name = 'Genomes Per uL',
                     sec.axis = sec_axis(~ ./scale_cos_sim_to_mean_ul_urine(urine_bc_levels, .y), name = 'Cosine Similarity')) +
  scale_color_manual(values = animal_colors) +
  facet_wrap(vars(animal), ncol = 1) +
  #labs(title = 'Cosine Similarity Between Subsequent Time Points vs Total Shed Virus') +
  theme_minimal()

ggsave('fig3_cos_sim_subsequent_days_all_line.pdf',
       height = 8,
       width = 7,
       path = '../plots/urine_tissue_manuscript2')

## -- Urine Ridge Plots ---------------------------------------------------

winner_barcode_levels = urine_bc_levels %>%
  group_by(animal, barcode) %>%      # Get highest level for each barcode across days
  mutate(max_level = max(bc_level)) %>%
  ungroup() %>%
  nest_by(animal, max_level) %>%     
  group_by(animal) %>%               
  slice_max(max_level, n = 10) %>%  
  ungroup() %>% 
  unnest(cols = data)              

plots_top30_ridge <- urine_bc_levels %>% group_by(animal, barcode) %>%
  top_barcodes_by_max(30) %>%
  split(.$animal) %>%
  imap(~ ggplot(.x, aes(x = Days_pi, y = barcode, height = bc_level)) +
         geom_ridgeline(fill = animal_colors[[.y]],
                        color = 'grey90',
                        scale = 0.0005) +
         scale_x_continuous(expand = c(0.01, 0)) +
         scale_y_discrete(expand = c(0.01, 0)) +
         labs(x = 'Days post-injection',
              y = 'Barcode',
              #title = 'Top 30 Barcodes Ridge Plot',
              subtitle = paste(.y)) +
         theme_ridges())

top30_ridge_plot_combined <- grid.arrange(
  plots_top30_ridge[["FL"]], plots_top30_ridge[["FR"]],
  plots_top30_ridge[["ML"]], plots_top30_ridge[["MR"]],
  ncol = 4, nrow = 1#,  # Two columns, two rows
  #top = "Top 30 Barcodes Ridge Plot"  # Add a title
)

ggsave("fig6_combined_ridge_plots.pdf", 
       plot = top30_ridge_plot_combined , 
       path = "../plots/urine_tissue_manuscript2", 
       height = 12, 
       width = 18)


## -- GC Plot -------------------------------------------------------------

mean_gc_content_plot <- tibble(animal = c('FL', 'FR', 'ML', 'MR'),        # Adding stock barcodes as timepoint 0
                               Days_pi = c(0, 0, 0, 0)) %>%
  mutate(stock = map(.x = animal, ~ cutoff_99pct_stock_barcodes)) %>%
  unnest(., stock) %>%
  select(animal, Days_pi, barcode, pct) %>%
  bind_rows(urine_bc_levels %>%
              select(animal, Days_pi, barcode, frac_weighted_count) %>%
              rename(pct = frac_weighted_count)) %>%
  mutate(GC_pct = str_count(barcode, '[GC]') / str_length(barcode)) %>%
  mutate(wt_GC_pct = GC_pct * pct) %>%          # multiply %GC of barcode by %barcodes in sample
  group_by(animal, Days_pi) %>%
  summarize(timepoint_GC_pct = sum(wt_GC_pct)) %>%
  ggplot(aes(x = Days_pi, y = timepoint_GC_pct)) +
  geom_line(aes(color = as_factor(animal)),
            linewidth = 0.8,
            show.legend = FALSE) +
  facet_wrap(vars(animal), ncol = 1) +
  labs(title = 'A',
       x = 'Days Post-Injection',
       y = 'GC Content') +
  scale_color_manual(values = animal_colors) +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent_format()) +
  theme_minimal()+
  theme(plot.title = element_text(face = "bold"))


GC_top10barcodes_boxplot <- urine_bc_levels %>%
  top_barcodes_by_max(10) %>%
  distinct(animal, barcode) %>%
  mutate(GC_pct = str_count(barcode, '[GC]') / str_length(barcode)) %>%
  bind_rows(., gc_pct_all_barcodes) %>%
  mutate(animal = fct_relevel(animal, c('FL', 'FR', 'ML', 'MR', 'All Barcodes'))) %>%
  ggplot(aes(animal, GC_pct)) +
  geom_boxplot(aes(color = animal), linewidth = 0.8, show.legend = FALSE) +
  labs(title = 'B',
       x = 'Mouse',
       y = 'GC Content') +
  scale_color_manual(values = c(animal_colors, `All Barcodes` = 'gray70')) +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent_format()) +
  theme_minimal()+
  theme(plot.title = element_text(face = "bold"))

gc_combined_plots <- grid.arrange(mean_gc_content_plot , GC_top10barcodes_boxplot, ncol = 2)
gc_combined_plots 
ggsave('fig3s_GC_content_urine.pdf',
       plot = gc_combined_plots,
       path = '../plots/urine_tissue_manuscript2')


urine_bc_levels %>%
  top_barcodes_by_max(10) %>%
  distinct(animal, barcode) %>%
  mutate(length = str_length(barcode)) %>%
  select(animal, barcode, length) %>%
  ## count(animal, length) %>%
  ggplot(aes(as_factor(length))) +
  geom_bar(aes(fill = animal),
           ## binwidth = 1,
           ## breaks = seq(4, 12),
           show.legend = FALSE) +
  facet_wrap(vars(animal), ncol = 1) +
  labs(title = '',
       x = 'Barcode Length',
       y = 'Count') +
  #ylim(0, 10) +
  scale_y_discrete(limits = seq(0, 10, 2)) +
  scale_fill_manual(values = animal_colors) +
  theme_minimal()

ggsave('fig4s_length_urine_winners_histogram.pdf',
       path = '../plots/urine_tissue_manuscript2',
       height = 8)



## -- Rank Urine Winners In Stock -----------------------------------------

urine_bc_levels %>%
  top_barcodes_by_max(10) %>%
  filter(bc_level == max_level) %>%
  split(.$animal) %>%
  imap(~ mutate(.x, urine_rank = desc(max_level) %>% min_rank() %>% as_factor())) %>%
  bind_rows() %>% 
  left_join(cutoff_99pct_stock_barcodes %>%
              mutate(stock_rank = desc(count) %>% min_rank())) %>%
  select(animal, max_level, barcode, urine_rank, count, stock_rank) %>%
  ggplot(aes(x = fct_relevel(urine_rank, rev), y = stock_rank)) +
  geom_point(aes(color = as_factor(animal)), size=3, show.legend=FALSE) +
  facet_wrap(vars(animal)) +
  labs(title = '',
       x = "Rank Of Winner In Urine",
       y = "Rank In Stock") +
  scale_y_continuous(limits = c(1, dim(cutoff_99pct_stock_barcodes)[[1]])) +
  scale_x_discrete(breaks = seq(1, 10)) +
  scale_color_manual(values=animal_colors) +
  coord_flip() +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(linetype = 'dotted',
                                          color = 'gray70'),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(2, 'lines'),
        strip.background = element_rect(fill = 'gray90',
                                        color = 'gray90')
  )

ggsave('fig8_rank_urine_winners_in_stock_dot.pdf',
       path = '../plots/urine_tissue_manuscript2',
       height=8,
       width=8)

### Calculate median or mean rank of all 40 top barcodes (10 per animal) in stock
urine_bc_levels %>%
  top_barcodes_by_max(10) %>%
  distinct(barcode) %>%
  left_join(cutoff_99pct_stock_barcodes %>%
              mutate(stock_rank = min_rank(desc(count)))) %>%
  pull(stock_rank) %>%
  median()

### Tables of median rank of urine winners in stock
urine_bc_levels %>%
  top_barcodes_by_max(10) %>%
  distinct(animal, barcode) %>%
  left_join(cutoff_99pct_stock_barcodes %>%
              mutate(stock_rank = min_rank(desc(count)))) %>%
  split(.$animal)

## -- Diversity Of Genomes In Urine ---------------------------------------

# make it relative to Day1

urine_complexity_table = urine_bc_levels %>%
  group_by(animal, Days_pi) %>%
  arrange(desc(frac_weighted_count), .by_group = TRUE) %>%
  mutate(cum_frac = cumsum(frac_weighted_count)) %>%
  summarize(entropy = entropy(frac_weighted_count),
            num_to_50pct = sum(cum_frac < 0.5) + 1,
            num_to_75pct = sum(cum_frac < 0.75) + 1,
            .groups = 'drop') %>%
  group_by(animal) %>%
  mutate(entropy_percent = 100 * entropy / first(entropy[Days_pi == 1]),
         num_to_50pct_percent = 100 * num_to_50pct / first(num_to_50pct[Days_pi == 1]),
         num_to_75pct_log = log(num_to_75pct),
         num_to_75pct_percent = 100 * num_to_75pct_log/ first(num_to_75pct_log[Days_pi == 1])) %>%
  ungroup()

urine_complexity_table %>%
  ggplot(aes(x = Days_pi)) +
  geom_line(aes(y = num_to_75pct_percent, color = 'Number To 75% (Log-Scaled)')) +
  geom_line(aes(y = entropy_percent, color = 'Entropy'), linewidth = 1) +
  facet_wrap(vars(animal), ncol = 1, scales = 'free') +
  labs(title = bquote(''),#subtitle = 'Percentage Scale',
       x = 'Days Post-Injection',
       y = 'Percentage Relative to Day 1') +
  scale_color_manual(name = 'Measure', values = complexity_colors[c('Entropy', 'Number To 75% (Log-Scaled)')]) +
  theme_minimal() +
  theme(axis.text.y = element_text(),
        axis.title.y = element_text())

ggsave('Fig5_real_diversity_entropy_num75_log_real.pdf',  ### modify needed
       path='../plots/urine_tissue_manuscript2',
       width = 8,
       height = 6)


## -- Urine Diversity Donuts ----------------------------------------------


# Function to generate the diversity donut plot
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

animals <- c("FL", "FR", "ML", "MR")
donut_plot <- map(animals, ~ {
  bc_table <- urine_bc_levels %>% filter(animal == .x)
  barcode_highlights <- winner_barcodes_by_animal[[.x]]
  
  plot_diversity_donut_w_colored_barcodes(bc_table, barcode_highlights, Days_pi) +
    labs(title = .x)
})

#gridExtra::grid.arrange(grobs = donut_plot, ncol = 2)

combined_donut_plot <- gridExtra::grid.arrange(grobs = donut_plot, ncol = 2)


combined_donut_plot <- gridExtra::grid.arrange(
  grobs = donut_plot, 
  ncol = 2,
  widths = c(1, 1.5),  # Make the second column (plot 4) bigger
  heights = c(1, 1)    # Adjust the row heights to equal
)


ggsave(
  filename = "fig4_diversity_changes_over_time.pdf", 
  plot = combined_donut_plot, 
  path = "../plots/urine_tissue_manuscript2", 
  width = 8, height = 8
)

ggsave(
  filename = "fig4_diversity_changes_over_time_check.pdf", 
  plot = combined_donut_plot, 
  path = "../plots/urine_tissue_manuscript2", 
  width = 8, height = 8
)
#width 8 and height 11 worked


#### tissue stuff



tissue_pcr <- read_excel('organ_pcr.xlsx')

write_tsv(tissue_pcr, 'organ_pcr.tsv')

tissue_pcr = read_excel('organ_pcr.xlsx') %>%
  rename(sample = `Sample Name`, 
         animal = Animal, 
         organ = Organ,
         quantity_mean_ug_dna = `Quantity Mean/µg DNA`,
         quantity_mean_ul_blood = `Quantity Mean/µl blood or virus stock`) %>%
  select(sample, animal, organ, quantity_mean_ug_dna, quantity_mean_ul_blood) %>%
  mutate(blood = str_detect(organ, '.*blood.*')) %>%
  mutate(organ = recode(organ,
                        Gonades = 'Gonads',
                        `Plasma*` = 'Plasma',
                        Testicules = 'Testicles',
                        `Whole blood*` = 'Whole blood'))

colnames(tissue_pcr)


tissue_pcr2 <- tissue_pcr %>%
  mutate(sample2 = paste0(animal, "_", str_replace_all(organ, " ", "_")))

head(tissue_pcr2 %>% select(sample, sample2))

tissue_pcr2<- tissue_pcr2 %>%
  rename(
    sample_backup = sample,  
    sample = sample2       
  )



tissue_bc_levels = tissue_to_stock_barcodes %>%
  group_by(sample) %>%
  mutate(frac_weighted_count = weighted_count / sum(weighted_count)) %>%
  ungroup() %>%
  left_join(tissue_pcr2,by = 'sample') %>%
  mutate(bc_level = if_else(blood,
                            frac_weighted_count * as.numeric(quantity_mean_ul_blood) * 1000,
                            frac_weighted_count * as.numeric(quantity_mean_ug_dna))) %>%
  select(sample, animal, organ, stock_bc, bc_level, frac_weighted_count, quantity_mean_ug_dna, quantity_mean_ul_blood, blood) %>%
  rename(barcode = stock_bc) %>%
  filter(sample != 'MR_Plasma') 

# Rank of tissue winners in stock

tissue_bc_levels %>%
  top_barcodes_by_max(10) %>%
  filter(bc_level == max_level) %>%
  split(.$animal) %>%
  imap(~ mutate(.x, tissue_rank = desc(max_level) %>% min_rank() %>% as_factor())) %>%
  bind_rows() %>% 
  left_join(cutoff_99pct_stock_barcodes %>%
              mutate(stock_rank = desc(count) %>% min_rank())) %>%
  select(animal, max_level, barcode, tissue_rank, count, stock_rank) %>%
  ggplot(aes(x = fct_relevel(tissue_rank, rev), y = stock_rank)) + 
  geom_point(aes(color = as_factor(animal)), size=3, show.legend=FALSE) +
  facet_wrap(vars(animal)) +
  labs(title = '',
       x = "Rank Of Winner In Tissue",
       y = "Rank In Stock") +
  scale_y_continuous(limits = c(1, dim(cutoff_99pct_stock_barcodes)[[1]])) +
  scale_x_discrete(breaks = seq(1, 10)) +
  scale_color_manual(values=animal_colors) +
  coord_flip() +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(linetype = 'dotted',
                                          color = 'gray70'),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(2, 'lines'),
        strip.background = element_rect(fill = 'gray90',
                                        color = 'gray90')
  ) 
ggsave('fig9_rank_tissue_winners_in_stock_dot.pdf',
       path = '../plots/urine_tissue_manuscript2',
       height = 8,
       width = 8)

# rank of organ winners in last day urine

tissue_10_from_each_organ_by_animal = tissue_bc_levels %>%
  split(.$animal) %>%
  imap(~ group_by(.x, organ) %>%
         slice_max(bc_level, n = 10) %>%
         mutate(organ_rank = desc(bc_level) %>% min_rank() %>% as_factor()) %>%
         select(-quantity_mean_ul_blood, -frac_weighted_count)) 

urine_last_day_rank_by_animal = urine_bc_levels %>%
  filter(bc_level > 0) %>%
  group_by(animal) %>%
  mutate(last_day = max(Days_pi)) %>%
  filter(Days_pi == last_day) %>% 
  ungroup() %>%
  split(.$animal) %>%
  imap(~ mutate(.x, last_day_rank = desc(bc_level) %>% min_rank())) %>%
  imap(~ select(.x, barcode, last_day_rank))

combined_data <- map2_dfr(tissue_10_from_each_organ_by_animal, urine_last_day_rank_by_animal, left_join)

combined_plot <- ggplot(combined_data, aes(x = last_day_rank, y = fct_relevel(organ_rank, rev))) +
  geom_point(aes(color = animal), size = 2, show.legend = FALSE) + 
  facet_grid(animal ~ organ) +  
  scale_color_manual(values = animal_colors) + 
  scale_x_continuous(limits = c(1, 350)) +
  labs(title = '',
       x = "Rank of Winners in Last Day Urine",
       y = "Rank Of Winner In Organ") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(linetype = 'dotted',
                                          color = 'gray70'),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(2, 'lines'),
        strip.background = element_rect(fill = 'gray90',
                                        color = 'gray90')
  )

combined_plot

ggsave("fig7s_organ_winners_in_last_day_urine_by_organ_dot.pdf", combined_plot,path = '../plots/urine_tissue_manuscript2',
       width = 20, height = 12)


## -- Urine Winners In Stock And Tissue -----------------------------------

## Get median rank of urine in stock
urine_rank_in_stock_table = urine_bc_levels %>%
  split(.$animal) %>%
  imap(~ top_barcodes_by_max(.x, 10)) %>%
  imap(~ arrange(.x, bc_level)) %>%
  imap(~ left_join(.x, mutate(cutoff_99pct_stock_barcodes, stock_ranking = desc(count) %>% min_rank()))) %>%
  bind_rows() %>%
  group_by(animal) %>%
  summarize(median_stock_rank = median(stock_ranking)) %>%
  rename(Animal = animal, Urine = median_stock_rank)

## Get median rank of tissue in stock
tissue_rank_in_stock_table = tissue_bc_levels %>%
  split(.$animal) %>%
  imap(~ top_barcodes_by_max(.x, 10)) %>%
  imap(~ arrange(.x, bc_level)) %>%
  imap(~ left_join(.x, mutate(cutoff_99pct_stock_barcodes, stock_ranking = desc(count) %>% min_rank()))) %>%
  bind_rows() %>%
  group_by(animal) %>%
  summarize(median_stock_rank = median(stock_ranking)) %>%
  rename(Animal = animal, Tissue = median_stock_rank)

left_join(urine_rank_in_stock_table, tissue_rank_in_stock_table) %>%
  kbl(caption = 'Median Rank Of Urine And Tissue Winners In Stock') %>%
  kable_styling(bootstrap_options = 'striped', full_width = FALSE) %>%
  column_spec(2, width = '8em') %>%
  column_spec(3, width = '8em') %>%
  save_kable('../plots/urine_tissue_manuscript2/urine_and_tissue_winner_median_stock_ranks_table.html')

tissue_bc_levels %>%
  slice_max(n = 10, order_by = bc_level, by = c(animal, organ)) %>%
  group_by(animal) %>%
  distinct(barcode) %>%
  mutate(GC_pct = str_count(barcode, '[GC]') / str_length(barcode)) %>%
  ungroup() %>%
  bind_rows(., gc_pct_all_barcodes) %>%
  mutate(animal = fct_relevel(animal, c('FL', 'FR', 'ML', 'MR', 'All Barcodes'))) %>%
  ggplot(aes(animal, GC_pct)) +
  geom_boxplot(aes(color = animal), linewidth = 0.8, show.legend = FALSE) +
  labs(title = '',
       x = 'Mouse',
       y = 'GC Content') +
  scale_color_manual(values = c(animal_colors, `All Barcodes` = 'gray70')) +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent_format()) +
  theme_minimal()

ggsave('fig6s_GC_tissue_winners_boxplot.pdf',
       path = '../plots/urine_tissue_manuscript2',
       width = 8,
       height = 6)




### Spearman correlation within tissues per animal

generate_correlation_plot <- function(data, animal_id) {
  animal_data <- data %>%
    filter(animal == animal_id) %>%
    select(organ, barcode, frac_weighted_count) %>%
    pivot_wider(names_from = organ, values_from = frac_weighted_count, values_fill = 0)
  
  organ_matrix <- animal_data %>% select(-barcode) %>% as.matrix()
  
  cor_matrix <- cor(organ_matrix, method = "spearman")
  
  cor_long <- as.data.frame(as.table(cor_matrix))
  colnames(cor_long) <- c("Organ1", "Organ2", "Spearman_Correlation")
  
  cor_plot <- ggplot(cor_long, aes(x = Organ1, y = Organ2, fill = Spearman_Correlation)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1)) +
    geom_text(aes(label = sprintf("%.2f", Spearman_Correlation)), color = "black", size = 3) +
    labs(
      title = paste("Spearman Correlation of Barcode Abundance"),
      subtitle = paste("Animal:", animal_id),
      x = "Tissue 1",
      y = "Tissue 2",
      fill = "Spearman\nCorrelation"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
  
  return(cor_plot)
}

animals <- c("FL", "FR", "MR", "ML")
correlation_plots <- lapply(animals, function(animal) {
  generate_correlation_plot(tissue_bc_levels, animal)
})

library(gridExtra)
gridExtra::grid.arrange(grobs = correlation_plots, ncol = 2)

ggsave("../plots/tissue_correlation_plots.pdf", 
       gridExtra::marrangeGrob(grobs = correlation_plots, ncol = 2, nrow = 2),
       width = 14, height = 10)


#try2

# Load this for `textGrob`

compute_spearman_correlation <- function(data, animal_id) {
  animal_data <- data %>%
    filter(animal == animal_id) %>%
    select(organ, barcode, frac_weighted_count) %>%
    pivot_wider(names_from = organ, values_from = frac_weighted_count, values_fill = 0) %>%
    column_to_rownames(var = "barcode")  
  
  cor_matrix <- cor(as.matrix(animal_data), method = "spearman", use = "pairwise.complete.obs")
  
  colnames(cor_matrix) <- gsub(" ", "\n", colnames(cor_matrix))
  rownames(cor_matrix) <- gsub(" ", "\n", rownames(cor_matrix))
  
  return(cor_matrix)
}

generate_correlation_plot <- function(cor_matrix, animal_id, animal_color) {
  corrplot.mixed(
    cor_matrix,
    upper = 'circle',                  
    lower = 'number',                 
    tl.col = "black",
    tl.font = 2,                       
    tl.cex = .6,                     
    tl.font = 2,
    #tl.pos = "lt",                   
    tl.srt = 45,  
    number.cex = 0.5,                 
    lower.col = "black",              
    is.corr = TRUE,                    
    upper.col = colorRampPalette(c("gray", "white", animal_color))(400), 
    col.lim = c(-1, 1), #
    #title = paste(animal_id)             
    # Add the title closer to the plot
    #title(main = paste("Animal:", animal_id), line = -1, cex.main = 1.2)
  )
}

generate_correlation_plot 


pdf("../plots/urine_tissue_manuscript2/Fig10_tissue_correlation_by_animal.pdf", onefile = TRUE, height = 9, width = 10)


par(mfrow = c(2, 2), mar = c(2, 2, 4, 2))  

for (animal_id in names(animal_colors)) {
  cor_matrix <- compute_spearman_correlation(tissue_bc_levels, animal_id)
  generate_correlation_plot(cor_matrix, animal_id, animal_colors[animal_id])
  title(main = paste(animal_id), line = 3, cex.main = 1.2)
}

par(mfrow = c(1, 1))

dev.off()



# Fig s5 
# total amounts of muPyV DNA in organs

#write.table(tissue_pcr2, file = "../plots/urine_tissue_manuscript2/tissue_pcr2.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

tissue_pcr2_filtered <- tissue_pcr2 %>%
  filter(!organ %in% c("PLMR", "Lung tumor")) %>%
  mutate(across(starts_with("quantity"), ~ ifelse(is.na(.), "N.D.", .))) %>%
  mutate(genome_equivalent = ifelse(blood == FALSE, quantity_mean_ug_dna, quantity_mean_ul_blood)) %>%
  mutate(genome_equivalent = as.numeric(genome_equivalent))

tissue_pcr2_filtered <- tissue_pcr2_filtered %>%
  mutate(organ = ifelse(organ == "Gonads", "Ovaries", organ))

fig5s_dna.amount.in.organ <- ggplot(tissue_pcr2_filtered %>%
                                      filter(!is.na(genome_equivalent) & genome_equivalent > 0), 
                                    aes(x = organ, y = log10(genome_equivalent), fill = animal)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = animal_colors) + 
  labs(
    title = "",
    x = "Organ",
    y = expression("log"[10] * "(Genome Equivalent)"),
    fill = "Mouse"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    panel.grid.major.x = element_blank(),  # Remove vertical grid lines
    panel.grid.minor.x = element_blank(),
    legend.position = "right"  # Move legend to the right
  ) +
  scale_y_continuous(
    breaks = c(2, 3, 4, 5, 6, 7, 8),  # Custom breaks for log10 scale
    labels = scales::comma  # Format labels with commas
  )

ggsave("fig5s_genome_equivalents_in_organ.pdf", 
       plot = fig5s_dna.amount.in.organ, 
       path= "../plots/urine_tissue_manuscript2/",
       width = 10, height = 8)




tissue_bc_levels %>% filter(bc_level>0) %>%
  group_by(animal) %>%
  summarise(unique_barcodes = n_distinct(barcode))


tissue_bc_levels %>% filter(bc_level>0) %>%
  +   group_by(animal) %>%
  +   summarise(unique_barcodes = n_distinct(barcode))


##########################
#fisher
filtered_urine_bc_levels <- urine_bc_levels %>% filter(bc_level >0)
filtered_tissue_bc_levels <- tissue_bc_levels %>% filter(bc_level >0)

analyze_top5_overlap <- function(urine_data, tissue_data, animal_id) {
  # get top 5% barcodes in urine
  top_urine <- urine_data %>%
    filter(animal == animal_id) %>%
    group_by(barcode) %>%
    summarize(total_urine = sum(bc_level)) %>%
    slice_max(total_urine, prop = 0.05)
  
  #get top 5% barcodes in tissues
  top_tissue <- tissue_data %>%
    filter(animal == animal_id) %>%
    group_by(barcode) %>%
    summarize(total_tissue = sum(bc_level)) %>%
    slice_max(total_tissue, prop = 0.05)
  
  # contingency table banabo
  contingency <- tibble(
    category = c("Top Urine & Top Tissue", "Top Urine Only",
                 "Top Tissue Only", "Neither"),
    count = c(
      sum(top_urine$barcode %in% top_tissue$barcode),
      nrow(top_urine) - sum(top_urine$barcode %in% top_tissue$barcode),
      nrow(top_tissue) - sum(top_urine$barcode %in% top_tissue$barcode),
      nrow(clustered_combined_stock_barcodes) - 
        (sum(top_urine$barcode %in% top_tissue$barcode) + 
           (nrow(top_urine) - sum(top_urine$barcode %in% top_tissue$barcode)) + 
           (nrow(top_tissue) - sum(top_urine$barcode %in% top_tissue$barcode)))
    )
  )  
  
  # Fisher's exact test
  mat <- matrix(contingency$count[1:4], nrow = 2)
  ftest <- fisher.test(mat)
  
  return(list(
    animal = animal_id,
    p_value = ftest$p.value,
    contingency_table = contingency
  ))
}


#library(knitr)  # For kable()
#library(kableExtra) 

animals <- c("FL", "FR", "ML", "MR")
#results <- map(animals, ~ analyze_top5_overlap(urine_bc_levels, tissue_bc_levels, .x))
results <- map(animals, ~ analyze_top5_overlap(filtered_urine_bc_levels, filtered_tissue_bc_levels, .x))

# Print results
walk(results, function(res) {
  cat("\n", res$animal, ":\n")
  cat("--------------------------------\n")
  print(kable(res$contingency_table))
  cat("Fisher's exact test p-value:", format.pval(res$p_value, digits = 3), "\n")
})

#plot_fisher


results_df <- map_dfr(results, ~ {
  .x$contingency_table %>%
    mutate(animal = .x$animal,
           p_value = .x$p_value)
})

# overlap_plot <- results_df %>%
#   group_by(animal) %>%
#   mutate(total = sum(count),
#          prop = count/total) %>%
#   ggplot(aes(x = animal, y = prop, fill = category)) +
#   geom_col() +
#   geom_text(aes(label = scales::percent(prop, accuracy = 0.1)), 
#             position = position_stack(vjust = 0.5),
#             size = 3) +
#   scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728")) +
#   labs(title = "Overlap Between Top 5% Barcodes in Urine and Tissues",
#        x = "Mouse",
#        y = "Proportion",
#        fill = "Category") +
#   theme_minimal() +
#   theme(legend.position = "bottom")

pvalue_plot <- results_df %>%
  distinct(animal, p_value) %>%
  mutate(signif = case_when(
    p_value < 0.00001 ~ "*****",
    p_value < 0.0001 ~ "****",
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ "NS"
  )) %>%
  ggplot(aes(x = animal, y = -log10(p_value), fill = animal)) +
  geom_col(aes(fill = animal), show.legend = FALSE) +
  geom_text(aes(label = signif), vjust = -0.5, size = 5) +
  scale_fill_manual(values = animal_colors) +
  labs(title = "",
       #"Statistical Significance of Overlap Between Top 5% Barcodes in Urine and Tissues",
       x = "Mouse",
       y = "-log10(p-value)") +
  theme_minimal()
pvalue_plot
ggsave("../plots/urine_tissue_manuscript2/urine_tissue_overlap_pval_plot.pdf", 
       pvalue_plot,
       width = 8, height = 8)
combined_plot <- overlap_plot / pvalue_plot + 
  plot_layout(heights = c(2, 1))
combined_plot

ggsave("../plots/urine_tissue_manuscript2/urine_tissue_overlap_analysis.pdf", 
       combined_plot,
       width = 10, height = 16)


