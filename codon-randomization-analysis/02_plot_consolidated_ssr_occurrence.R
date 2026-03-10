library(tidyverse)
library(ggplot2)

dir.create("plots", recursive = TRUE, showWarnings = FALSE)

min_ssr_length = 1
max_ssr_length = 12

base_pairs = c("GC", "AT")
ssr_lengths = min_ssr_length:max_ssr_length

# Get all the base names so we can construct filenames from them
files <- list.files(
  path = "processed_data/output_gene_gc_content",
  pattern = "\\.csv$",
  full.names = TRUE
)
base_names <- tools::file_path_sans_ext(basename(files))
#base_names <- basename(files)


metadata = read_csv("raw_data/Basel_deduplicated_metadata.csv",  comment = "#")

matched_names = metadata$GenBank[metadata$GenBank %in% base_names]

#Each phage on its own
consolidated_df = data.frame()
for (b in matched_names) {
  
  #b = "MH751506"
  # Load GC content
  gc = read_csv(file.path("processed_data", "output_gene_gc_content", paste0(b, ".csv")))
  gc_fraction = gc$gc_percent[1]*0.01
  at_fraction = 1 - gc_fraction
  total_bases = gc$total_bases
  
  # Calculate expected numbers - must have n of base and then not that base
  expected_gc_ssrs = data.frame(length=min_ssr_length:max_ssr_length, base_pair="GC")
  expected_gc_ssrs$ssrs_per_kb = 2 * 1000 * (gc_fraction/2)^expected_gc_ssrs$length * (1-gc_fraction/2)^2
  expected_at_ssrs = data.frame(length=min_ssr_length:max_ssr_length, base_pair="AT")
  expected_at_ssrs$ssrs_per_kb = 2 * 1000 * (at_fraction/2)^expected_at_ssrs$length * (1-at_fraction/2)^2
  expected_ssrs = expected_gc_ssrs %>% bind_rows(expected_at_ssrs)
  
  #Load Actual SSRs
  actual_ssrs = read_csv(file.path("processed_data", "output_gene_ssrs_actual",paste0(b, ".csv")))
  actual_ssrs_by_length = actual_ssrs %>% group_by(length, base_pair) %>% summarize(ssrs=sum(count))
  actual_ssrs_by_length = actual_ssrs_by_length %>% ungroup() %>% complete(length=ssr_lengths, base_pair=base_pairs, fill = list(ssrs = 0))
  actual_ssrs_by_length = actual_ssrs_by_length %>% mutate(ssrs_per_kb = ssrs / total_bases * 1000)
  
  #graph on log scale
  ggplot(expected_ssrs, aes(x=length,y=ssrs_per_kb, shape=base_pair)) +
    geom_point(color="black") +
    geom_point(data=actual_ssrs_by_length, color="red") +
    scale_y_log10() +
    theme_bw()
  
  # Load genome-randomized counts
  genome_randomized_ssrs = read_csv(file.path("processed_data", "output_gene_ssrs_codon_resampled_by_genome", paste0(b, ".csv")))
  num_genome_replicates = length(unique(genome_randomized_ssrs$replicate))
  genome_randomized_ssrs_by_length = genome_randomized_ssrs %>% group_by(length, base_pair) %>% summarize(ssrs=sum(count)/num_genome_replicates)
  genome_randomized_ssrs_by_length = genome_randomized_ssrs_by_length %>% ungroup() %>% complete(length=ssr_lengths, base_pair=base_pairs, fill = list(ssrs = 0))
  genome_randomized_ssrs_by_length = genome_randomized_ssrs_by_length %>% mutate(ssrs_per_kb = ssrs / total_bases * 1000)
  
  #graph on log scale
  ggplot(expected_ssrs, aes(x=length,y=ssrs_per_kb, shape=base_pair)) +
    geom_point(color="black") + geom_line(color="black") +
    geom_point(data=actual_ssrs_by_length, color="red") + geom_line(data=actual_ssrs_by_length, color="red") +
    geom_point(data=genome_randomized_ssrs_by_length, color="blue") + geom_line(data=genome_randomized_ssrs_by_length, color="blue") +
    scale_y_log10() +
    theme_bw()
  
  #consolidate
  this_consolidated_df = expected_ssrs %>% 
    rename(expected_ssrs_per_kb=ssrs_per_kb) %>%
    left_join(actual_ssrs_by_length %>% select(length, base_pair, ssrs_per_kb) %>% rename(actual_ssrs_per_kb=ssrs_per_kb), by=c("length","base_pair")) %>%
    left_join(genome_randomized_ssrs_by_length %>% select(length, base_pair, ssrs_per_kb) %>% rename(genome_randomized_ssrs_per_kb=ssrs_per_kb), by=c("length","base_pair")) %>% mutate(seq_id=b)
  
  
  consolidated_df = consolidated_df %>% bind_rows(this_consolidated_df)
}

write_csv(consolidated_df, file.path("processed_data", "consolidated_ssrs_by_length.csv"))

# Recalculate overall averages across all phages together
num_phages = length(unique(consolidated_df$seq_id))
summary_consolidated_df = consolidated_df %>% group_by(length, base_pair) %>% 
  summarize(
    expected_ssrs_per_kb = sum(expected_ssrs_per_kb)/num_phages, 
    actual_ssrs_per_kb = sum(actual_ssrs_per_kb)/num_phages,
    gene_randomized_ssrs_per_kb = sum(gene_randomized_ssrs_per_kb)/num_phages,
    genome_randomized_ssrs_per_kb = sum(genome_randomized_ssrs_per_kb)/num_phages
    )

#graph on log scale
ggplot(summary_consolidated_df, aes(x=length,y=expected_ssrs_per_kb, shape=base_pair, color=base_pair)) +
#  geom_point(color="black") + geom_line(color="black") +
  scale_color_manual(values = c("AT"="red", "GC"="blue")) +
  geom_point(aes(y=actual_ssrs_per_kb)) + geom_line(aes(y=actual_ssrs_per_kb)) +
#  geom_point(aes(y=gene_randomized_ssrs_per_kb), color="green") + geom_line(aes(y=gene_randomized_ssrs_per_kb), color="green") +
  #geom_point(aes(y=genome_randomized_ssrs_per_kb)) + 
  geom_line(aes(y=genome_randomized_ssrs_per_kb), linetype=3) +
  scale_y_log10(lim=c(1E-4, 1E2), breaks=c(1E-4, 1E-3, 1E-2, 1E-2, 1E-1, 1E0, 1E1, 1E2)) +
  scale_x_continuous(lim = c(2, 8), breaks=2:8) +
  theme_bw()

ggsave(file.path("plots", "consolidated_ssrs_by_length_dedup.pdf"))

# graph actual and observed on linear scale


