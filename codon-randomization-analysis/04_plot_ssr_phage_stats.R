library(tidyverse)
library(ggplot2)

#Include unknown/hypothetical category
use_unknown = T

my_theme= theme_bw() + theme(
  panel.grid = element_blank()
)

my_theme_no_x_labels = my_theme + theme(
  axis.text.x = element_blank(),
  axis.title.x = element_blank()
)

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

metadata = metadata %>% filter(GenBank %in% matched_names)
metadata$Genus[is.na(metadata$Genus)] = "Unknown"

metadata %>% group_by(Genus) %>% summarize(n=n())

### Load phage genome sizes and GC content
phage_stats = data.frame()
for (b in matched_names) {
  this_stats = read_csv(file.path(main_output_path, "output_gene_gc_content",paste0(b, ".csv")))
  phage_stats=bind_rows(phage_stats, this_stats)
}
phage_stats$seq_id = sub("\\..*", "", phage_stats$seq_id)
phage_stats = phage_stats %>% 
  mutate(
    GC=gc_percent*total_bases/100, 
    AT=(100-gc_percent)*total_bases/100
    )
phage_stats_longer = phage_stats %>% pivot_longer(cols=c(GC,AT), names_to="base_type", values_to="bases")

### Load phages actual SSR
ssrs = data.frame()
for (b in matched_names) {
  actual_ssrs = read_csv(file.path(main_output_path, "output_gene_ssrs_actual",paste0(b, ".csv")))
  ssrs=bind_rows(ssrs, actual_ssrs)
}
ssrs$seq_id = sub("\\..*", "", ssrs$seq_id)
ssrs_per_phage = ssrs %>% filter(length>=6) %>% group_by(seq_id) %>% summarize(n=sum(count))
ssrs_per_phage_by_base_length = ssrs %>% filter(length>=6) %>% group_by(seq_id, base_pair, length) %>% summarize(n=sum(count))
ssrs_per_phage = ssrs_per_phage %>% left_join(metadata %>% select(GenBank, Genus) %>% rename(seq_id=GenBank), by="seq_id")
actual_ssrs_by_genus = ssrs_per_phage %>% group_by(Genus) %>% summarize(n=sum(n))

genera  = unique(ssrs_per_phage$Genus)

### Load results of randomization tests

randomized_ssrs = data.frame()
for (b in matched_names) {
  this_ssrs = read_csv(file.path(main_output_path, "output_gene_ssrs_codon_resampled_by_genome",paste0(b, ".csv")))
  randomized_ssrs=bind_rows(randomized_ssrs, this_ssrs)
}
randomized_ssrs$seq_id = sub("\\..*", "", randomized_ssrs$seq_id)
randomized_ssrs_per_phage = randomized_ssrs %>% filter(length>=6) %>% group_by(seq_id, replicate) %>% summarize(n=sum(count))

### Load category information

categories = c(
  "assembly",
  #  "immune",
  "infection",
  "integration",
  "lysis",
  "packaging",
  "regulation",
  "replication"
  #  "tRNA_related"
)

if (use_unknown) {
  categories = c(categories, "unknown")
}

okabe_ito_categories =setNames(
  c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00","#CC79A7", "#EEEEEE"),
  categories
  #  "grey"         = "#999999"
)

phage_proteins = read_csv("phage_protein_accessions.csv")
phage_proteins$seq_id = sub("\\..*", "", phage_proteins$seq_id)
phage_proteins = phage_proteins %>% filter(seq_id %in% matched_names)

proteins_with_categories = read_csv("proteins_with_categories.csv")

phage_proteins_total_gene_length = phage_proteins %>% left_join(proteins_with_categories, by="protein_id") %>% filter(!is.na(nt_length)) %>% group_by(seq_id) %>% summarize(total_gene_length = sum(nt_length))

# only show 6, 7, 8+ as categories
ssrs_per_phage_by_base_length_capped = ssrs_per_phage_by_base_length

ssrs_per_phage_by_base_length_capped$length[ssrs_per_phage_by_base_length_capped$length>8] = 8

ssrs_per_phage_by_base_length_capped = ssrs_per_phage_by_base_length_capped %>% 
  mutate(base_pair_length = paste0(base_pair, "-", length)) %>% group_by(seq_id,base_pair_length) %>% summarize(n=sum(n)) %>% ungroup()

ssrs_per_phage_by_base_length_capped = ssrs_per_phage_by_base_length_capped %>% 
  complete(seq_id=unique(phage_stats$seq_id), base_pair_length, fill=list(n=0))

ssrs_per_phage_by_base_length_capped = ssrs_per_phage_by_base_length_capped %>% left_join(phage_proteins_total_gene_length, by="seq_id") 

ssrs_per_phage_by_base_length_capped = ssrs_per_phage_by_base_length_capped %>% mutate(ssrs_per_kb = n/total_gene_length*1000)

phage_proteins = phage_proteins %>% left_join(proteins_with_categories, by="protein_id") %>% filter(!is.na(nt_length))

proteins_with_categories = proteins_with_categories %>%
separate_rows(protein_classification, sep = ";") %>%
  mutate(protein_classification = str_trim(protein_classification)) %>% 
  filter(!is.na(protein_classification)) %>%
  filter(protein_classification != "")

#if a protein has n classifications, only count each as 1/n
protein_category_weights = proteins_with_categories %>% group_by(protein_id) %>% summarize(category_weight =1/n())

phage_proteins = phage_proteins %>% select(seq_id, protein_id) %>% left_join(proteins_with_categories, by="protein_id") %>% left_join(protein_category_weights, by="protein_id")

phage_proteins$protein_classification <- factor(
  phage_proteins$protein_classification,
  levels = rev(categories)
)

#### Ordering seq_id factor

phage_stats = phage_stats %>% arrange(-total_bases)
seq_id_ordered_levels = phage_stats$seq_id

### Panel: Actual relative to resampled

## Note pseudocount of 0.5 added to actual and resamples
resampled_relative = randomized_ssrs_per_phage %>% left_join(ssrs_per_phage %>% rename(actual_n=n), by="seq_id") %>% replace_na(list(actual_n = 0)) %>% mutate(relative_n = (actual_n+0.5)/(n+0.5))

resampled_relative_quantiles <- resampled_relative %>% select(seq_id, relative_n) %>% group_by(seq_id) %>%
  summarise(
    q025 = quantile(relative_n, 0.025),
    q975 = quantile(relative_n, 0.975)
  )

metadata = metadata %>% rename(seq_id=GenBank)
metadata$seq_id=factor(metadata$seq_id, levels=seq_id_ordered_levels)
metadata = metadata %>% arrange(seq_id)

resampled_relative$seq_id=factor(resampled_relative$seq_id, levels=seq_id_ordered_levels)

## version that shows violin plot
resampled_relative_panel = ggplot(resampled_relative, aes(x=seq_id,y=relative_n)) + 
  stat_summary(fun = mean, geom="bar") + geom_violin() + theme_bw()  + geom_hline(yintercept=1) + my_theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + scale_x_discrete(labels = metadata$ID) 

## version that shows 95% CL
resampled_relative_quantiles

resampled_relative_panel = ggplot(resampled_relative, aes(x=seq_id,y=relative_n)) + 
  stat_summary(fun = mean, geom="bar") + geom_errorbar(aes(ymin=q025, ymax=q975, y=NA), data=resampled_relative_quantiles, width=0) + theme_bw()  + geom_hline(yintercept=1) + my_theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + scale_x_discrete(labels = metadata$ID) 

### Panel: genome size

phage_stats_longer$seq_id=factor(phage_stats_longer$seq_id, levels=seq_id_ordered_levels)

base_pair_palette = list(
  "GC" = "#0000EE",
  "AT" = "#EE0000"
)

genome_size_panel = 
  
  ggplot(phage_stats_longer, aes(x=seq_id,y=bases, fill=base_type)) + 
  geom_col() + theme_bw()  +
  scale_fill_manual(values=base_pair_palette) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_continuous(limits=c(0,180000), expand=c(0,0)) + my_theme_no_x_labels

### Panel gene categories
  
phage_proteins$seq_id=factor(phage_proteins$seq_id, levels=seq_id_ordered_levels)

phage_proteins$category_weight_by_length = phage_proteins$category_weight * phage_proteins$nt_length

protein_category_panel =
ggplot(phage_proteins, aes(x=seq_id, fill=protein_classification, weight = category_weight_by_length)) + 
  geom_bar(position = "fill") + theme_bw()  +
  scale_fill_manual(values=okabe_ito_categories) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  scale_y_continuous(limits=c(0,1), expand=c(0,0)) + my_theme_no_x_labels


phage_protein_categories_by_length = phage_proteins %>% group_by(seq_id, protein_classification) %>% summarize(category_weight_by_length = sum(category_weight_by_length))

phage_protein_categories_by_length = phage_protein_categories_by_length %>% pivot_wider(names_from=protein_classification, values_from = category_weight_by_length, values_fill=0)

write_csv(phage_protein_categories_by_length, "processed_data/phage_protein_categories_by_length.csv")


### Panel SSRs per kb

ssrs_per_phage_by_base_length_capped$seq_id=factor(ssrs_per_phage_by_base_length_capped$seq_id, levels=seq_id_ordered_levels)

base_pair_length_palette = list(
  "GC-6" = "#8888FF",
  "GC-7" = "#0000EE",
  "GC-8" = "#000099",
  "AT-6" = "#FF8888",
  "AT-7" = "#EE0000",
  "AT-8" = "#990000"
)

ssrs_per_kb_panel =
  ggplot(ssrs_per_phage_by_base_length_capped, aes(x=seq_id, y=ssrs_per_kb, fill=base_pair_length)) + 
  geom_col() + theme_bw()  +
  scale_fill_manual(values=base_pair_length_palette) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  scale_y_continuous(limits=c(0,2.5), expand=c(0,0)) + my_theme_no_x_labels

### Put the panels together

library(patchwork)
wrap_plots(ssrs_per_kb_panel, genome_size_panel, protein_category_panel, resampled_relative_panel, ncol = 1)

write_csv(resampled_relative %>% group_by(seq_id) %>% summarize(relative_n = mean(relative_n)) %>% mutate(id = metadata$ID) , "processed_data/relative_abundance_ssrs_per_phage.csv")
ggsave("plots/ss_phage_stats.pdf", height=8, width=8)
  
