library(tidyverse)

main_output_path = "plots"
#Include unknown/hypothetical category
use_unknown = T

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

proteins_with_categories = read_csv("raw_data/phage_protein_categories.csv")

ssrs = data.frame()
for (b in matched_names) {
  actual_ssrs = read_csv(file.path("processed_data", "output_gene_ssrs_actual",paste0(b, ".csv")))
  ssrs=bind_rows(ssrs, actual_ssrs)
}

ssrs_per_gene = ssrs %>% filter(length>=6) %>% group_by(seq_id, locus_tag) %>% summarize(n=sum(count)) %>% rename(protein_id=locus_tag)

proteins_with_categories_and_ssrs = ssrs_per_gene %>% left_join(proteins_with_categories, by="protein_id")

proteins_with_categories_and_ssrs = proteins_with_categories_and_ssrs %>%
  filter(!is.na(protein_classification))

total_ssrs = sum(proteins_with_categories_and_ssrs$n)

count_ssrs_by_protein_category <- function(df) {
  
  # We need to split apart on protein classification by ;
  ret_df = df %>%
    separate_rows(protein_classification, sep = ";") %>%
    mutate(protein_classification = str_trim(protein_classification)) %>% 
    filter(!is.na(protein_classification)) %>% 
    filter(protein_classification != "") %>%
    group_by(protein_classification) %>%
    summarize(n=sum(n))
}

actual_ssrs_by_protein_category_counts = count_ssrs_by_protein_category(proteins_with_categories_and_ssrs)
sum(actual_ssrs_by_protein_category_counts$n)

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


randomized_ssrs = data.frame()
for (b in matched_names) {
  this_ssrs = read_csv(file.path("processed_data", "output_gene_ssrs_codon_resampled_by_genome",paste0(b, ".csv")))
  randomized_ssrs=bind_rows(randomized_ssrs, this_ssrs)
}
randomized_ssrs_summary = randomized_ssrs %>% filter(length>=6) %>% group_by(replicate,locus_tag) %>% summarize(n=sum(count)) %>% rename(protein_id=locus_tag)

randomized_ssrs_with_categories = randomized_ssrs_summary %>% left_join(proteins_with_categories, by="protein_id")

randomized_ssrs_with_categories = randomized_ssrs_with_categories %>%
  filter(!is.na(protein_classification))

num_randomized_datasets = length(unique(randomized_ssrs$replicate))

resampled_randomizations = data.frame()
for (i in 1:num_randomized_datasets) {
  cat("Resampling: ", i, "\n")
  
  rdf = randomized_ssrs_with_categories %>% filter(replicate==i)
  
  # Merge to set categories in case some weren't sampled
  rcdf = count_ssrs_by_protein_category(rdf) %>% mutate(replicate = i)
  
  resampled_randomizations = resampled_randomizations %>% bind_rows(rcdf)
}


library(ggplot2)


resampled_randomizations_medians = resampled_randomizations %>% group_by() %>% group_by(protein_classification) %>% summarize(median_n=median(n))

resampled_randomization_distributions_log2 = resampled_randomizations %>% left_join(resampled_randomizations_medians) %>% mutate(log2_relative_n = log2(n/median_n))

actual_ssrs_by_protein_category_counts_log2 = actual_ssrs_by_protein_category_counts %>% left_join(resampled_randomizations_medians) %>% mutate(log2_relative_n = log2(n/median_n))

resampled_randomizations_relative_log2 = resampled_randomizations %>% left_join(actual_ssrs_by_protein_category_counts %>% rename(actual_n=n), by="protein_classification") %>% mutate(log2_relative_n = log2(actual_n/n))

resampled_randomizations_relative_log2_median = resampled_randomizations_relative_log2 %>% 
  group_by(protein_classification) %>%
  summarize(log2_relative_median_n = median(log2_relative_n))

ggplot(resampled_randomizations_relative_log2, aes(x=protein_classification,y=log2_relative_n)) + 
  geom_violin() + theme_bw() + stat_summary(fun = mean, geom="crossbar") + geom_hline(yintercept=0)

ggsave("plots/category_relative_abundance_violin_log2.pdf")

resampled_randomizations_relative = resampled_randomizations %>% left_join(actual_ssrs_by_protein_category_counts %>% rename(actual_n=n), by="protein_classification") %>% mutate(relative_n = actual_n/n)

ggplot(resampled_randomizations_relative, aes(x=protein_classification,y=relative_n, fill=protein_classification)) + 
  geom_violin() + theme_bw() + stat_summary(fun = mean, geom="crossbar") + geom_hline(yintercept=1) + scale_y_continuous(limits=c(0,1.1), breaks = 0.1*(0:11), expand=c(0,0)) + scale_fill_manual(values=okabe_ito_categories) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + theme(panel.grid = element_blank())

ggsave("plots/category_relative_abundance.pdf", width=6, height=4)

#Here we calculate the letters for significant differences between categories
#For each test calculate the difference in log2 between all categories
category_differences = data.frame()
for (c1 in categories) {
  
  set1 = resampled_randomizations_relative_log2 %>% filter(protein_classification==c1)
  for (c2 in categories) {
    
    if (c1 == c2) {
      next
    }
    set2 = resampled_randomizations_relative_log2 %>% filter(protein_classification==c2)

    category_differences = category_differences %>% 
      bind_rows(
        data.frame(
          category1 = set1$protein_classification, 
          category2 = set2$protein_classification, 
          replicate=set1$replicate, 
          log2_relative_diff = set1$log2_relative_n - set2$log2_relative_n
          )
        )
  }
}

category_comparisons = category_differences %>% group_by(category1, category2) %>% summarize(tail_upper = sum(log2_relative_diff>0), tail_lower = sum(log2_relative_diff<0), p_value = 2*(min(tail_lower, tail_upper))/(num_randomized_datasets))

category_comparisons$adj_p_value = p.adjust(category_comparisons$p_value, method="holm")
pvals = category_comparisons$adj_p_value
names(pvals) <- paste(category_comparisons$category1, category_comparisons$category2, sep = "-")

library(multcompView)
multcompLetters(pvals, threshold = 0.1)

