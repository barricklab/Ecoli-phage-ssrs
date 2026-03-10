library(ggplot2)
library(tidyverse)
library(stringi)


my_theme = theme_bw() + theme(
  panel.grid.minor.x = element_blank(),
  panel.grid.minor.y = element_blank()
)

X = read_csv("processed_data/ssr_variation.csv")

X$repeat_base[X$repeat_base=="A"] = "A/T"
X$repeat_base[X$repeat_base=="T"] = "A/T"
X$repeat_base[X$repeat_base=="G"] = "G/C"
X$repeat_base[X$repeat_base=="C"] = "G/C"

###

all_phage_GC_6bp = X %>% filter(repeat_base=="G/C") %>% filter(ssr_length==6)
all_phage_GC_6bp = all_phage_GC_6bp %>% filter(!grepl("intergenic", gene_position))

all_phage_7bp = X  %>% filter(ssr_length==7)
all_phage_7bp = all_phage_7bp %>% filter(!grepl("intergenic", gene_position))


phages = c("T2", "T4", "T6")

all_overall_rates_df = data.frame()
for (this_phage in phages) {
  cat(this_phage, "\n")
  
  this_phage_AT_6bp = X %>% filter(phage==this_phage & repeat_base=="A/T") %>% filter(ssr_length==6)
  
  #no intergenic
  this_phage_AT_6bp = this_phage_AT_6bp %>% filter(!grepl("intergenic", gene_position))
  
  if (nrow(this_phage_AT_6bp) == 0) {
    next
  }
  
  # Add across replicates
  this_phage_AT_6bp_nr = this_phage_AT_6bp %>% group_by(phage, ssr_position, ssr_length, repeat_base, gene_position, gene, gene_product) %>% summarize(total_read_count = sum(total_read_count), ref_read_count = sum(ref_read_count), mut_read_count = sum(mut_read_count))
  
  this_phage_AT_6bp_nr_all = this_phage_AT_6bp_nr %>% group_by(phage) %>% summarize(total_read_count = sum(total_read_count), ref_read_count = sum(ref_read_count), mut_read_count = sum(mut_read_count))
  
  all_overall_rates_df = all_overall_rates_df %>% bind_rows(this_phage_AT_6bp_nr_all)
}

per_phage_rate_model = glm(cbind(mut_read_count, ref_read_count) ~ phage, data = all_overall_rates_df, family = binomial)
overal_rate_model = glm(cbind(mut_read_count, ref_read_count) ~ 1, data = all_overall_rates_df, family = binomial)
anova_res = anova(overal_rate_model, per_phage_rate_model, test="LRT")
anova_res

#The rates are not distinguishable for different phages: use one overall rate
overall_rate = exp(coef(overal_rate_model))

overall_rate = sum(all_overall_rates_df$mut_read_count) /(sum(all_overall_rates_df$ref_read_count) + sum(all_overall_rates_df$mut_read_count))

all_binomial = binom.test(x=sum(all_overall_rates_df$mut_read_count), n = sum(all_overall_rates_df$mut_read_count) + sum(all_overall_rates_df$ref_read_count))

all_phages_AT_6bp_nr = data.frame()

#for (this_phage in unique(X$phage)) {
for (this_phage in phages) {
  cat(this_phage, "\n")
  
  this_phage_AT_6bp = X %>% filter(phage==this_phage & repeat_base=="A/T") %>% filter(ssr_length==6)
  
  #no intergenic
  this_phage_AT_6bp = this_phage_AT_6bp %>% filter(!grepl("intergenic", gene_position))
  
  if (nrow(this_phage_AT_6bp) == 0) {
    next
  }
  
  # Add across replicates
  
  this_phage_AT_6bp_nr = this_phage_AT_6bp %>% group_by(phage, ssr_position, ssr_length, repeat_base, gene_position, gene, gene_product) %>% summarize(total_read_count = sum(total_read_count), ref_read_count = sum(ref_read_count), mut_read_count = sum(mut_read_count))

  
  # Add across everything to figure out the overall rate (assumes even coverage)
  
  this_phage_AT_6bp_nr_all = this_phage_AT_6bp_nr %>% group_by(phage) %>% summarize(total_read_count = sum(total_read_count), ref_read_count = sum(ref_read_count), mut_read_count = sum(mut_read_count))
  
  #ALT: per phage overall rate
  #overall_rate = this_phage_AT_6bp_nr_all$mut_read_count / (this_phage_AT_6bp_nr_all$ref_read_count + this_phage_AT_6bp_nr_all$mut_read_count)
  
# Loop through and assign p-values to observations

this_phage_AT_6bp_nr = this_phage_AT_6bp_nr %>%
  rowwise() %>%
  mutate(
    binomial = list(binom.test(x=mut_read_count, n = mut_read_count + ref_read_count, p = overall_rate)),
    p_value = binomial$p.value,
    estimate = binomial$estimate,
    conf_low = binomial$conf.int[1],
    conf_high = binomial$conf.int[2]
  ) %>%
  ungroup() %>%
  select(-binomial)

this_phage_AT_6bp_nr$adjusted_p_value = p.adjust(this_phage_AT_6bp_nr$p_value, method="BH")

this_phage_AT_6bp_nr = this_phage_AT_6bp_nr %>% arrange(adjusted_p_value)

all_phages_AT_6bp_nr = all_phages_AT_6bp_nr %>% bind_rows(this_phage_AT_6bp_nr)

write_csv(this_phage_AT_6bp_nr, paste0("processed_data/", this_phage, "_mutagenic_ssrs_in_genes.csv"))
}


write_csv(all_phages_AT_6bp_nr, "processed_data/all_phages_mutagenic_ssrs_in_genes.csv")
