library(ggplot2)
library(tidyverse)
library(stringi)

my_theme = theme_bw() + theme(
  panel.grid.minor.x = element_blank(),
  panel.grid.minor.y = element_blank()
)

use_unknown = TRUE


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
  c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00","#CC79A7", "#999999"),
  categories
#  "grey"         = "#999999"
)


###
phages = c("T2", "T4", "T6")

categories_df = read_csv("raw_data/ssr_gene_category_metadata.csv")

phage_ssr_rates = read_csv("processed_data/all_phages_mutagenic_ssrs_in_genes.csv")

# Calculate overall rate
phage_overall_rates = phage_ssr_rates %>% summarize(ref_read_count = sum(ref_read_count), mut_read_count = sum(mut_read_count))
phage_overall_rates = phage_overall_rates %>% mutate(estimate = mut_read_count / (mut_read_count + ref_read_count))

#Also calculate overall rate holding out the outliers
phage_overall_rates_holdout = phage_ssr_rates %>% filter(adjusted_p_value > 0.05) %>% summarize(ref_read_count = sum(ref_read_count), mut_read_count = sum(mut_read_count))
phage_overall_rates_holdout = phage_overall_rates_holdout %>% mutate(estimate = mut_read_count / (mut_read_count + ref_read_count))

# Prepare for plotting significant ones
phage_ssr_rates = phage_ssr_rates %>% filter(adjusted_p_value <= 0.05) 

phage_ssr_rates = phage_ssr_rates %>% mutate(phage_pos_gene = paste0(phage, "-", ssr_position, "-", gene))
phage_ssr_rates = phage_ssr_rates %>% arrange(phage, estimate)
phage_ssr_rates$phage_pos_gene = factor(phage_ssr_rates$phage_pos_gene, levels=phage_ssr_rates$phage_pos_gene)

phage_ssr_rates = phage_ssr_rates %>% left_join(categories_df, by=c("phage","ssr_position"))

phage_ssr_rates$category[phage_ssr_rates$category=="immune"] = "unknown"

phage_ssr_rates$category = factor(phage_ssr_rates$category, levels=categories)


library(ggh4x)

phage_ssr_rates <- phage_ssr_rates %>%
  group_by(phage) %>%
  mutate(phage_pos_gene = factor(phage_pos_gene, levels = unique(phage_pos_gene))) %>%
  ungroup()

ggplot(phage_ssr_rates, aes(x = phage_pos_gene, y = estimate)) +
  geom_errorbar(aes(ymin = conf_low, ymax = conf_high),  # replace 'se' with your error column
                width = 0, color = "grey") +
  geom_point(aes(color=category), size=2.2) +  
  scale_color_manual(values = okabe_ito_categories) +
# consistent bar width
  coord_flip() +
  ggh4x::facet_grid2(
    phage ~ .,
    scales = "free_y",                # each facet only shows its positions
    space = "free"                    # panel widths proportional to number of bars
  ) +
  scale_x_discrete(drop = TRUE) +
  #scale_y_continuous(trans = "log10", limits=c(1E-4, 1E-2), expand = c(0, 0)) + 
  scale_y_continuous( limits=c(0, 3.5E-3), expand = c(0, 0)) + 
  
  geom_hline(yintercept=phage_overall_rates$estimate, color="red", linetype=2) +
  geom_hline(yintercept=phage_overall_rates_holdout$estimate, color="red") +
  
  my_theme +
#  geom_rect(aes(xmin = -Inf, xmax=Inf, ymin=conf_low, ymax=conf_high), data=phage_overall_rates, fill="red", alpha=0.5) +
  theme(
    panel.grid.major.y = element_blank(),  # remove horizontal major grid lines
    panel.grid.minor.y = element_blank(),   # remove horizontal minor grid lines
    panel.grid.major.x = element_blank(),  # remove horizontal major grid lines
    panel.grid.minor.x = element_blank(),   # remove horizontal minor grid lines
  )

ggsave("plots/all_phages_mutagenic_ssrs_in_genes.pdf", width=6, height=10)

ggplot(phage_ssr_rates, aes(x = phage_pos_gene, y = estimate)) +
  #geom_errorbar(aes(ymin = conf_low, ymax = conf_high),  # replace 'se' with your error column
               # width = 0.6, color = "grey") +
  geom_point(aes(color=category)) +  
  scale_color_manual(values = okabe_ito_categories) +
  # consistent bar width
  coord_flip() +
  ggh4x::facet_grid2(
    phage ~ .,
    scales = "free_y",                # each facet only shows its positions
    space = "free"                    # panel widths proportional to number of bars
  ) +
  scale_x_discrete(drop = TRUE) +
  #scale_y_continuous(trans = "log10", limits=c(1E-4, 1E-2), expand = c(0, 0)) + 
  scale_y_continuous( limits=c(0, 2.2E-3), expand = c(0, 0)) + 
  
  geom_hline(yintercept=phage_overall_rates$estimate, color="red", linetype=2) +
  my_theme +
  #  geom_rect(aes(xmin = -Inf, xmax=Inf, ymin=conf_low, ymax=conf_high), data=phage_overall_rates, fill="red", alpha=0.5) +
  theme(
    panel.grid.major.y = element_blank(),  # remove horizontal major grid lines
    panel.grid.minor.y = element_blank(),   # remove horizontal minor grid lines
    panel.grid.major.x = element_blank(),  # remove horizontal major grid lines
    panel.grid.minor.x = element_blank(),   # remove horizontal minor grid lines
  )

## Create plots of the total number of each type
X = read_csv("processed_data/ssr_variation.csv")

X = X %>% filter(phage %in% phages)
X$repeat_base[X$repeat_base=="A"] = "AT"
X$repeat_base[X$repeat_base=="T"] = "AT"
X$repeat_base[X$repeat_base=="G"] = "GC"
X$repeat_base[X$repeat_base=="C"] = "GC"

X =
  X %>%
  group_by(phage, ssr_position) %>%
  slice_head(n = 1) %>%
  ungroup()

all_phage= X %>% filter(ssr_length>=6)
all_phage = all_phage %>% filter(!grepl("intergenic", gene_position))

all_phage$ssr_length[all_phage$ssr_length>8] = 8

all_phage = all_phage %>% mutate(base_pair_length=paste0(repeat_base, "-", ssr_length))

all_phage = all_phage %>% group_by(phage,base_pair_length) %>% summarize(n=n()) %>% ungroup()

base_pair_length_palette = list(
  "GC-6" = "#8888FF",
  "GC-7" = "#0000EE",
  "GC-8" = "#000099",
  "AT-6" = "#FF8888",
  "AT-7" = "#EE0000",
  "AT-8" = "#990000"
)

all_phage$base_pair_length = factor(all_phage$base_pair_length, levels = rev(names(base_pair_length_palette)))

ggplot(all_phage, aes(x=phage, y=n, fill=base_pair_length)) + 
  geom_col() + theme_bw()  +
  coord_flip()  +
  scale_fill_manual(values=base_pair_length_palette) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  scale_y_continuous(limits=c(0,200), expand=c(0,0)) + theme(panel.grid = element_blank())

ggsave("plots/all_phages_num_ssrs_in_genes.pdf", height=4, width=6)
