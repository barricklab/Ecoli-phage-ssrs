library(tidyverse)

dir.create("processed_data", recursive = TRUE, showWarnings = FALSE)

file_metadata = read_csv("raw_data/sample_metadata.csv")

X = tibble(
  seq_id = character(),
  start = numeric(),
  end = numeric(),
  gene_position = character(),
  gene = character(),
  gene_product = character(),
  `1-bp` = numeric(),
  `2-bp` = numeric(),
  `3-bp` = numeric(),
  `4-bp` = numeric(),
  `5-bp` = numeric(),
  `6-bp` = numeric(),
  `7-bp` = numeric(),
  `8-bp` = numeric(),
  `9-bp` = numeric(),
  `10-bp` = numeric(),
  `11-bp` = numeric(),
  `12-bp` = numeric(),
  `13-bp` = numeric(),
  `14-bp` = numeric(),
  `15-bp` = numeric(),
  repeat_base = character(),
  )

for (i in 1:nrow(file_metadata)) {
  filename = file.path("raw_data/breseq_tabulate_CL_output", paste0(file_metadata$File_name[i],".csv"))
  cat(filename, "\n")
  Y = read_csv(filename)
  Y$phage = file_metadata$Phage[i]
  Y$replicate = file_metadata$Replicate[i]
  X = X %>% bind_rows(Y)
}
  
#None with ≥13 bp
unique(X$`13-bp`)

X$ssr_position = X$start
#X$ssr_ref_length = X$end - X$start + 1

X = X %>% select(-seq_id, -start, -end, -`13-bp`, -`14-bp`, -`15-bp`)
X = X %>% mutate(across(4:15, ~ replace_na(., 0)))


#Some SSRs have a diff length in this sample than in the reference sequence, so we need the max count
X <- X %>%
  mutate(ssr_length = max.col(across(4:15), ties.method = "first"))

# Refilter
X = X %>% filter(ssr_length>=3) #%>% select(-ssr_ref_length)

write_csv(X, "all_ssrs.csv")

X = X %>% mutate(total_read_count = rowSums(across(4:15), na.rm = TRUE))

X = X %>%
  rowwise() %>%
  mutate(ref_read_count = c_across(4:15)[ssr_length]) %>%
  mutate(mut_read_count = c_across(4:15)[ssr_length-1] + c_across(4:15)[ssr_length+1]) %>%
  ungroup()

X = X %>% select(phage, replicate, ssr_position, ssr_length, repeat_base, gene_position,	gene,	gene_product, total_read_count, ref_read_count, mut_read_count)

# mutation rate by ssr length
X %>% group_by(ssr_length) %>% summarize(mut_read_count = sum(mut_read_count), total_read_count = sum(total_read_count)) %>% mutate(mutation_frequency=mut_read_count/total_read_count)

#X = X %>% filter(ssr_length>=6)

write_csv(X, "processed_data/ssr_variation.csv")

