library(ggplot2)
library(tidyverse)
library(stringi)

dir.create("plots", recursive = TRUE, showWarnings = FALSE)

my_theme = theme_bw() + theme(
  panel.grid.minor.x = element_blank(),
  panel.grid.minor.y = element_blank()
)

X = read_csv("processed_data/ssr_variation.csv")

X$repeat_base[X$repeat_base=="A"] = "A/T"
X$repeat_base[X$repeat_base=="T"] = "A/T"
X$repeat_base[X$repeat_base=="G"] = "G/C"
X$repeat_base[X$repeat_base=="C"] = "G/C"

by_phage_replicate_repeat_base = X %>% group_by(phage, replicate, repeat_base, ssr_length) %>% summarize(total_read_count = sum(total_read_count), ref_read_count = sum(ref_read_count))
by_phage_replicate_repeat_base = by_phage_replicate_repeat_base %>% mutate(mut_read_count = total_read_count-ref_read_count, mut_frequency = mut_read_count/total_read_count)

by_phage_repeat_base = X %>% group_by(phage, repeat_base, ssr_length) %>% summarize(total_read_count = sum(total_read_count), ref_read_count = sum(ref_read_count))
by_phage_repeat_base = by_phage_repeat_base %>% mutate(mut_read_count = total_read_count-ref_read_count, mut_frequency = mut_read_count/total_read_count)

by_phage_repeat_base = by_phage_repeat_base %>%
  rowwise() %>%
  mutate(
    ci = list(binom.test(mut_read_count, total_read_count)$conf.int),
    lower_95_mut_frequency = ci[[1]],
    upper_95_mut_frequency = ci[[2]]
  ) %>%
  ungroup() %>%
  select(-ci)

write_csv(by_phage_repeat_base, "processed_data/ssrs_by_phage_repeat_base.csv")

#plot each replicate separately with error bars
ggplot(by_phage_replicate_repeat_base, aes(y=mut_frequency, x=ssr_length, color=repeat_base)) +
  scale_y_log10() +
  scale_x_continuous(
    breaks = 3:12,
    labels = 3:12
  ) +
  geom_point(aes(group=replicate)) +
  geom_line(data=by_phage_repeat_base) +
  geom_errorbar(data=by_phage_repeat_base, aes(ymin=lower_95_mut_frequency, ymax=upper_95_mut_frequency)) +
  facet_wrap(~ phage, nrow = 2) +
  my_theme

ggsave("plots/ssrs_by_repeat_base_and_ssr_length.pdf", width=5, height=3.5)

# just T2
ggplot(by_phage_replicate_repeat_base %>% filter(phage=="T2"), aes(y=mut_frequency, x=ssr_length, color=repeat_base)) +
  scale_y_log10() +
  scale_x_continuous(
    breaks = 3:12,
    labels = 3:12
  ) +
  geom_point(aes(group=replicate)) +
  geom_line(data=by_phage_repeat_base%>% filter(phage=="T2")) +
  geom_errorbar(data=by_phage_repeat_base%>% filter(phage=="T2"), aes(ymin=lower_95_mut_frequency, ymax=upper_95_mut_frequency), width=0.2) +
  my_theme +
  labs(
    x = "SSR length (bp)",
    y = "Mutation frequency",
    color = "SSR base"
  )

ggsave("plots/T2_by_repeat_base_and_ssr_length.pdf", width=5, height=3.5)

#plot only averages
ggplot(by_phage_repeat_base, aes(y=mut_frequency, x=ssr_length, color=repeat_base)) +
  scale_y_log10() +
  scale_x_continuous(
    breaks = 3:12,
    labels = 3:12
  ) +
  geom_point() +
  geom_line(data=by_phage_repeat_base) +
  geom_errorbar(data=by_phage_repeat_base, aes(ymin=lower_95_mut_frequency, ymax=upper_95_mut_frequency)) +
  facet_wrap(~ phage, nrow = 2) +
  my_theme



#plot on top of one another
ggplot(by_phage_repeat_base, aes(y=mut_frequency, x=ssr_length, color=phage)) +
  scale_y_log10() +
  scale_x_continuous(
    breaks = 3:12,
    labels = 3:12
  ) +
  geom_point() +
  geom_line() +
  #geom_errorbar(data=by_phage_repeat_base, aes(ymin=lower_95_mut_frequency, ymax=upper_95_mut_frequency), width=0.2) +
  facet_wrap(~repeat_base) +
  my_theme

## Number of tracts of each type

ggsave("plots/all_phages_by_repeat_base_and_ssr_length.pdf", width=5, height=3.5)

