# Picard gene coverage plot

cov_dirs <- c('multiqc_data_chick_1000', 'multiqc_data_chick_2000',
              'multiqc_data_mouse_500', 'multiqc_data_mouse_1000',
              'multiqc_data_mouse_5000')

pic.list <- lapply(cov_dirs, function(x) {
  tibble(parse_number(t(read_tsv(
          here('picard', x, 'mqc_picard_rna_coverage_1.txt'))
          ))) %>%
          mutate(pct = c(0, seq(0, 100, 1))) %>%
          mutate(sample = x) %>%
          rename('cov' = 'parse_number(t(read_tsv(here("picard", x, "mqc_picard_rna_coverage_1.txt"))))') %>%
          filter(!is.na(cov))
})

pic <- bind_rows(pic.list)
pic$sample <- factor(pic$sample, levels = c('multiqc_data_chick_1000', 'multiqc_data_chick_2000',
                                            'multiqc_data_mouse_500', 'multiqc_data_mouse_1000',
                                            'multiqc_data_mouse_5000'))

ggplot(pic, aes(x = pct, y = cov), size = 1.5) +
  geom_line(aes(color = factor(sample))) +
  scale_color_discrete(name = 'Sample', labels = c('Chick 1K cells',
                                                   'Chick 2K cells',
                                                   'Mouse 500 cells',
                                                   'Mouse 1K cells',
                                                   'Mouse 5K cells')) +
  xlab('Percentage (%)') +
  ylab('Gene coverage') +
  theme_minimal() +
  theme(text = element_text(size = 15), legend.position = 'bottom')

# Save as 375 x 850 px to get 10 x 20 cm
