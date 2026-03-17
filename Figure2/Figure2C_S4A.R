
library(tidyverse)
library(patchwork)
library(stringr)


# merge all cnvkit .cns files

# all_cns_files <- list.files('cnv_calling_result/', pattern = 'recal\\.cns$')
# sam_names <- all_cns_files %>% str_remove('_recal.cns')
# 
# cnvkit_cns <- tibble(Sample = sam_names) %>% 
#   mutate(Data = map(Sample, ~read_tsv(str_c('cnv_calling_result/', ., '_recal.cns')))) %>% 
#   unnest(cols = Data)
# 
# save(cnvkit_cns, file = 'cnv_calling_result_merged/cnvkit_cns.Rdata')


# merge all cnvkit .call.cns files

# all_call_cns_files <- list.files('cnv_calling_result/', pattern = 'call\\.cns$')
# sam_names <- all_call_cns_files %>% str_remove('_recal.call.cns')
# 
# cnvkit_call_cns <- tibble(Sample = sam_names) %>% 
#   mutate(Data = map(Sample, ~read_tsv(str_c('cnv_calling_result/', ., '_recal.call.cns')))) %>% 
#   unnest(cols = Data)
# 
# save(cnvkit_call_cns, file = 'cnv_calling_result_merged/cnvkit_call_cns.Rdata')


load('cnvkit_call_cns.Rdata')


# Only genes used in inferCNV ---------------------------------------------


Sys.setenv("VROOM_CONNECTION_SIZE" = 5000000)

infercnv_mat <- read_tsv('infercnv.observations.txt', col_names = FALSE)
genes <- infercnv_mat$X1[-1] %>% str_extract("^[^\\s]+")

gene_pos_dat <- read_tsv('gencode_v21_gene_pos.txt', col_names = c('Hugo_Symbol', 'chr', 'start', 'end'))

infercnv_gene_pos <- gene_pos_dat %>% filter(Hugo_Symbol %in% genes)

save(infercnv_gene_pos, file = 'infercnv_gene_pos.Rdata')

load("cnvkit_call_cns.Rdata")
load('gene_hm_plot_dat.Rdata')
load("infercnv_gene_pos.Rdata")

cnvkit_gene_cnv_dat <- cnvkit_call_cns %>% 
  mutate(CNV = 2^log2) %>% 
  select(Sample, gene, CNV) %>% 
  filter(str_detect(gene, paste(infercnv_gene_pos$Hugo_Symbol, collapse = "|"))) %>% 
  separate_rows(gene, sep = ",") %>%
  mutate(gene = str_trim(gene)) %>%
  filter(gene %in% infercnv_gene_pos$Hugo_Symbol) %>% 
  distinct(Sample, gene, .keep_all = TRUE)

gene_hm_plot_dat <- infercnv_gene_pos %>% 
  left_join(cnvkit_gene_cnv_dat, by = c('Hugo_Symbol' = 'gene')) %>% 
  filter(!is.na(CNV)) %>% 
  mutate(sam_label = str_extract(Sample, "^[^-]+"))

save(gene_hm_plot_dat, file = 'CNV_gene_heatmap_plotting/gene_hm_plot_dat.Rdata')

# plot_mat <- gene_hm_plot_dat %>% 
#   select(Hugo_Symbol, Sample, CNV) %>% 
#   pivot_wider(names_from = Hugo_Symbol, values_from = CNV)
# 
# sam_info <- cnvkit_call_cns %>% distinct(Sample) %>% mutate(sam_label = str_extract(Sample, "^[^-]+"))
# gene_chr_label <- infercnv_gene_pos %>% semi_join(gene_hm_plot_dat, by = 'Hugo_Symbol') %>% select(Hugo_Symbol, chr)
# 
# all_sam_label <- sam_info %>% distinct(sam_label) %>% arrange(sam_label) %>% pull(sam_label)

sam_order <- c('5-C1micro', '5-C1', '5-C2', '5-C3', 
               '09-C1micro', '09-C1', '09-C3', 
               '32-C1', '32-C2', '32-C3', '32-C6', 
               '36-C3micro', '36-C3', '36-C4', '36-C5', 
               '37-C7', '37-C8', 
               'Ix8889-C1micro', 'Ix8889-C1')
chr_order <- as.character(1:22)
gene_order <- infercnv_gene_pos$Hugo_Symbol


plot_gene_heatmap <- function(sam_lab){
  
  # sam_lab <- '5'
  
  sin_plt <- gene_hm_plot_dat %>% 
    filter(sam_label %in% sam_lab) %>% 
    mutate(chr = str_remove(chr, 'chr'), 
           chr = factor(chr, levels = chr_order), 
           Hugo_Symbol = factor(Hugo_Symbol, levels = gene_order),
           Sample = factor(Sample, levels = rev(sam_order)))
  sin_plt$chr <- factor(sin_plt$chr, levels = paste0(chr_order))
  sample_order <- sin_plt %>% 
    distinct(Sample) %>% 
    arrange(Sample) %>% 
    mutate(y = as.numeric(factor(Sample, levels = rev(unique(Sample)))))
  
  row_lines <- sample_order$y + 0.5
  
  sin_plot <- sin_plt %>% 
    ggplot(aes(Hugo_Symbol, Sample, fill = CNV)) +
    geom_tile() +
    geom_hline(yintercept = row_lines, color = "black", linewidth = 0.2) +
    scale_fill_gradientn(
      colors = c('#16375f', '#005a8d', '#359dc7', '#c7dcea', '#f8f8f9',
                 '#edbdb8', '#d96f59', '#982b1f', '#5f1a0c'),
      limits = c(0.7, 1.3),midpoint = 1,
      # colors = rev(cnv.cols),
      oob = scales::squish,
      name = "CNV"
    ) +
    labs(x = NULL, y = NULL) +
    facet_grid(.~chr, space = 'free', scales = 'free') +
    theme(axis.ticks = element_blank(), 
          axis.text.x = element_blank(), 
          panel.grid = element_blank(), 
          panel.background = element_blank(), 
          panel.spacing = unit(0, "pt"), 
          strip.background = element_blank(), 
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3)) +
    scale_y_discrete(expand = c(0, 0))
  
  
   pdf('sam_5_logR_heatmap_genes.pdf', width = 14, height = 1.5)
  
  sin_plot
  
   dev.off()
  
}

hm_plt_5 <- plot_gene_heatmap('5')
hm_plt_09 <- plot_gene_heatmap('09')
hm_plt_32 <- plot_gene_heatmap('32')
hm_plt_36 <- plot_gene_heatmap('36')
hm_plt_37 <- plot_gene_heatmap('37')
hm_plt_Ix8889 <- plot_gene_heatmap('Ix8889')



pdf('all_sam_CNV_heatmap_genes.pdf', width = 14, height = 7)

hm_plt_5 + hm_plt_09 + hm_plt_32 + hm_plt_36 + hm_plt_37 + hm_plt_Ix8889 + plot_layout(ncol = 1, guides = "collect", heights = c(4/19, 3/19, 4/19, 4/19, 2/19, 2/19))

dev.off()

