
my_comparisons <- list(
  c("Liver micrometastasis tumor",  "Liver micrometastasis stroma"),
  c("Liver micrometastasis tumor","Liver macrometastasis tumor"),
  
  
  c("Liver micrometastasis tumor",  "Liver macrometastasis stroma"),
  c("Liver micrometastasis tumor",  "Hepatic lobule")
  
  
)



ggplot(crc_merge_filt_meta_liver, aes(x = anno,
                                       y = six_gene, 
                                       fill = anno)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.7) +  
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA) +  
  scale_fill_manual(values =my_cols3 ) +  
  theme_classic() + 
  facet_wrap(vars(sample), ncol = 25)+
  labs(x = "", y = "six_gene") + 
  stat_compare_means(aes(group = anno),
                     comparisons = my_comparisons,
                     # label.y = 1.2,
                     label = "p.signif",
                     method = "wilcox.test" ,
                     # method.args = list(var.equal = T),
                     paired =F)+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12, color = "black")
  )
