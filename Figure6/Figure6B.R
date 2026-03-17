

ggboxplot(dat_box_filt, x = "cluster", y = "value",
          color = "region", palette = "jco",
          add = c("jitter", "sd"), 
          line.color = "gray", 
          line.size = 0.4, 
          outlier.colour = NULL)+ 
  facet_wrap(vars(dataset), ncol = 2,
             scales = "free")+
  stat_compare_means(aes(group = region), label = "p.format",
                     method = "wilcox.test")+
    scale_color_manual(
    values=c( "#6A3D9A","#E31A1C"))+
  theme(legend.position="right")
