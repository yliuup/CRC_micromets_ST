
my_comparisons <- list(c("colon","liver"),
                       c("colon","lung"),
                       c("lung","liver"))
ggboxplot(mean_nb_score, x = "organ", y = "value",
          color = "organ", palette = "jco",
          add = c("jitter", "sd"),
          line.color = "gray",
          line.size = 0.4,
          outlier.colour = NULL) +
  facet_wrap(vars(ref), ncol = 5,
             scales = "free"
  )+
  stat_compare_means(aes(group = organ),
                     comparisons = my_comparisons, 
                     label = "p.format",
                     method = "wilcox.test") +
  # coord_flip() +
  theme(legend.position = "right") +scale_color_manual(name = '',values = c("#FF95A8FF","#FED976","#C6DBEF"))+
  # labs(title = i)+
  theme(axis.text.x=element_text(angle = 45,hjust = 1,vjust = 1),legend.position="right")+
  theme(
    axis.title.y.right = element_blank(),  #
    axis.ticks.y = element_blank(),      #
    axis.text.y = element_text(margin = margin(r = 0),colour = 'black',size = 10), #
    strip.text.y.left = element_text(angle = 0),
    axis.text.x = element_text(colour = 'black',size = 10),
    panel.spacing = unit(-0.2, "mm"),
    strip.background = element_blank(),
    panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid")
    # strip.background = element_rect(size = 1, colour = "black",fill = NA, linetype = "solid")
  )


ggboxplot(mean_dis_score, x = "organ", y = "value",
          color = "organ", palette = "jco",
          add = c("jitter", "sd"),
          line.color = "gray",
          line.size = 0.4,
          outlier.colour = NULL) +
  facet_wrap(vars(ref), ncol = 5,
             scales = "free"
  )+
  stat_compare_means(aes(group = organ),
                     comparisons = my_comparisons, 
                     label = "p.format",
                     method = "wilcox.test") +
  # coord_flip() +
  theme(legend.position = "right") +scale_color_manual(name = '',values = c("#FF95A8FF","#FED976","#C6DBEF"))+
  # labs(title = i)+
  theme(axis.text.x=element_text(angle = 45,hjust = 1,vjust = 1),legend.position="right")+
  theme(
    axis.title.y.right = element_blank(),  #
    axis.ticks.y = element_blank(),      #
    axis.text.y = element_text(margin = margin(r = 0),colour = 'black',size = 10), #
    strip.text.y.left = element_text(angle = 0),
    axis.text.x = element_text(colour = 'black',size = 10),
    panel.spacing = unit(-0.2, "mm"),
    strip.background = element_blank(),
    panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid")
    # strip.background = element_rect(size = 1, colour = "black",fill = NA, linetype = "solid")
  )

