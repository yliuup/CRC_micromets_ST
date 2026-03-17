
dat_plot$col2 <- max(dat_plot$col)-dat_plot$col
dat_plot$row2 <- max(dat_plot$row)-dat_plot$row
p <-ggplot(dat_plot, aes(x = col2, y = row2,
                         color = anno,
                         fill = signature1_1)) +
  geom_point(shape = 21, size = 2) +  # shape 21 has both color and fill
  scale_color_manual(values = my_cols3) +  # Custom border colors
  scale_fill_gradient2(low = "#4DBBD5FF", mid = "white",high = "#E64B35FF"
                       , midpoint = 0.3
  ) +  # Gradient fill colors
  theme_classic2() +
  labs(title = unique(dat_plot$sample),
       x = "",
       y = "")

print(p)

