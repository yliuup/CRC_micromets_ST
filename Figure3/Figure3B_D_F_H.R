
#Figure 3B
library(viridis)
cols <- turbo(30)[2:28]
# obj <- dat_list[[j]]
size = 1.5
p <- SpatialFeaturePlot(dat, 
                        features = j,
                        ncol = 1,stroke = NA,image.alpha = 1,
                        pt.size.factor = size)+scale_fill_gradientn(colours = cols)
theme(legend.position = "top")

#Figure 3D

library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)
library(ggsci)
library(RColorBrewer)
library(scales)

clusterPlot(sce,color = NA, label = group,
                  palette = my_cols
                  
)+ labs(title=paste0(i))+theme(legend.position = "none")



#Figure 3F
png(paste0("liver_hd", "_anno_dimplot.png"), width = 8, height = 7, res = 400, units = "in")
DimPlot(
  liver_hd_filt3,
  group.by = "annotation",
  pt.size = 0.1,
  label = TRUE,
  label.size = 5,
  raster = FALSE,
  ncol = 1,
  cols = my_cols
)

dev.off()


#Figure 3H
ggplot(meta, aes(x = x, y = y,color = annotation)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = liver_color) +
  coord_fixed() +
  labs(title = "", x = "X", y = "Y") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),         
    panel.background = element_blank(),   
    axis.ticks = element_blank(),         
    axis.text = element_blank(),          
    panel.border = element_rect(color = "black", fill = NA, size = 1)  
  )+guides(color = guide_legend(override.aes = list(size = 4)))
