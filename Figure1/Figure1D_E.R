
load("crc_merge_meta.RData")
aa <- crc_merge_meta[,c("layer1","layer2","layer3")]

library(dplyr)
library(plyr)
library(ggsankey)
library(dplyr)
library(ggplot2)
df <- aa %>%  make_long(layer1,layer2,layer3)

df <- as.data.frame(df)

ggplot(df, aes(x = x, next_x = next_x, 
               node = node,       
               next_node = next_node,  
               fill = node,     
               label = node)) + 
  scale_fill_manual(values = my_cols3)+
  geom_sankey(flow.alpha = 0.5, node.color = 1) +  
  geom_sankey_label(size = 4, color = 1) + 
  theme_sankey(base_size = 16)



levels2 <- data.frame(key = names(table(aa$layer3)),
                      values = table(aa$layer3),
                      subject = 2)

library(ggplot2)
library(ggpubr)
ggplot(data=levels2, aes(x=values.Var1, y = values.Freq,fill = values.Var1)) +
  geom_bar(stat = "identity",width = 0.8)+
  scale_fill_manual(values=my_cols3)+
  geom_text(aes(label=values.Freq_raw), vjust=0.3, size=3.5)+
  # facet_grid( .~subject, scales="free")+
  # facet_grid(rows = vars(subject))+
  coord_flip()+labs(y='Counts (log2)')+
  theme_classic2()+
  geom_hline(yintercept = median(as.numeric(levels2$values.Freq)), linetype = "dashed", color = "red", linewidth = 0.5)


crc_merge_meta$nFeature_Spatial2 <- log2(crc_merge_meta$nFeature_Spatial+1)
crc_merge_meta$nCount_Spatial2 <- log2(crc_merge_meta$nCount_Spatial+1)

ggplot(data = crc_merge_meta, aes(x = layer3, 
                                      y = nFeature_Spatial2, 
                                      fill = layer3, 
                                      color = layer3)) + 
  geom_boxplot(width = 0.6, 
               outlier.shape = NA, 
               alpha = 0.5, 
               size = 1) + 
  scale_color_manual(values = my_cols3) + 
  scale_fill_manual(values = my_cols3) + 
  coord_flip() + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size = 10, colour = "black"), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.position = 'top') + 
  labs(y = 'nFeature_Spatial (log2)') + 
  theme_classic()



ggplot(data = crc_merge_meta, aes(x = layer3, 
                                      y = nCount_Spatial2, 
                                      fill = layer3, 
                                      color = layer3)) + 
  geom_boxplot(width = 0.6, 
               outlier.shape = NA, 
               alpha = 0.5, 
               size = 1) +  
  scale_color_manual(values = my_cols3) + 
  scale_fill_manual(values = my_cols3) + 
  coord_flip() + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size = 10, colour = "black"), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.position = 'top') + 
  labs(y = 'nCount_Spatial (log2)') + 
  theme_classic()