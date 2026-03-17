
###figure 4A
library(ggalluvial)
library(ggplot2)

ggplot(dat,aes(x = anno,y = re_freq,fill= cluster,
               stratum = cluster, alluvium = cluster))+
  geom_col(width = 0.5,color= "black")+
  geom_flow(width = 0.5,alpha = 0.8,knot.pos = 0)+
  theme_classic()+
  scale_fill_manual(values = nmf_color)+labs(x = "",y = "")+
  theme(axis.text.x = element_text(size = 12,  color = "black"))+
  theme(axis.text.y = element_text(size = 11,  color = "black"))+
  coord_flip()

###figure 4D
dat_liver_tum_filt %>%
  ggplot(aes(x = cluster, y = re_freq, fill = color))+
  geom_bar(stat = "identity")+
  coord_flip()+
  geom_text(aes(label = round(re_freq, 2)), # 
            
            color = "black", # 
            size = 3) + # 
  scale_fill_manual(values = c(cols))+
  # scale_y_continuous(breaks = breaks_values,
  #                    labels = abs(breaks_values))+
  theme_minimal()+ theme(
    panel.grid = element_blank(), # 
    panel.border = element_rect(color = "black", fill = NA, size = 0.8) # 添加边框
  )


