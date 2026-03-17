
library('GSVA')
library("limma")
library(GSEABase)


gmt_file2 <- getGmt("./merged_fixed.gmt")
list_gmt <- gmt_file2@.Data

new_list2 <- list()

for (i in c(1:160)){
  a <- list(list_gmt[i][[1]]@geneIds  )
  names(a) <- list_gmt[i][[1]]@setName
  new_list2 <- append(new_list2,a)
  
}


for (i in c(1:length(new_list2))){
  
  dat_both <- AddModuleScore(dat_both,
                             features = new_list2[i] ,
                             name = names(new_list2)[i])
}

dat_both_meta <- dat_both@meta.data
save(dat_both_meta, file = "dat_both_meta.RData")


dat_both_meta_melt <- melt(dat_both_meta)
ggboxplot(dat_both_meta_melt,
          x="anno",y="value",
              color = "anno",   line.color = "gray", line.size = 0.4,outlier.colour=NULL) +
 
  facet_wrap(~variable, scales = "free_y",ncol =3)+
  theme(axis.text.x=element_text(angle = 45,hjust = 1,vjust = 1),legend.position="right")+
  stat_compare_means(
    aes(group = region), 
    comparisons = my_comparisons, 
    # label = "p.signif", 
     method = "t.test", #wilcox.test
    paired = FALSE
  )+theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), 
          legend.position = "right")+
  scale_color_manual(name = '',values = my_cols3)

