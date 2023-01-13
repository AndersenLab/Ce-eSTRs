setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/"))
#setwd("../")

### load R packages####

library(tidyverse)

### ggplot theme ####

theme_cust <- theme_bw() + 
  theme(plot.title = ggplot2::element_text(size=12,  color = "black"),
        legend.text = ggplot2::element_text(size=12,  color = "black"),
        legend.title =  ggplot2::element_text(size=12,  color = "black"),
        axis.title =  ggplot2::element_text(size=12,  color = "black"),
        axis.text =  ggplot2::element_text(size=12,  color = "black"),
        strip.text = ggplot2::element_text(size=12, vjust = 1,  color = "black"),
        strip.background = ggplot2::element_blank(), 
        panel.grid = ggplot2::element_blank(),
        text = ggplot2::element_text(family="Helvetica"))


ancestry.colours <- c('red',"gold2","plum4",     "lightskyblue2", 
                      "springgreen4", "lightpink2",  "deepskyblue4", 
                      "yellow3",  "yellow4",  
                      'black', 'cornflowerblue', 'magenta', 'darkolivegreen4', 
                        'tan4', 'darkblue', 'yellowgreen', "tan1",
                      'darkgray', 'wheat4', '#DDAD4B', 'chartreuse','seagreen1',
                      'moccasin', 'mediumvioletred', 'cadetblue1',"darkolivegreen1" ,"#7CE3D8",
                      "gainsboro","#E69F00","#009E73", "#F0E442", "sienna4", "#0072B2", 
                      "mediumpurple4","#D55E00", "burlywood3","gray51","#CC79A7","gray19", "firebrick") 



############# Figure  1    ###############
#            local eSTRs                 #
##########################################

###### fig_1a ######
load("../processed_data/local_eSTRs_LRT.RData") 

factoral_eSTR_lrt_log <- factoral_eSTR_lrt_supp2 %>% 
  dplyr::group_by(str_genotype) %>% 
  dplyr::add_count(name = "N_tests") %>% 
  dplyr::mutate(real_lrt_p=-log10(real_lrt_p),
                perm_lrt_p=-log10(perm_lrt_p),
                distance=distance/1e6,
                bf_thres =  -log10(0.05/N_tests))


## plot 
factoral_eSTR_lrt_df <- factoral_eSTR_lrt_log %>% 
  dplyr::select(transcript,pSTR,distance,real_lrt_p,perm_lrt_p,str_genotype,bf_thres,nSTR_pheno) %>% 
  tidyr::gather(type,padj,-transcript,-pSTR,-distance,-str_genotype,-bf_thres,-nSTR_pheno) %>% 
  dplyr::mutate(type=ifelse(type=="real_lrt_p","Real STR","Permuted STR")) %>% 
  dplyr::ungroup()  

perc_eSTR <- factoral_eSTR_lrt_log %>% 
  dplyr::select(transcript,pSTR,bin_center,real_lrt_p,str_genotype,bf_thres ) %>% 
  dplyr::mutate(sig=ifelse(real_lrt_p > bf_thres,"pass","no")) %>% 
  dplyr::group_by(bin_center,sig) %>% 
  dplyr::add_count(name = "sig_count") %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(bin_center ) %>% 
  dplyr::add_count(name = "bin_test") %>% 
  dplyr::mutate(sig_perc= 1000*sig_count/bin_test,
                distance=bin_center/1e6 ) %>% 
  dplyr::filter(abs(bin_center)>10 ) %>%  
  dplyr::ungroup() %>% 
  dplyr::select(distance,sig,sig_perc,str_genotype,bf_thres) %>% 
  dplyr::distinct() %>% 
  tidyr::spread(sig,sig_perc,fill = 0) %>% 
  dplyr::rename( padj=pass)  

fig_1a <- ggplot() +
  geom_point(data=factoral_eSTR_lrt_df, alpha=0.8,size=0.2,shape=20, aes(x=distance,y=padj,color=type)) +
  geom_line(data=perc_eSTR,  aes(x=distance,y=padj),color="darkorange1",size=1)+
  scale_color_manual(values=c("darkgray","lightskyblue2")) + 
  geom_hline(data=perc_eSTR,aes(yintercept = bf_thres),
             color = "red",
             linetype=3 )+
  theme_cust+
  theme(legend.position = "none") +
  labs(x="Distance to TSS (Mb)", 
       color="")+
  theme(legend.position = "none", panel.spacing = unit(2,"line")) + 
  scale_x_continuous(expand = c(0 , 0.01) ) + 
  scale_y_continuous( expand = c(0, 0) ,name = expression(-log[10](italic(p))), 
                      sec.axis = sec_axis( trans=~./10, name="Mean percentile of\nsignificant STRs / 20 kb") )+
  facet_grid(.~str_genotype,scales="free")
  


perc_sig <- factoral_eSTR_lrt_log %>% 
  dplyr::filter(str_genotype=="STR genotype" & nSTR_pheno==6)
  dplyr::ungroup() %>% 
  dplyr::mutate(sig=ifelse(real_lrt_p > bf_thres,"pass","no")) %>% 
  dplyr::group_by(nSTR_pheno,sig) %>% 
  dplyr::count() %>% 
  tidyr::spread(sig,n) %>% 
  dplyr::mutate(sig_perc=100*pass/(pass+no))







 
fig_S1_new <- ggplot() +
  geom_point(data=factoral_eSTR_lrt_log, alpha=0.8,size=0.2,shape=20, aes(x=distance,y=real_lrt_p ,color=factor(nSTR_pheno)))  +
  geom_hline(data=perc_eSTR,aes(yintercept = bf_thres),
             color = "black",
             linetype=3 )+
    theme_cust+
  theme(legend.position = "none") +
  labs(x="Distance to TSS (Mb)", 
       y=expression(-log[10](italic(p))) )+
  theme(legend.position = "none", panel.spacing.x  = unit(2,"line")) + 
  scale_x_continuous(expand = c(0 , 0.01) ) + 
    facet_grid(nSTR_pheno~str_genotype,scales="free")+
  scale_color_manual(values=c("2"="gray","3"="#D55E00", "4"="mediumpurple4", "5"="#0072B2","6"="red","7"="lightskyblue2"))


ggsave(fig_S1_new, filename = paste("../figures/supp_fig1_new.png",sep = ""), units = "mm",height = 160, width = 170)



###### fig_1b ######

Lrt_localeSTR_eQTL_varexp <- data.table::fread("../processed_data/Lrt_localeSTR_eQTL_varexp.tsv")

Lrt_localeSTR_eQTL_varexp$LD_level2<- factor(Lrt_localeSTR_eQTL_varexp$LD_level,levels = c("High LD","Moderate LD","Low LD"))
 
fig_1b <- ggplot( ) +
  geom_point(data=subset(Lrt_localeSTR_eQTL_varexp, nSTR_pheno %in% c(2,3)),
             aes(x=eQTL_var_exp,y=STR_var_exp,color=factor(nSTR_pheno)),size=1,alpha=0.5)+
  geom_point(data=subset(Lrt_localeSTR_eQTL_varexp, nSTR_pheno %in% c(4,5,6)),
             aes(x=eQTL_var_exp,y=STR_var_exp,color=factor(nSTR_pheno)),size=1)+
  theme_cust+
  theme(panel.spacing.x = unit(2,"line"),
        panel.spacing.y = unit(1,"line"),
        legend.position = "bottom")+
  geom_abline(intercept=0,slope=1,colour="black",linetype=2)+
  scale_y_continuous(breaks=c(0,0.5,  1),#expand = c(0, 0), 
                     limits = c(0,1))  +
  scale_x_continuous(breaks=c(0, 0.5, 1),#expand = c(0, 0), 
                     limits = c(0,1)) +
  facet_grid( str_genotype ~LD_level2 ) +
  ylab("Variance explained by eSTRs")+
  xlab("Variance explained by local eQTL") +
  labs(color="Number of STR alleles") +
  scale_color_manual(values=c("2"="gray","3"="#D55E00", "4"="mediumpurple4", "5"="#0072B2","6"="black"))

###### fig_1  ######

fig_1 <- cowplot::plot_grid(fig_1a, fig_1b,  
                           labels = c('A', 'B' ), 
                           label_size = 12, 
                           label_fontfamily="Helvetica",
                           rel_heights = c(1,1.5),
                           axis = "lr",
                           #    align = "v",
                           nrow = 2)

ggsave(fig_1, filename = paste("../figures/fig_1.png",sep = ""), units = "mm",height = 160, width = 170)



############# Figure  2    ###############
#            STR_24584                  #
##########################################

###### fig_2a ######
str_exp_pxg <- data.table::fread("../processed_data/str_exp_pxg.tsv")

data_fig_2a <- str_exp_pxg %>% 
  dplyr::filter(pSTR=="STR_24584")

fig_2a <- ggplot(data_fig_2a,aes(x=as.factor(genotype),y=expression ))+
  geom_jitter(shape=21,position=position_jitter(0.2),aes(fill=STRallele) ) +
  geom_boxplot( outlier.shape = NA,alpha=0.5) +
  theme_cust+
  theme(strip.text = ggplot2::element_text(size=12, vjust = 1,  color = "black",face = "italic"),
        legend.position = "none")+
  ggplot2::scale_fill_manual(values = c("REF"="orange","ALT"="blue") )+
  xlab("Length of STR_24584") + ylab("Expression") +  
  facet_grid(.~transcript ,scales = "free") +
  scale_y_continuous(breaks=c(-1,1,3,5),expand = c(0, 0), limits = c(-1,6))


ggsave(fig_2a, filename = paste( "../figures/fig_2a.png",sep = ""), units = "mm",height = 60, width = 170)






############# Figure  3    ###############
#          distant eSTRs                 #
##########################################

###### fig_3a ######

data_fig_3a <- data.table::fread("../processed_data/Lrt_distanteSTR_eQTL_varexp.tsv") %>% 
  dplyr::mutate(str_genotype="STR length")

data_fig_3a$LD_level2<- factor(data_fig_3a$LD_level,levels = c("High LD","Moderate LD","Low LD"))

fig_3a <- ggplot()+
  geom_point(data = subset(data_fig_3a, nSTR_pheno %in% c(2,3)),
             aes(x=eQTL_var_exp,y=STR_var_exp,color=factor(nSTR_pheno)), size=1,alpha=0.5)+
  geom_point(data = subset(data_fig_3a, nSTR_pheno %in% c(4,5 )),
             aes(x=eQTL_var_exp,y=STR_var_exp,color=factor(nSTR_pheno)), size=1)+
  theme_cust+
  theme(panel.spacing.x = unit(2,"line"),
        panel.spacing.y = unit(1,"line"),
        legend.position = "bottom")+
  geom_abline(intercept=0,slope=1,colour="black",linetype=2)+
  scale_y_continuous(breaks=c(0,0.5,  1), limits = c(0,1))  +
  scale_x_continuous(breaks=c(0, 0.5, 1), limits = c(0,1)) +
  facet_grid( str_genotype ~LD_level2 ) +
  ylab("Variance explained\nby eSTRs")+
  xlab("Variance explained by distant eQTL") +
  labs(color="Number of STR alleles") +
  scale_color_manual(values=c("2"="gray","3"="#D55E00", "4"="mediumpurple4", "5"="#0072B2" ))

###### fig_3b ######
data_fig_3b <- data.table::fread("../processed_data/Lrt_hotspot_eSTRs.tsv") %>% 
  dplyr::filter(n_med_STR>4)


hotspot_pos <- data_fig_3b %>% 
  dplyr::distinct(hotspot_Chr, hotspot_cM_center, Hotspot, merged_Hotspot_QTL_count)  

newrows_size <- data.frame(p1=0,
                           hotspot_Chr=c("I","II","III","IV","V","X"  ),
                           hotspot_cM_center=c(51.5,42.5,49,44,45.5,45.5) )

fig_3b <- ggplot() + 
  geom_bar(data=hotspot_pos,aes(x=hotspot_cM_center,y=merged_Hotspot_QTL_count/3), stat='identity',color="black")+
  geom_point(data=na.omit(data_fig_3b),aes(x=hotspot_cM_center,y=p1),size=1,shape=25,color="lightskyblue2")+
  geom_segment(data=newrows_size,aes(x = 0, y = p1, xend = hotspot_cM_center,yend = p1 ), size = 2.5, alpha = 0) +
  facet_grid(.~hotspot_Chr,scales = "free_x" ) +
  xlab("Hotspot position (cM)")  + 
  theme_cust  +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        panel.spacing = unit(0.1,"line") )+
  scale_y_continuous( name = "Percentage of\ndistant eQTL (%)", sec.axis = sec_axis( trans=~.*3, name="Total number of\ndistant eQTL\nin each hotspot"), expand = c(0, 0),   limits = c(0,70) ) +
  scale_x_continuous(expand = c(0, 0) ) 




fig_3 <- cowplot::plot_grid(fig_3a, fig_3b,   
                            labels = c('A', 'B' ), 
                            label_size = 12, 
                            label_fontfamily="Helvetica",
                            rel_heights = c(1.4,1),
                            axis = "l",
                            nrow = 2)

ggsave(fig_3, filename = paste("../figures/fig_3.png",sep = ""), units = "mm",height = 120, width = 170)



############# Figure  4    ###############
#       STR mutation variation           #
##########################################


###### fig_4a ######
 


med  <- data.table::fread("../processed_data/STR_nema_med207.tsv")

med_sig <- med %>% dplyr::filter(!is.na(q99_mediator) ) %>% 
  dplyr::arrange(desc(multi_abs_est))

med_sig$q99_mediator2 <- factor(  med_sig$q99_mediator, levels = unique(med_sig$q99_mediator) )

med_other <- med %>% dplyr::filter(is.na(q99_mediator))  

fig_4a <- ggplot() +
  geom_point(data=med_other, aes( x=e_peak/1e6,
                                  y=multi_abs_est ), color =  "gray80" ,size=0.3) +
  geom_point(data=med_sig, aes( x=e_peak/1e6,
                                y=multi_abs_est, color = q99_mediator2 ),size=1 ) +
  geom_hline( data=med_other, aes(yintercept = q99) , color = "grey")+
  geom_point(data=med_sig, aes( x=e_peak/1e6,
                                y=multi_abs_est, color = q99_mediator2 ),size=1 ) +
  labs(x = "Genomic position (Mb)", 
       y = "Mediation estimate",
       color = "Mediator gene" )+
  facet_grid(.~gwchr,scales = "free_x") +
  theme_cust +
  theme(#legend.text =  ggplot2::element_text(size=12,  color = "black",face = "italic"), 
        legend.title = element_blank(),
        legend.margin = margin(0,0,0,-20),
        legend.text = element_text(size=12,  color = "black",face = "italic" ,
                                   margin = margin(r = 0.1,l = 0.1, t = 0.1,b = 0.1,unit = "pt")),
       # legend.spacing.x = unit(0.01, "cm"), 
      #  legend.spacing.x = unit(0.001, 'cm'),
        legend.key = element_rect(size = 3, fill = "white", colour = NA), 
       # legend.key.size = unit(0.3, "cm"),
        legend.position = "bottom",
        panel.spacing = unit(0.7,"line"))+
  scale_color_manual(values=ancestry.colours,guide = guide_legend(nrow = 6) )  + 
  scale_y_continuous(expand = c(0, 0.01) )   + 
  scale_x_continuous(breaks   = c(0, 5,10,15) )  



###### fig_4b ######
#
top_mediators_cor <- data.table::fread("../processed_data/STR_mutation_trait_mediators_cor.tsv") %>% 
  dplyr::mutate(transcript=paste(ext_gene,"\n",transcript))%>% 
 # dplyr::left_join(strain_ALT_frac )%>% 
  dplyr::mutate(pp="p ",ppp="r ") %>% 
  dplyr::filter(ext_gene=="ctl-1")

top_mediators_cor$transcript2 <- factor(  top_mediators_cor$transcript, levels = c("ctl-1 \n Y54G11A.6.1","ctl-1 \n Y54G11A.6.2","F59E12.15 \n F59E12.15.1", "F59E12.15 \n F59E12.15.2") )


fig_4b <- ggplot(top_mediators_cor,aes(y=Total_mutation,x=exp,color=ext_gene ))+
  geom_point(shape=19,alpha=0.8 ,size=0.3)+
  #  scale_fill_gradient(high = "#D7263D", low = "#0072B2",name="Expression" )+
#  scale_color_manual(values = c("#E7B800", "#FC4E07")) +
  #  scale_color_manual(values = c("F59E12.15"="gold2","ctl-1"="plum4")) +
  facet_grid(.~transcript2 ) +
  theme_cust+
  theme(legend.position = "none") +
  xlab("Expression")+
  ylab("STR variation")+
  ggplot2::theme( strip.text.x = ggplot2::element_text(size=12, vjust = 1,  color = "black",face = "italic"),
                  strip.text.y = element_blank(), panel.spacing = unit(1,"line")) +
  geom_text(data= subset(top_mediators_cor, strain=="AB1"), aes(label = paste0(" : ",pearson_cor), x= 3.8, y=33 ),color="gray6" )+ 
  geom_text(data= subset(top_mediators_cor, strain=="AB1" ), aes(label = paste0(" : ",pvalue ), x=4, y=25  ),color="gray6" )+ 
  geom_text(data= subset(top_mediators_cor, strain=="AB1" ), aes(label = pp, x=3.1, y=25 ,fontface=3 ),color="gray6"  )+
  geom_text(data= subset(top_mediators_cor, strain=="AB1" ), aes(label = ppp, x=3.1, y=33 ,fontface=3 ),color="gray6"  )+
  scale_x_continuous(breaks=c(-1,1,3,5 )  )+
  scale_color_manual(values=ancestry.colours ) 

 
###### fig_4 ######





fig_4 <- cowplot::plot_grid(  fig_4a,fig_4b,  
                            labels = c('A', 'B'  ), 
                            label_size = 12, 
                            label_fontfamily="Helvetica",
                            rel_heights =  c( 1.6,1 ),
                            axis = "l",
                            # align = "v",
                            nrow = 2)

ggsave(fig_4, filename = paste( "../figures/fig_4.png",sep = ""), units = "mm",height = 150, width = 170)





############# Figure  4    #############
#           MA lines mutation rate     #
#########################################


data_fig_5 <- data.table::fread("../processed_data/MA_u.tsv") 


data_fig_5_stats <- ggpubr::compare_means( mutation_rate ~ strain, 
                                            data= data_fig_5 ,
                                            group.by = c( "mutation" ),  
                                            p.adjust.method = "bonferroni", 
                                            ref.group = "mev-1",
                                            label = "p.signif", 
                                            method = "wilcox.test" ) %>% 
  #dplyr::filter(p.adj<0.05)%>% 
  dplyr::select(-p.format,-p.signif) %>% 
  dplyr::mutate(group_factor="mutation",
                group_factor_catogory=mutation,
                method="two-sided Wilcoxon test",
                padjustment="BF") %>% 
  dplyr::select(method,group_factor,group_factor_catogory,'.y.', group1,group2,p,padjustment,p.adj)


fig_5 <- ggpubr::ggboxplot(data_fig_5, x="strain",y="mutation_rate",outlier.shape = NA,
                            color="strain"  ) +
  geom_point( position = position_jitterdodge(jitter.width = 1) ,aes(color=strain), size=0.5, alpha=0.8)+
  facet_grid(.~mutation,scales="free")+
  theme_cust +
  theme(legend.position = "none")+
  scale_color_manual(values = c("orange","#007e2f","plum4","#721b3e") ) +
  labs(#x="Mutations",
    y="Mutation rate" ) + 
  ggpubr::stat_compare_means( 
    label = "p.signif",   
    ref.group = "mev-1",
    symnum.args = list(cutpoints = c(0, 0.00001, 0.0001, 0.001, 1), 
                       symbols = c("****","***", "**",  "ns")),
    size = 4,
    method = "wilcox.test") +
  scale_x_discrete("MA lines", labels=expression(N2,PB306, italic('mev-1')))


ggsave(fig_5, filename = paste( "../figures/fig_5.png",sep = ""), units = "mm",height = 70, width = 170)


##############################################
  
############# Figure  S1    #############
#              STR_13795                #
#########################################

str_exp_pxg <- data.table::fread("../processed_data/str_exp_pxg.tsv")

list2 <- str_exp_pxg %>% 
  dplyr::filter(pSTR %in% c("STR_13795", "STR_13083" )) %>% 
  dplyr::distinct(pSTR,transcript) %>% 
  dplyr::group_by(transcript) %>% 
  dplyr::count() %>% 
  dplyr::filter(n>1)

data_fig_S1 <- str_exp_pxg %>% 
  dplyr::filter(pSTR=="STR_13795") %>% 
  dplyr::mutate(overlap=ifelse(transcript %in% list2$transcript,"yes","no")) %>% 
  dplyr::mutate(transcript = paste(ext_gene, transcript, sep="\n")) 

data_fig_S1$genotype2<- factor(data_fig_S1$genotype,levels = unique(data_fig_S1$genotype))

fig_S1a <- ggplot(subset(data_fig_S1,ext_gene=="cls-2"),aes(x=as.factor(genotype2),y=expression ))+
  geom_jitter(shape=21,position=position_jitter(0.2),aes(fill=STRallele) ) +
  geom_boxplot( outlier.shape = NA,alpha=0.5) +
  theme_cust+
  theme(strip.text = ggplot2::element_text(size=12, vjust = 1,  color = "black",face = "italic"),
        legend.position = "none")+
  ggplot2::scale_fill_manual(values = c("REF"="orange","ALT"="blue"  ) )+
  xlab(paste("Length of STR_13795" )) + 
  ylab("Expression") +  
  facet_grid(.~transcript ,scales = "free") 

fig_S1b <- ggplot(subset(data_fig_S1,ext_gene!="cls-2"),aes(x=as.factor(genotype2),y=expression ))+
  geom_jitter(position=position_jitter(0.2),aes(fill=STRallele) ,shape=21) +
  geom_boxplot( outlier.shape = NA,alpha=0.5,aes(color=overlap)) +
  theme_cust+
  theme(strip.text = ggplot2::element_text(size=11, vjust = 1,  color = "black",face = "italic"),
        legend.position = "none") +
  ggplot2::scale_fill_manual(values = c("REF"="orange","ALT"="blue"  ) ) +
  ggplot2::scale_color_manual(values = c("yes"="red","no"="black"  ) ) +
  xlab(paste("Length of STR_13795" )) + 
  ylab("Expression") +  
  facet_wrap(.~transcript ,scales = "free",nrow=2) 


fig_S1a2 <- cowplot::plot_grid(fig_S1a,  NULL,  
                               label_size = 12, 
                               label_fontfamily="Helvetica",
                               nrow = 1)




fig_S1  <- cowplot::plot_grid(fig_S1a2,fig_S1b,
                              labels = c('A', 'B' ), 
                              label_size = 12, 
                              label_fontfamily="Helvetica",
                              rel_heights = c(0.5,1),
                              axis = "lr",
                              nrow = 2)


ggsave(fig_S1, filename = paste("../figures/supp_fig1_STR13795.png",sep = ""), units = "mm",height = 160, width = 170)



############# Figure  S2    #############
#              STR_13083                #
#########################################

data_fig_S2 <- str_exp_pxg %>% 
  dplyr::filter(pSTR=="STR_13083")%>% 
  dplyr::mutate(overlap=ifelse(transcript %in% list2$transcript,"yes","no"))%>% 
  dplyr::mutate(transcript = paste(ext_gene, transcript, sep="\n"))

data_fig_S2$genotype2<- factor(data_fig_S2$genotype,levels = unique(data_fig_S2$genotype))

fig_S2a <- ggplot(subset(data_fig_S2,ext_gene=="polq-1"),aes(x=as.factor(genotype2),y=expression ))+
  geom_jitter( position=position_jitter(0.2),aes(fill=STRallele) ,shape=21 ) +
  geom_boxplot( outlier.shape = NA,alpha=0.5,aes(color=overlap)) +
  theme_cust+
  theme(strip.text = ggplot2::element_text(size=12, vjust = 1,  color = "black",face = "italic"),
        legend.position = "none")+
  ggplot2::scale_fill_manual(values = c("REF"="orange","ALT"="blue"  ) )+
  ggplot2::scale_color_manual(values = c("yes"="red","no"="black"  ) ) +
  xlab(paste("Length of STR_13083" )) + 
  ylab("Expression") +  
  facet_grid(.~transcript ,scales = "free") 

fig_S2b <- ggplot(subset(data_fig_S2,ext_gene!="polq-1"),aes(x=as.factor(genotype2),y=expression ))+
  geom_jitter( position=position_jitter(0.2),aes(fill=STRallele )  ,shape=21) +
  geom_boxplot( outlier.shape = NA,alpha=0.5,aes(color=overlap)) +
  theme_cust+
  theme(strip.text = ggplot2::element_text(size=11, vjust = 1,  color = "black",face = "italic"),
        legend.position = "none")+
  ggplot2::scale_fill_manual(values = c("REF"="orange","ALT"="blue"  ) )+
  ggplot2::scale_color_manual(values = c("yes"="red","no"="black"  ) ) +
  xlab(paste("Length of STR_13083" )) + 
  ylab("Expression") +  
  facet_wrap(.~transcript ,scales = "free",nrow=2) 


fig_S2a2 <- cowplot::plot_grid(fig_S2a,  NULL,  NULL,  NULL,   
                               label_size = 12, 
                               label_fontfamily="Helvetica",
                               nrow = 1)


fig_S2  <- cowplot::plot_grid(fig_S2a2,fig_S2b,
                              labels = c('A', 'B' ), 
                              label_size = 12, 
                              label_fontfamily="Helvetica",
                              rel_heights = c(0.5,1),
                              axis = "lr",
                              nrow = 2)

ggsave(fig_S2, filename = paste("../figures/Supp_fig2_STR13083.png",sep = ""), units = "mm",height = 160, width = 170)



############# Figure  S3    #############
#           manha                       #
#########################################

###### fig_S3 ######

str_trait_st207_total <- data.table::fread("../processed_data/STR_mutation_trait.tsv")  


data_fig_S3a <-str_trait_st207_total %>% 
  dplyr::mutate(p1= ifelse(strain=="MY23","deletions",NA),
                p2= ifelse(strain=="MY23","insertions",NA),
                p3= ifelse(strain=="MY23","substitutions",NA)) 

fig_S3a <- ggplot(data_fig_S3a) + 
  geom_bar(stat='identity',aes( x=fct_reorder(strain, Total_mutation),y = Total_mutation ),fill="gray69",color="gray69")+
  theme_cust +
  theme(axis.text.x = element_blank(),
        legend.position =  c(0.35,0.7),
        # axis.title.y =  ggplot2::element_text(size=12,  color = "black",hjust =  -3),
        legend.title = element_blank(),
        legend.background = element_rect(  fill = NA ),
        axis.ticks.x=element_blank())+
  labs(x="Strains",
       y="STR variation",
       fill="Strains",
       color="Strains")+
  geom_point(aes(x=fct_reorder(strain, Total_mutation),y = deletion ),color="purple",size=0.1,alpha=0.8)+ 
  geom_point(aes(x=fct_reorder(strain, Total_mutation),y = insertion ),color="black",size=0.1,alpha=0.8)+ 
  geom_point(aes(x=fct_reorder(strain, Total_mutation),y = substitution ),color="red",size=0.1,alpha=0.8)+
  geom_text(  aes(label = p1, x="MY23", y=13 ),color="purple",size = 10*5/14 )+
  geom_text(  aes(label = p2, x="MY23", y=2 ),color="black" ,size = 10*5/14)+
  geom_text(  aes(label = p3, x="MY23", y=7 ),color="red" ,size = 10*5/14)



 




###### fig_S3b ######

data_fig_S3b <- data.table::fread("../processed_data/STR_nema_manha207.tsv") 
 
fig_S3b <- nema_manha_plot(data_fig_S3b)

###### fig_S3c ######
 

data_fig_S3c <- data.table::fread("../processed_data/STR_nema_manha_reg206.tsv") 
 
fig_S3c <- nema_manha_plot(data_fig_S3c)

  
### cow fig_S3 ####
fig_S3 <- cowplot::plot_grid(fig_S3a , fig_S3b, fig_S3c,
                            labels = c('A', 'B' ,'C' ), 
                            rel_heights = c(1.5,2,2),
                            label_size = 12, 
                            label_fontfamily="Helvetica",
                            axis = "lr",
                            nrow =  3)

ggsave(fig_S3, filename = paste( "../figures/Supp_fig3_manha.png",sep = ""), units = "mm",height = 180, width = 170)
 

 




#######