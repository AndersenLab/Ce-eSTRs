nema_manha_plot <- function(processed_mapping_data){
  
  processed_mapping_data$algorithm2<- factor(processed_mapping_data$algorithm,levels = c("LOCO","INBRED"))
  
  pmax<-max(processed_mapping_data$log10p)+0.5
  
  
  test <- processed_mapping_data %>% 
    dplyr::select(BF,EIGEN) %>% 
    dplyr::distinct() %>% 
    tidyr::gather(name,value)
  
  peak_markers <- na.omit(processed_mapping_data) %>% 
    dplyr::distinct(CHROM,POS,log10p,sig)
  
  fig_manha <- ggplot() + 
    geom_point(data = processed_mapping_data, 
               mapping = aes(x = POS/1e6, 
                             y = log10p,
                             colour = sig ),size=0.5) +
    scale_colour_manual(  values = c(  "BF" = "red", "EIGEN" = "hotpink3","NONSIG" = "black"))   + 
    scale_x_continuous(expand = c(0, 0), breaks = c(0, 5, 10, 15, 20),limits = c(0,NA)) +
    geom_hline(data = test, aes(yintercept = value, linetype = name),  color = "gray", 
               alpha = .75,  
               size = 1) + 
    scale_linetype_manual(values = c("BF" = 1, "EIGEN" = 2 )) + 
    labs(x = "Genomic position (Mb)",
         y = expression(-log[10](italic(p)))) +
    theme_cust +
    theme(legend.position = "none" , 
          axis.text.x = element_blank()) + 
    facet_grid(algorithm2 ~ CHROM, scales = "free_x" ) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0,pmax)) 
  
  return(fig_manha)
  
}

