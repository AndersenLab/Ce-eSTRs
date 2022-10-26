library(tidyverse)
library(lmtest)


args <- commandArgs(trailingOnly = TRUE)

# load trait 
exptrait <-  args[1]

#pair
transcript25849_Pos1Mb_pSTR <- data.table::fread(args[2])  
 
#exp
expression_207strains <- data.table::fread(args[3])

lgmtpm <- expression_207strains %>% 
  tidyr::gather(transcript,exp,-strain)

#str
 
SAMPLE_strGT_all <- data.table::fread(args[4])


SAMPLE_strLength_all <- data.table::fread(args[5])


 ### function for LRT test 
perm_lrtest_func <- function(str_geno){
 
  
  str_geno_perm <- str_geno %>% 
    dplyr::mutate(genotype_perm = base::sample(genotype, length(genotype)))
  
 ## true data  
  
  #fit full model
  model_full <- lm(exp ~ genotype, data = str_geno)
   #fit reduced model
  model_reduced <- lm(exp ~ 1, data = str_geno)
   #perform likelihood ratio test for differences in models
  lrt_test <- lmtest::lrtest(model_full, model_reduced)
  
   
  true_p <-  lrt_test$`Pr(>Chisq)`[2] 
  ## permuted data  
  
  
  #fit full model
  model_full_perm <- lm(exp ~ genotype_perm, data = str_geno_perm)
  
  #fit reduced model
  model_reduced_perm <- lm(exp ~ 1, data = str_geno_perm)
  
  #perform likelihood ratio test for differences in models
  lrt_test_perm <- lmtest::lrtest(model_full_perm, model_reduced_perm)
  

 # perm_p <- round(lrt_test_perm$`Pr(>Chisq)`[2], digits = 5)
  perm_p <-  lrt_test_perm$`Pr(>Chisq)`[2] 
  
  
  return( paste(true_p, perm_p,  sep = "_")) 
 
  
}


 
 ### input data for lrt test #####3
tx_pSTR_pair <- transcript25849_Pos1Mb_pSTR %>% dplyr::filter(transcript==exptrait)
  
tx_exp <- lgmtpm %>% dplyr::filter(transcript == exptrait)
  

  ####common geno ####
   
  SAMPLE_strGT <- SAMPLE_strGT_all %>% 
    dplyr::filter(pSTR %in% tx_pSTR_pair$pSTR) %>% 
    dplyr::left_join(tx_exp) %>% 
    na.omit() %>% 
    dplyr::group_by(pSTR) %>% 
    dplyr::add_count(name = "str_nst") %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(pSTR,genotype) %>% 
    dplyr::add_count(name="str_allele_n") %>% 
    dplyr::mutate(al=str_allele_n/str_nst) %>% 
    dplyr::distinct(pSTR,genotype,str_nst,str_allele_n,al) %>% 
    dplyr::mutate(major_allele=ifelse(al>0.05, "Major","Minor")) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(pSTR,major_allele) %>% 
    dplyr::add_count(name="major_allele_n") %>% 
    dplyr::filter(major_allele=="Major" & major_allele_n>=2) %>% 
    dplyr::ungroup() %>% 
    dplyr::distinct(pSTR,genotype) %>% 
    dplyr::left_join(SAMPLE_strGT_all) %>% 
    dplyr::left_join(tx_exp) %>% 
    na.omit()
  
 
# test genotype as factoral data
  if(nrow(SAMPLE_strGT)>0){
   
  str_FACgeno_lrt <- SAMPLE_strGT %>% 
     dplyr::mutate(genotype=as.factor(genotype)) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(transcript, pSTR) %>% 
    dplyr::do(data.frame(lrt_pvalue=perm_lrtest_func(.))) %>% 
    tidyr::separate(lrt_pvalue, into = c("real_lrt_p","perm_lrt_p" ),sep="_") %>% 
    na.omit()   %>% 
    dplyr::mutate(str_genotype="STR_genotype")
  
  }else{
    
     str_FACgeno_lrt <- data.frame()
  }
  
### common length ####
  
  SAMPLE_strLength <- SAMPLE_strLength_all %>% 
    dplyr::filter(pSTR %in% tx_pSTR_pair$pSTR) %>% 
    dplyr::left_join(tx_exp) %>% 
    na.omit() %>% 
    dplyr::group_by(pSTR) %>% 
    dplyr::add_count(name = "str_nst") %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(pSTR,genotype) %>% 
    dplyr::add_count(name="str_allele_n") %>% 
    dplyr::mutate(al=str_allele_n/str_nst) %>% 
    dplyr::distinct(pSTR,genotype,str_nst,str_allele_n,al) %>% 
    dplyr::mutate(major_allele=ifelse(al>0.05, "Major","Minor")) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(pSTR,major_allele) %>% 
    dplyr::add_count(name="major_allele_n") %>% 
    dplyr::filter(major_allele=="Major" & major_allele_n>=2) %>% 
    dplyr::ungroup() %>% 
    dplyr::distinct(pSTR,genotype) %>% 
    dplyr::left_join(SAMPLE_strLength_all)%>% 
    dplyr::left_join(tx_exp) %>% 
    na.omit()
  
  # test length as factoral data
  
  if(nrow(SAMPLE_strLength) >0 ){
  
  
  str_FACgenoL_lrt <- SAMPLE_strLength %>% 
    dplyr::mutate(genotype=as.factor(genotype)) %>% 
    dplyr::group_by(transcript, pSTR) %>% 
    dplyr::do(data.frame(lrt_pvalue=perm_lrtest_func(.))) %>% 
    tidyr::separate(lrt_pvalue, into = c("real_lrt_p","perm_lrt_p" ),sep="_") %>% 
    na.omit() %>% 
    dplyr::mutate(str_genotype="STR_length")
  
  }else{
    
    str_FACgenoL_lrt <- data.frame()
  }
  
  
  
  ### summarize ####
  str_lrt <- dplyr::bind_rows(str_FACgeno_lrt, str_FACgenoL_lrt)
  
  if( nrow(str_lrt) > 0 ){

write.table(str_lrt, paste(exptrait, "LrtD_eSTRs.tsv", sep = "_"), sep = "\t", row.names = F,  quote = FALSE)

  }
