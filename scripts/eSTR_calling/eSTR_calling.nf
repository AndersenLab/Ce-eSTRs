#!/usr/bin/env nextflow



date = new Date().format( 'yyyyMMdd' )
params.out = "eSTR_${date}"



 
/*
~ ~ ~ > *  STR_transcript pair path file
*/


params.fqs = "$baseDir/transcript25849_Pos1Mb_pSTR_test.csv"




File fq_file = new File(params.fqs)

fq_handle = fq_file.getAbsolutePath()
  

Channel
     .fromPath(params.fqs)
     .splitCsv(header: true)
     .map{ row -> "${row.transcript}" }
     .unique()
     .into { tx; 
            tx2}


/*
~ ~ ~ > *  Expression data
*/


params.exp_files = "$baseDir/bin/expression_207strains.tsv"


File exp = new File("${params.exp_files}")
exp_handle = exp.getAbsolutePath()





/*
~ ~ ~ > *  STR genotype
*/


params.strGT_files = "$baseDir/bin/STR_genotype_540strains.tsv"


File strGT = new File("${params.strGT_files}")
strGT_handle = strGT.getAbsolutePath()




/*
~ ~ ~ > *  STR length
*/


params.strL_files = "$baseDir/bin/STR_length_540strains.tsv"

File strL = new File("${params.strL_files}")
strL_handle = strL.getAbsolutePath()


 


/*
~ ~ ~ > *  call eSTRs
*/


  process eSTR {


     
    tag "${transcript}"

    cpus 4
    memory '32 GB'
    
    input:

      val(transcript) from tx
      
    output:
      
      set val(transcript), file("${transcript}_LrtD_eSTRs.tsv") optional true into all_eSTRs 
 


    """

    Rscript --vanilla `which LRT_nf_perm.R` ${transcript} ${fq_handle} ${exp_handle} ${strGT_handle} ${strL_handle}
       
    """

}


/*
~ ~ ~ > *  summarize eSTRs
*/

  process eSTR_summ {

    publishDir "${params.out}", mode: 'copy', pattern: "eSTR_lrt.tsv"
 
    cpus 1
    memory '32 GB'
    
    input:

      file('*') from all_eSTRs.collect()
      
    output:
      
      file("eSTR_lrt.tsv")


    """
 

      cat *_LrtD_eSTRs.tsv | grep "transcript" | sort | uniq  > eSTR_lrt.tsv

      cat *_LrtD_eSTRs.tsv | grep -v "transcript"  >> eSTR_lrt.tsv


 

    """
    

}







