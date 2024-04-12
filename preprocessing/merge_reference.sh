#!/bin/sh

###################### data paths ######################
qc_dir="azure/june21_atlas_update/final_qc"
output_dir="azure/previous_reference"
simons_path="/opt/genomics/IPHatlasreleases/SimonsDiversityData/hg38/annotated"
hgdp_path="/opt/genomics/IPHatlasreleases/HGDP/annotated"
thousandG_path="/opt/genomics/IPHatlasreleases/1000_genomes/phased_vcfs/annotated_rsid"
atlas_path="final_qc"

for ((i=22; i>=1; i--))
do

###################### subset chromosomes ##########################
# make ATLAS data individual chromosomes for merging and IBD calling 

    plink --bfile $qc_dir"/final_qc_chr"$i --recode vcf --out $qc_dir"/final_qc_chr"$i

    bgzip -f $qc_dir"/final_qc_chr"$i".vcf"
    tabix -f $qc_dir"/final_qc_chr"$i".vcf.gz"


#################### put to same dbsnp build #######################
# i want to normalize all the reference data and ATLAS to the same DbSNP 
# this helps not to miss any SNPs when merging together. In IBD, more SNPS = better  
# at some point i did this for the reference data as well but you don't need to do again 
# if you use the sme DbSNP build, otherwise add this step for reference data 

    bcftools annotate \
         -a dbsnp/All_20180418.vcf.gz  \
         -c ID $qc_dir"/final_qc_chr"$i".vcf.gz" \
         -O z \
         -o $qc_dir"/final_qc_chr"$i"_annotated.vcf.gz"

    gunzip -f $qc_dir"/final_qc_chr"$i"_annotated.vcf.gz" # unzip for next step 

################## flip strand and norm ############################
# in my version of ATLAS, there was some strand issues so I normalize to the the same genome build
# as the reference data. Unfortunately requires manually adding "chr" in front of numeric chromosomes 

    awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' $qc_dir"/final_qc_chr"$i".vcf" > $qc_dir"/final_qc_chr"$i"_annotated_chr.vcf"
     
    bgzip -f  $qc_dir"/final_qc_chr"$i"_annotated_chr.vcf"
    tabix -f  $qc_dir"/final_qc_chr"$i"_annotated_chr.vcf.gz"

    bcftools norm \
        $qc_dir"/final_qc_chr"$i"_annotated_chr.vcf.gz" \
        -m-any --check-ref s -f genome_references/hg38.fa \
        -Ov -o $"/final_qc_chr"$i"_fixref_chr.vcf"

    awk '{gsub(/^chr/,""); print}' $qc_dir"/final_qc_chr"$i"_fixref_chr.vcf" > $qc_dir"/final_qc_chr"$i"_fixref2.vcf"

    bgzip -f  $qc_dir"/final_qc_chr"$i"_fixref2.vcf" # shitty name sorry 
    tabix -f  $qc_dir"/final_qc_chr"$i"_fixref2.vcf.gz"

################## subset reference data  ############################

    mkdir -p $output_dir"/reference"
       
    # take final ATLAS sites to subset reference data to 
    awk '{print $1,"\t", $2,"\t",$2}' \
    $qc_dir"/final_qc_chr"$i"_annotated.vcf" > $output_dir"/reference/chrom"$i"_positions.txt"
        
############################# 1000G 
    # the hack-y way i figured out how to very quickly subset 
    # reference data to a bedfile is to use tabix -R 
    # it takes a lot less time to first subset and then merge 
    # however, for some reason newer versions of tabix don't have this functionality 
    # I am using tabix 1.10.2  

    tabix -R $output_dir"/reference/chrom"$i"_positions.txt" \
        $thousandG_path/"chr"$i"_1000G_phased_hg38_annotated_redone.vcf.gz" > \
        $output_dir"/reference/chrom"$i"_1000g.txt"
    
    # add the header back to vcf 
    tabix -H  $thousandG_path/"chr"$i"_1000G_phased_hg38_annotated_redone.vcf.gz" > \
        $output_dir"/reference/chrom"$i"_header_1000g" 
        
    cat $output_dir"/reference/chrom"$i"_header_1000g" $output_dir"/reference/chrom"$i"_1000g.txt" > \
        $output_dir"/reference/chrom"$i"_1000g_header.vcf"
        
    bgzip -f $output_dir"/reference/chrom"$i"_1000g_header.vcf"
    tabix -f $output_dir"/reference/chrom"$i"_1000g_header.vcf.gz"
          
############################# HGDP 
    tabix -R $output_dir"/reference/chrom"$i"_positions.txt" \
        $output_dir"/reference/hgdp_only_samples_chr"$i".statphase_annotated_redone2.vcf.gz" \
        > $output_dir"/reference/chrom"$i"_hgdp.txt"
        
    tabix -H  $output_dir"/reference/hgdp_only_samples_chr"$i".statphase_annotated_redone2.vcf.gz" > \
        $output_dir"/reference/chrom"$i"_header_hgdp" 
        
    cat  $output_dir"/reference/chrom"$i"_header_hgdp"  $output_dir"/reference/chrom"$i"_hgdp.txt" > \
        $output_dir"/reference/chrom"$i"_hgdp_header.vcf"
        
    bgzip -f $output_dir"/reference/chrom"$i"_hgdp_header.vcf"
    tabix -f $output_dir"/reference/chrom"$i"_hgdp_header.vcf.gz"
          
############################# SGDP

    tabix -R $output_dir"/reference/chrom"$i"_positions.txt" \
        $simons_path/"sgdp_phased_chr"$i"_hg38_annotated_redone.vcf.gz" > \
        $output_dir"/reference/chrom"$i"_sgdp.txt"
        
    tabix -H  $simons_path/"sgdp_phased_chr"$i"_hg38_annotated_redone.vcf.gz" > \
        $output_dir"/reference/chrom"$i"_sgdp_header" 
        
    cat $output_dir"/reference/chrom"$i"_sgdp_header" $output_dir"/reference/chrom"$i"_sgdp.txt" > \
        $output_dir"/reference/chrom"$i"_sgdp_header.vcf"
    
    bgzip -f $output_dir"/reference/chrom"$i"_sgdp_header.vcf"
    tabix -f $output_dir"/reference/chrom"$i"_sgdp_header.vcf.gz"

################## merge finally!  ############################
# on lapgnomap01, you can use a lot of cores if you need - > 70 

    bcftools merge $qc_dir"/final_qc_chr"$i"_fixref2.vcf.gz" \
        $output_dir"/reference/chrom"$i"_1000g_header.vcf.gz" \
        $output_dir"/reference/chrom"$i"_hgdp_header.vcf.gz" \
        $output_dir"/reference/chrom"$i"_sgdp_header.vcf.gz" \
        --threads 32 -Oz -o $output_dir"/chr"$i"merged.vcf.gz"

done 