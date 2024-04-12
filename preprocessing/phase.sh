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

###################### preprocess ##########################
# remove missing sites and add constant FID for PLINK 

    plink --vcf $output_dir"/chr"$i"merged.vcf.gz" \
        --const-fid 0  --allow-extra-chr  \
        --geno 0.01 --recode vcf  \
        --out $output_dir"/chr"$i"no_missing"
                
    bgzip -f $output_dir"/chr"$i"no_missing.vcf"
    tabix -f $output_dir"/chr"$i"no_missing.vcf.gz"

###################### run shapeit ##########################

    mkdir -p $output_dir"/shapeit"

    shapeit4 --input $output_dir"/chr"$i"no_missing.vcf.gz" \
        --map "shapeit4/maps/chr"$i".b38.gmap" --region $i --output \
        $output_dir"/shapeit/chr"$i"no_missing.vcf.gz" --thread 12

    tabix -f  $output_dir"/shapeit/chr"$i"no_missing.vcf.gz"


###################### post-process for iLASH #####################
# very annoyingly, iLASH requires a phased PED file for input, 
# but PLINK will not preserve the phasing information natively, 
# so i wrote a relatively simple script to keep this important info 

# also PLINK does not like multi alleles        
    bcftools view --max-alleles 2 \
       --exclude-types indels  $output_dir"/shapeit/chr"$i"no_missing.vcf.gz"  \
       -Oz -o  $output_dir"/chr"$i"no_multi.vcf.gz" 

    #used PLINK2 here because it's faster 
    plink2 --vcf $output_dir"/chr"$i"no_multi.vcf.gz"   \
        --const-fid 0  --allow-extra-chr \
        --maf 0.01 --geno 0.00 \
        --make-bed \
        --export vcf id-delim="." \
        --out $output_dir"/chr"$i"_filtered"    

    # use PLINK to generate the cM map for iLASH, 
    # important to get accurate IBD estimates 
    # as far as i can tell PLINK2 can't do this the same  
    # the ped file will be A,T,C,G and not numeric 
    # which is what iLASH expects 
    plink --bfile $output_dir"/chr"$i"_filtered"   \
        --cm-map shapeit4/maps/chr@.b38.gmap \
        --keep-allele-order \
        --allow-extra-chr --recode \
        --tab --out $output_dir"/chr"$i"_filtered" 

    # the output of this script is a ped and a map file that has the same name
    # this is a little unneccessary but it's just so that I have ped and map files with the same name
    cp $output_dir"/chr"$i"_filtered.map" $output_dir"/chr"$i"_converted.map"

    # my script just takes this as a way of intitalizing an array 
    num_positions=$(bcftools query -f '%POS\n' $output_dir"/chr"$i"_filtered.vcf" | wc -l)
    num_samples=$(bcftools query -l $output_dir"/chr"$i"_filtered.vcf" | wc -l)

    # take VCF file and use it to make a phased PED file 
    # correctly formatted for iLASH 
    python convert_to_plink.py \
        $output_dir"/chr"$i"_filtered.vcf" \
        $output_dir"/chr"$i"_converted" \
        $num_positions \
        $num_samples \
        $i



