import sys
import os
import numpy as np
import gc
import vcfpy
import csv 


def get_sample_headers(vcf, num_individuals):
    """
    take the VCF file and get the basic information about each 
    individual that would be the beginning of PED file 
    i am just filling it in with NA since it isn't needed for IBD calling 
    but is needed for the file to be formatted correctly 
    """
    
    sample_ids = vcf.header.samples.names
        
    family_ids = ["NA"] * num_individuals
    individual_ids = [" " + s for s in sample_ids]
    maternal_ids = ["NA"] * num_individuals
    paternal_ids = ["NA"] * num_individuals
    sex = ["NA"] * num_individuals
    phenotypes = ["NA"] * num_individuals
    
    individual_headers = np.array([family_ids, individual_ids, maternal_ids, paternal_ids, sex, phenotypes]).T
    
    return individual_headers

def convert_gt(genotype):
    """
    ped has a different convention for reporting alleles 
    generally major allele is 1 and minor allele is 2 
    in VCF 2 is a second allele that shouldn't be included for iLASH so 
    i make it 0 for missing instead 
    """
    genotype_map = {"1":"2", "0":"1", "2":"0"}
    return genotype_map[genotype]


if __name__ == "__main__":
    
    infile = sys.argv[1]
    outfile = sys.argv[2]
    num_sites = int(sys.argv[3])
    num_individuals = int(sys.argv[4])
    chrom = int(sys.argv[5]) 
    
    # ped format is number of indivduals x number of sites 
    # except that each site is represented as two integers 
    #(one for maternal, one for paternal), so in practice 
    # it is the (number of sites)x2 
    genotypes = np.empty([num_individuals, num_sites*2], dtype="int")
    rsids = []
    positions = []
    
    vcf = vcfpy.Reader.from_path(infile)
            
    print("reading vcf") 
    print()
    
    genotype_idx = 0 
    for record in vcf:
        
        rsids.append(record.ID[0])
        positions.append(int(record.POS))
    
        individual_idx = 0 

        for call in record.calls:
            
            if "|" in call.data["GT"]: 
                genotype0, genotype1 = call.data["GT"].split("|")  # phased SNP and preserve order 
                
            else: 
                genotype0, genotype1 = call.data["GT"].split("/")
            
            genotypes[individual_idx, genotype_idx] = convert_gt(genotype0)
            genotypes[individual_idx, genotype_idx+1] = convert_gt(genotype1)
                        
            individual_idx += 1 
            
        genotype_idx += 2 

    individual_headers = get_sample_headers(vcf, num_individuals)

    
    print("writing " + ".ped") 
    print()
    np.savetxt(outfile + ".ped", genotypes, delimiter=" ", fmt='%s') 
    np.savetxt(outfile + ".header", individual_headers, delimiter=" ", fmt='%s') 

    