import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
import sys 



###### input files 
input_file = sys.argv[1] # per chromosome IBD file
output_file = sys.argv[2] # file with ibd segments overlapping removed 
blacklist_file = sys.argv[3] # path to bed file of regions you want to remove (i.e. centromeres)
hla_file = sys.argv[4] # path to bed file of HLA regions to remove 
gap_file = sys.argv[5] # path to bed file of gaps in hg38 genome  

all_sites = pd.read_csv(input_file, header=None)
all_sites.columns = ["chrom", "start", "count"]
all_sites = all_sites.sort_values(by=["chrom", "start"])


#### calculate the mean and standard deviation of IBD depth 
mean = all_sites["count"].mean()
sd = all_sites["count"].std()

above_value = mean + (3*sd) 
below_value = mean - (3*sd) 


### define outliers on IBD coverage
outliers = all_sites[(all_sites["count"] >= above_value) | (all_sites["count"] <= below_value)]


### format outliers correctly  
outliers = outliers[["chrom", "start"]]
outliers["end"] = outliers["start"].astype(int) + 1 
outliers["chrom"] = outliers["chrom"].astype(str)


### read in and format the regions to remove correctly 
blacklist = pd.read_csv(blacklist_file, header=None, sep="\t" )
blacklist = blacklist[[0, 1, 2]]
blacklist.columns = ["chrom", "start", "end"]
blacklist["chrom"] = blacklist["chrom"].str[3:]

hla = pd.read_csv(hla_file, header=None, sep="\t")
hla.columns = ["chrom", "start", "end"]
hla["chrom"] = hla["chrom"].str[3:]

gaps = pd.read_csv(gap_file, header=None, sep="\t")
gaps.columns = ["chrom", "start", "end"]
gaps["chrom"] = gaps["chrom"].str[3:]


### group outliers and regions to remove together 
all_outliers = pd.concat([outliers, blacklist, hla, gaps])

### output a bed file 
all_outliers = all_outliers.sort_values(by=["chrom", "start", "end"])
all_outliers.to_csv(output_file, header=False, index=False, sep="\t")