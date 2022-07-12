import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

all_sites = pd.read_csv(input_file, header=None)
all_sites.columns = ["chrom", "start", "count"]

all_sites = all_sites.sort_values(by=["chrom", "start"])

mean = all_sites["count"].mean()
sd = all_sites["count"].std()

above_value = mean + (1*sd)
below_value = mean - (1*sd)

outliers = all_sites[(all_sites["count"] >= above_value) | (all_sites["count"] <= below_value)]

print(len(outliers))
outliers = outliers[["chrom", "start"]]
outliers["end"] = outliers["start"].astype(int) + 1

blacklist = pd.read_csv("hg38-blacklist.v2.bed", header=None, sep="\t" )

blacklist = blacklist[[0, 1, 2]]
blacklist.columns = ["chrom", "start", "end"]

hla = pd.read_csv("hla/hla_genes_regions_collapsed.bed", header=None, sep="\t")
hla.columns = ["chrom", "start", "end"]

gaps = pd.read_csv("gap.bed", header=None, sep="\t")
gaps.columns = ["chrom", "start", "end"]

all_outliers = pd.concat([outliers, blacklist, hla, gaps])
all_outliers = outliers

all_outliers = all_outliers.sort_values(by=["chrom", "start", "end"])
all_outliers.to_csv(output_file, header=False, index=False, sep="\t")
