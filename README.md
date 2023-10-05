# IBD Analysis

## Description

This is code to accompany the manuscript "Health care utilization of fine-scale identity by descent clusters in a Los Angeles biobank."

The goal of this research is to identify clusters of patients based on their genetic relatedness. We then examine why and how these clusters come to the hospital, including diagnoses and what zip codes they visit doctors in.

This repository contains code for identifying these clusters, which was done through identity-by-descent mapping along with Louvain clustering. It also contains code for some population genetics analyses, along with calculating the cluster-disease associations.

For a visualization of our results, please see our website at www.ibd.la and the accompanying github repository at https://github.com/misingnoglic/ibd.la

<b>Note, that this repository is under active construction. Check back in for updates!</b>

## Repository Structure

`admixture`: code for generating admixture plots with the software package [SCOPE](https://github.com/sriramlab/SCOPE)

`associations`: code for calculating the phecode, specialty, and zipcode associations

`fst`: script for calculating Hudson's FST from summary statistic data

`IBD-calling`: code for generating IBD calls with [iLASH](https://github.com/roohy/iLASH) and identifying clusters with the Louvain algorithm

## External Software and Dependencies

Data processing pipelines were written in bash. Analysis and visualization was done in jupyter noteboks.

- Basic python 3 packages needed include `pandas`, `seaborn`, `networkx`, `statsmodels.api`, `jupyter`, and `numpy.`

- [iLASH](https://github.com/roohy/iLASH)

- [SCOPE](https://github.com/sriramlab/SCOPE)

- [PLINK](https://zzz.bwh.harvard.edu/plink/)

- [bedtools](https://bedtools.readthedocs.io/en/latest/)

- [bcftools](https://samtools.github.io/bcftools/bcftools.html)


## Contact

This code is a work in progress and will be updated to improve usability. If you have any questions please contact christa@g.ucla.edu

## Citation

Caggiano, C., Boudaie, A., Shemirani, R., Mefford, J., Petter, E., Chiu, A., Ercelen, D., He, R., Tward, D., Paul, K.C., Chang, T.S., Pasaniuc, B., Kenny, E.E., Shortt, J.A., Gignoux, C.R., Balliu, B., Arboleda, V.A., Belbin, G., Zaitlen, N., 2023. Disease risk and healthcare utilization among ancestrally diverse groups in the Los Angeles region. Nat Med 29, 1845–1856. https://doi.org/10.1038/s41591-023-02425-1
