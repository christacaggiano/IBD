import csv 
import pandas as pd
import sys
import pickle as pkl 
from bisect import * 

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def write_output(output_file, pair_dict, output_file_type): 
    
    if output_file_type == "pkl": 
        pkl.dump(pair_dict, open(f"{output_file}.pkl", "wb"))

    else: 
        with open(f"{output_file}.csv", "w") as csv_file:  
            writer = csv.writer(csv_file)

            for key, value in pair_dict.items():
                writer.writerow([key[0], key[1], value[0], value[1]])

    
def strip_char(value): 
        return value[:-2]

def find_ge(a, x):
    """Find leftmost item greater than or equal to x"""
    return bisect_left(a, x)


def find_le(a, x):
    """Find rightmost value less than or equal to x"""
    return bisect_right(a, x) - 1


def check_outlier(outlier_starts, outlier_ends, start, end): 
    """check if it is an outlier against a sorted list, using binary search """
    
    if len(outlier_starts) > 0 and len(outlier_ends) > 0: 
    
        if not outlier_starts is None: 
            idx = find_ge(outlier_starts, start)
            if idx != len(outlier_starts):
                if (outlier_starts[idx] <= start <= outlier_ends[idx]) or (outlier_starts[idx] <= end <= outlier_ends[idx]) :
                    return True

            idx = find_le(outlier_starts, start)
            if (outlier_starts[idx] <= start <= outlier_ends[idx]) or (outlier_starts[idx] <= end <= outlier_ends[idx]):
                return True

            idx = find_le(outlier_ends, end)
            if (outlier_starts[idx] <= start <= outlier_ends[idx]) or (outlier_starts[idx] <= end <= outlier_ends[idx]):
                return True
            return False

            idx = find_ge(outlier_ends, end)
            if idx != len(outlier_starts):
                if (outlier_starts[idx] <= start <= outlier_ends[idx]) or (outlier_starts[idx] <= end <= outlier_ends[idx]):
                    return True
                
    return False


def add_pair(line, outlier_starts, outlier_ends, pair_dict):
    """add all the IBD shared between two individuals on this chromosome """
    
    ind1 = strip_char(line[1])
    ind2 = strip_char(line[3]) 

    start = int(line[5])
    end = int(line[6])

    ibd = float(line[-2])
    prob = str(line[-1])
    
    if not(prob == "-nan" or prob == "nan"):

        if ind1 != ind2 and not check_outlier(outlier_starts, outlier_ends, start, end):

            p = [ind1, ind2]
            p.sort()

            ind1 = p[0]
            ind2 = p[1]


            if (ind1, ind2) in pair_dict: 

                pair_dict[(ind1, ind2)][0] += ibd
                pair_dict[(ind1, ind2)][1] += 1

            else: 
                pair_dict[(ind1, ind2)] = [ibd, 1]

            
if __name__ == "__main__":
    
    chrom = int(sys.argv[1])
    
    input_file = str(sys.argv[2])  

    output_file = sys.argv[3]
    outlier_file = sys.argv[4]
    output_file_type = sys.argv[5]

    pair_dict = {}

    outliers = pd.read_csv(outlier_file, delimiter="\t", header=None)
    chrom_outliers = outliers[outliers[0] == chrom]

    outlier_starts = list(outliers[1])
    outlier_ends = list(outliers[2])

    with open(input_file, "r") as input: 
        input_reader = csv.reader(input, delimiter="\t") 
        
        for line in input_reader: 
            add_pair(line, outlier_starts, outlier_ends, pair_dict) 

    write_output(output_file, pair_dict, output_file_type) 
   
