import pandas as pd 
import bisect
import csv
import os 
import sys 


def find_ge_idx(a, x):
    """ find segment start"""
    return bisect.bisect_left(a, x)

def find_le_idx(a, x):
    """ segment end"""
    return bisect.bisect_right(a, x)

if __name__ == "__main__": 
    
    map_file = sys.argv[1] 
    input_file = sys.argv[2]
    output_file = sys.argv[3]
    i = int(sys.argv[4])
       
    ## get positions from map file 

    if os.path.isfile(map_file):         
        with open(map_file, "r") as map_input:
            map_reader = csv.reader(map_input, delimiter="\t")

            positions = [int(line[-1]) for line in map_reader]
            occurences = [0 for _ in positions]
            chrom = [i for _ in positions]

            positions.sort()
            

        ## get ibd per position in map file 
        if os.path.isfile(input_file): 
            with open(input_file) as ibd_input:
                ibd_reader = csv.reader(ibd_input, delimiter="\t")

                for line in ibd_reader:

                    start = find_ge_idx(positions, int(line[5]))
                    end = find_le_idx(positions, int(line[6]))
                    for pos_idx in range(start, end):
                        if pos_idx < len(positions):
                            occurences[pos_idx] += 1
                        else:
                            break
                            
            with open(output_file, 'w') as f: 
                ### write file where the total number of segments per position are output 
                values = zip(chrom, positions, occurences)
                w = csv.writer(f)
                w.writerows(values)



