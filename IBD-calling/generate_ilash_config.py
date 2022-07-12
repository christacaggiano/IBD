import sys

chrom = int(sys.argv[1])
map_file = sys.argv[2]
ped_file = sys.argv[3]
output_file = sys.argv[4]

file_contents = f"""map {map_file}

ped {ped_file}

output {output_file}

slice_size 350

step_size 350

perm_count 20

shingle_size 15

shingle_overlap 0

bucket_count 5

max_thread 20

match_threshold 0.99

interest_threshold 0.70

min_length 2.9

auto_slice 1

slice_length 2.9

cm_overlap 1

minhash_threshold 55
"""

config_file_name = f"ilash_config_chr{chrom}"

text_file = open(config_file_name, "wt")
n = text_file.write(file_contents)
text_file.close()
