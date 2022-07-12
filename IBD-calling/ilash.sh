
ilash_output="ilash_output/all_reference"
mkdir -p $ilash_output
echo $ilash_output


for ((i=22; i>=1; i--))
  do
    echo $i

    file="merged_with_reference/chr"$i"_converted"

    mv $file".ped" $file"_temp"
    paste -d ' ' $file".header" $file"_temp" > $file".ped"

    awk '$2 = $2 FS "NA"' $file".ped" > temp
    mv temp $file".ped"

    python generate_ilash_config.py $i $file".map" $file".ped" $ilash_output"/chr"$i".match"
    .iLASH/build/ilash "ilash_config_chr"$i

    python calculate_ibd_density.py $file".map" $ilash_output"/chr"$i".match" $ilash_output"/chr"$i"_depth.csv" $i

done

rm $ilash_output"/all_depth.csv"
for file in $ilash_output"/chr"*"_depth.csv";
  do
  cat $file >> $ilash_output"/all_depth.csv";
done

python remove_outliers.py $ilash_output"/all_depth.csv" $ilash_output"/outliers.csv"

rm $ilash_output"/outliers_collapsed.csv"
bedtools merge -i $ilash_output"/outliers.csv" > $ilash_output"/outliers_collapsed.csv"


for ((i=1; i>=1; i--))
  do
  # echo $i

    ./combine_split_file.sh $i $ilash_output"/outliers_collapsed.csv" $ilash_output

done
