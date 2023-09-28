#!/bin/bash

# takes two position arguments
# $1: text file of paths to seeds file 1 per line
# $2: .fai file for the reference file used to make the seeds.

# Initialize the output and set input args
merged_seed_file="concat_seeds.bed"
output_file="seed_counts.bed"
bed_file_list="$1"
fai_file="$2"

# Remove the output file if it exists
[ -e "$output_file" ] && rm "$output_file"
[ -e "$merged_seed_file.tmp" ] && rm "$merged_seed_file.tmp"

# Loop through the list of BED files and merge them
while read -r bed_file; do
    if [ -e "$bed_file" ]; then
#        bedtools merge -i "$bed_file" >> "$merged_seed_file"
      cut -f 1-3 "$bed_file" >> "$merged_seed_file.tmp"
    else
        echo "File not found: $bed_file"
    fi
done < "$bed_file_list"

# Sort and merge the final output
sort -k1,1 -k2,2n "$merged_seed_file.tmp" > "$merged_seed_file"
bedtools genomecov -i "$merged_seed_file" -g "$fai_file" -bg > "$output_file"
#rm "$merged_seed_file" "$merged_seed_file.tmp"

echo "Seed counts: $output_file"
