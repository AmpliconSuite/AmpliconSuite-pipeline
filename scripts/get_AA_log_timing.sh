#!/bin/bash

# takes as input a text file of AA log files, one per line

input_file=$1
output_file="${input_file%.*}_timings.txt"

# Clear the output file if it exists
> "$output_file"

# Loop through each line in logs.txt
while IFS= read -r filename; do
	echo $filename
    if [ -f "$filename" ]; then
		last_line=`tail -n 1 ${filename} | cut -f 1 | cut -f 2 -d ' '`
		base_filename=$(basename "$filename")
		echo -e "${base_filename%.log}\t${last_line}" >> $output_file
	fi
done < "$input_file"
