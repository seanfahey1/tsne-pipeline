#! /bin/bash

# Help function
Help()
{
   echo "Pipeline to produce t-SNE plots from path of .fasta formatted files."
   echo
   echo "Syntax: ./k-mer-pipeline.sh [path name] [.fasta file 1] ... [.fasta file n]"
   echo
   echo "options:"
   echo "   -h     Print this Help."
   echo
}

# check flags passed
while getopts ":h" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      *) echo "incorrect option(s): $*"
       exit 1 ;;
   esac
done

# TODO: set output path as required arg

# combine path to single file
name="$1"
combined_fasta="$name-combined.fasta"
basename="${combined_fasta%.*}"

for fasta_file in "${@:2}"
do
  echo "Merging $fasta_file into $combined_fasta..."
  touch "$combined_fasta"
  cat "$fasta_file" >> "$combined_fasta"
  echo "" >> "$combined_fasta"
done

# cd-hit combined file
echo "Running CD-Hit on $combined_fasta..."
cd_hit_output="$combined_fasta.output"
cd-hit \
  -i "$combined_fasta" \
  -o "$cd_hit_output"\
  -c 0.4 \
  -n 2 \
  -M 0 \
  -T 0 \
  -d 100 &&
echo "CD-Hit Output File: $cd_hit_output"

# blastn combined file
echo "Running blastn on $combined_fasta..."
blast_outfile="$basename.out"
truncated_blast_output="$basename-truncated.out"

# blast command for tsv formatted output
blastn -db nt -query "$combined_fasta" -outfmt 0 -out "$blast_outfile"

# get header row plus all first matches for truncated output
head -4 "$blast_outfile" | tail -1 > "$truncated_blast_output" && \
  grep -A 1 "hits found" "$blast_outfile" | grep -v -- "^--$" | grep -v "^#" >> "$truncated_blast_output"

echo "$truncated_blast_output"

# call PathObject python file
echo "Running python file..."
./main.py -n "$name" -i "${@:2}" -b "$blast_outfile" -c "$cd_hit_output.clstr"
