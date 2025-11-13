#!/bin/bash
seqkit grep -n -r -p "chromosome" "${species_dir}/${species}2.1.fasta.gz" > "${species_dir}/${species}2.1.chromosome.fasta"
seqkit grep -n -r -p "contig|organelle" "${species_dir}/${species}2.1.fasta.gz" > "${species_dir}/${species}2.1.contig.fasta"
seqkit fx2tab --length --name --header-line "${species_dir}/${species}2.1.chromosome.fasta" > "${species_dir}/${species}2.1.sum.txt"
sum=$(awk '{sum+=$8} END {print sum}' "${species_dir}/${species}2.1.sum.txt")
y=$(zcat "${species_dir}/${species}2.2.fasta.gz" | grep ">" | wc -l)
z=$(cat "${species_dir}/${species}2.1.contig.fasta" | grep ">" | wc -l)
window_size=$(echo "($sum / ($y - $z))" | bc -l)
window_size=$(printf "%.0f" "$window_size")
seqkit sliding -s "$window_size" -W "$window_size" "${species_dir}/${species}2.1.chromosome.fasta" -o "${species_dir}/${species}2.1.chromosome.split.fasta"
cat "${species_dir}/${species}2.1.chromosome.split.fasta" "${species_dir}/${species}2.1.contig.fasta" > "${species_dir}/${species}2.1.split.fasta"
gzip "${species_dir}/${species}2.1.split.fasta"
rm "${species_dir}/${species}2.1.fasta.gz"
mv "${species_dir}/${species}2.1.split.fasta.gz" "${species_dir}/${species}2.1.fasta.gz"
rm "${species_dir}/${species}2.1.chromosome.fasta" "${species_dir}/${species}2.1.contig.fasta" "${species_dir}/${species}2.1.sum.txt" "${species_dir}/${species}2.1.chromosome.split.fasta"
unset sum
unset y
unset z
unset window_size
echo "${species}2.1 split successfully"

