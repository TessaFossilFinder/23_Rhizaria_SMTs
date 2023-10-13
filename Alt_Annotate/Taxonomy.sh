#!/bin/bash

# Input file containing taxids (tmp4)
input_file="tmp4"
output_file="NCBI_Taxonomy.csv"

# Fetch taxonomic information for each taxid in tmp4
while IFS=$'\t' read -r accession taxid; do
    taxonomy=$(efetch -db taxonomy -id "$taxid" -format xml | \
               xtract -pattern Taxon -first TaxId -element TaxId,Kingdom,Clade,Phylum,Class,Order,Family,Genus,Species  -block "*/Taxon" -unless Rank -equals "no rank" -tab ",")
    echo "$taxid,$taxonomy" >> "$output_file"
done < "$input_file"

echo "Taxonomic information fetched and stored in $output_file."
