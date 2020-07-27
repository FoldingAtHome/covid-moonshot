#!/bin/bash

# Fix CSV files so that SMILES and then CID occur first

for prefix in "primary_amine_enumeration_for_chodera_lab_FEP" "boronic_ester_enumeration_for_chodera_lab_FEP"; do
  echo "Generating $prefix-permuted.csv"
  awk -F',' '{ print $4 "," $1 "," $2 "," $3 }' $prefix.csv > $prefix-permuted.csv
  echo "Generating $prefix.pdf"
  mols2pdf.py -in $prefix-permuted.csv -out $prefix.pdf -cols 6 -rows 3
done
