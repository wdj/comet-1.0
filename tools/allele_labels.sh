#!/bin/bash
#==============================================================================
# Extract the allele labels used in each line of the SNP text file.
# Recognizes A, C, G, T
# To conform to assumptions of other codes used here, for each line
# the labels are output in (ascending) alphabetical order.
#==============================================================================

function main
{
  if [ "$2" = "" ] ; then
    echo "${0##*/}: extract per-line allele labels from SNP text file"
    echo "Usage: ${0##*/} <snp_text_file> <allele_label_file>"
    exit
  fi

  local infile="$1"
  local outfile="$2"
  #local infile="28M_permuted/28M.txt"
  #local outfile="28M_permuted/28M_allele_labels.txt"

  cut -f5- < "$infile" \
    | sed -e 's/^[T0\t]*$/tt/' -e 's/^[A0\t]*$/aa/' -e 's/^[G0\t]*$/gg/' -e 's/^[C0\t]*$/cc/' \
          -e 's/^[CT0\t]*$/CT/' -e 's/^[AG0\t]*$/AG/' -e 's/^[AC0\t]*$/AC/' \
          -e 's/^[CG0\t]*$/CG/' -e 's/^[AT0\t]*$/AT/' -e 's/^[GT0\t]*$/GT/' \
    | tr 'a-z' 'A-Z' \
    > "$outfile"
}

#------------------------------------------------------------------------------

main $@

#==============================================================================
