#!/bin/bash
# create input files for moessbauer benchmark defining equal inputs:
# standard.inp

# can be batch submitted with something like
# for inp in cal/*/*inp; do
#   (
#     cd $(dirname $i)
#     suborca $(basename $i .inp)
#   )
# done


charge_mult_line=(
  "" # empty because bash arrays start with 0
  "* xyzfile -2 5 geom.xyz" 
  "* xyzfile  -4 1 geom.xyz"
  "* xyzfile  -4 5 geom.xyz"
  "* xyzfile  -1 5 geom.xyz"
  "* xyzfile  -1 5 geom.xyz"
  "* xyzfile  1 2 geom.xyz"
  "* xyzfile  -1 6 geom.xyz"
  "* xyzfile  -3 6 geom.xyz"
  "* xyzfile  -3 2 geom.xyz"
  "* xyzfile  -3 6 geom.xyz"
  "* xyzfile 3 6  geom.xyz"
  "* xyzfile -2 4 geom.xyz"
  "* xyzfile  1 2 geom.xyz"
  "* xyzfile -1 6 geom.xyz"
  "* xyzfile 0 6 geom.xyz"
  "* xyzfile -1 5 geom.xyz"
  "* xyzfile 1 4  geom.xyz"
  "* xyzfile 2 3 geom.xyz"
  "* xyzfile -2 3 geom.xyz"
  "* xyzfile 0 1 geom.xyz"
  "* xyzfile 1 1 geom.xyz"
  "* xyzfile 0 2 geom.xyz"
  "* xyzfile 0 1 geom.xyz"
  "* xyzfile 0 1 geom.xyz"
)

top_input=template.inp
xyz_folder=xyz/
cal_folder=cal/
bot_input="
%eprnmr
nuclei = all Fe {fgrad, rho}
end
"

for i in $(seq -w 24); do 
  folder="$cal_folder$i/" # folder for each calculation
  input="${folder}moess.inp" # input file

  mkdir -p "$folder"
  cp "$xyz_folder$i.xyz" "$folder/geom.xyz" # copy xyz file

  i_=$(echo $i|bc) # $i will be interpreted as octal
  {
    cat  "$top_input" # 1. general input
    echo
    echo "${charge_mult_line[$i_]}" # specific charge and mult
    echo "$bot_input" # eprnmr block
  } > "$input"
done
