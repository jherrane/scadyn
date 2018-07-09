#!/bin/bash
base="./scadyn"
geom_fold="Geometries"
T_fold="T_socpol"

# Sizes (in um)
# 0.0500 0.0598 0.0715 0.0855 0.1022 0.1223 0.1462 0.17480.2090 0.2500
declare -a shapes=("ell" "ob" "pro" "sph")
declare -a compositions=("P-aC")

BEGIN=1
END=15

for is in $(seq 0 3); do
   for iw in $(seq 1 10); do
      for i in $(seq $BEGIN $END); do
         mesh=" --mesh "$geom_fold"/${shapes[$is]}-P-aC-$i.h5"
         let "j = i + is*END"

         T=" -T $T_fold/T-$j.h5"
         out=" --out -$j"
         singleT=" -S 1"
         which=" -w $iw"
         
         comm="$base$mesh$T$out$singleT$which"
         repl="sed -i '9s|.*|"$comm"|' script.sh"
         
         echo "$comm"
#         eval "$repl"
#         eval "sbatch script.sh"         
      done
   done
done
