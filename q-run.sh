#!/bin/bash
base="./scadyn"
geom_fold="Geometries"
T_fold="T"

size="1m7"
declare -a shapes=("ell" "ob" "pro" "sph")
declare -a compositions=("P-aC")
#declare -a compositions=("aCH-aC" "aCH-aC-aCH" "O-aC-aCH" "OFe-aC-aCH" "P-aC-aCH" "PFe-aC-aCH" "O-aC" "OFe-aC" "P-aC" "PFe-aC")

BEGIN=1
END=15
for is in $(seq 0 3); do
   for ic in $(seq 0 0); do
      for i in $(seq $BEGIN $END); do
         mesh=" --mesh "$geom_fold"/${shapes[$is]}-${compositions[$ic]}-$i.h5"
         T=" -T $T_fold/T-$size-${shapes[$is]}-${compositions[$ic]}-$i.h5"
         
         out=" -o -${shapes[$is]}-${compositions[$ic]}-$i"
         params=" -p q-calc.in"
         
         comm="$base$mesh$T$params$out"
         repl="sed -i '9s|.*|"$comm"|' script.sh"
         
         eval "$comm"
         #eval "$repl"
         #eval "sbatch script.sh"    
      done
   done
done

