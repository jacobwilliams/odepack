#!/bin/bash
#!(@) extract lines starting with !! and convert to html using pandoc

for DIRNAME in M_da1 M_matrix M_main
do
   for SOURCE in $(find src/$DIRNAME -name '*.inc')
   do
      BNAME=$(basename $SOURCE .inc)
      (
      echo "# $BNAME"
      grep '^!!' $SOURCE|sed -e 's/^!!//'
      echo '```fortran'
      grep -v '^!!' $SOURCE
      echo '```'
      ) |pandoc --from=commonmark --to=html --columns=80 > out.html
   done
done

