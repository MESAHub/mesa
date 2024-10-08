#!/bin/sh

f_files=*.f

for file in $f_files; do
   outfile=`echo $file | sed 's,\.f,\.f90,'`
   mv $file $outfile
done
