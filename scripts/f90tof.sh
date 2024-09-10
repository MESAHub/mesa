#!/bin/sh

f_files=*.f90

for file in $f_files; do
   outfile=`echo $file | sed 's,\.f90,\.f,'`
   mv $file $outfile
done
