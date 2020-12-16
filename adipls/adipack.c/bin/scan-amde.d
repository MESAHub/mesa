#!/bin/bash
#  script to call  $aprgdir/adiajobs/scan-amde.d.x
if [ $# -eq 0 ] || [ "$1" == "-help" ]
then
	echo "Usage: scan-amde.d -ltfile-gt [case]"
	echo "       case = 1: full set of eigenfunctions"
	echo "       case = 2: restricted set of eigenfunctions"
	echo "       default case: 2"
	exit 1
fi 

case=${$2:-2}

echo -e "$1\n$case" | "$aprgdir"/adiajobs/scan-amde.d.x
