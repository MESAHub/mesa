#!/bin/bash
#  script to call  $aprgdir/adiajobs/lsqfreq.d.x
if [ $# -eq 0 ]
then
	"$aprgdir"/adiajobs/lsqfreq.d.x
elif [ "$1" == "-help" ] then
	echo "Usage: lsqfreq.d [control file]"
	exit 1
else
	"$aprgdir"/bin/get-input "$1" | "$aprgdir"/adiajobs/lsqfreq.d.x
fi
