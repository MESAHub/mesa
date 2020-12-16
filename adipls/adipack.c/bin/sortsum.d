#!/bin/bash
#  script to call  $aprgdir/adiajobs/sortsum.d.x
if [ $# -eq 0 ]
	"$aprgdir"/adiajobs/sortsum.d.x
elif [ "$1" == "-help" ]
then
	echo "Usage: sortsum.d [control file]"
	exit 1
else
	$aprgdir/bin/get-input "$1" | "$aprgdir"/adiajobs/sortsum.d.x
fi
