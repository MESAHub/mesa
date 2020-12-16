#!/bin/bash
#  script to call  $aprgdir/adiajobs/freqdif.d.x
if [ $# -eq 0 ]
	"$aprgdir"/adiajobs/freqdif.d.x
elif [ "$1" == "-help" ]
then
	echo "Usage: freqdif.d [control file]"
	exit 1
else
	"$aprgdir"/bin/get-input "$1" | "$aprgdir"/adiajobs/freqdif.d.x
fi
