#!/bin/bash
#  script to call  $aprgdir/adiajobs/selsum.d.x
if [ $# -eq 0 ]
then
	"$aprgdir"/adiajobs/selsum.d.x
elif [ "$1" == "-help" ]
then
	echo "Usage: selsum.d [control file]"
	exit 1
else
	"$aprgdir"/bin/get-input "$1" | "$aprgdir"/adiajobs/selsum.d.x
fi
