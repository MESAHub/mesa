#!/bin/bash
#  script to call  $aprgdir/adiajobs/selamde.x

if [ $# -eq 0 ]
	"$aprgdir"/adiajobs/selamde.d.x
elif [ "$1" == "-help" ]
then
	echo "Usage: selamde.d [control file]"
	exit 1
else
	"$aprgdir"/bin/get-input "$1" | "$aprgdir"/adiajobs/selamde.d.x
fi
