#!/bin/bash
#  script to call  $aprgdir/adiajobs/setexec.d.x
if [ $# -eq 0 ]
then
	"$aprgdir"/adiajobs/setexec.d.x
elif [ "$1" == "-help" ]
then
	echo "Usage: setexec.d [control file]"
	exit 1
else
	"$aprgdir"/bin/get-input "$1" | "$aprgdir"/adiajobs/setexec.d.x
fi
