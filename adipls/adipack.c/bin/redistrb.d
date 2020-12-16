#!/bin/bash
#  script to call  $aprgdir/adiajobs/redistrb.d.x
if [ $# -eq 0 ]
then
	"$aprgdir"/adiajobs/redistrb.d.x
elif [ "$1" == "-help" ]
 then
	echo "Usage: redistrb.d [control file]"
	exit 1
else
	"$aprgdir"/bin/get-input "$1" | "$aprgdir"/adiajobs/redistrb.d.x
fi
