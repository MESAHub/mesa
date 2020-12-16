#!/bin/bash
#  script to call  $aprgdir/adiajobs/redistrb.c.d.x
if [ $# -eq 0 ]
then
	"$aprgdir"/adiajobs/redistrb.c.d.x
elif [ "$1" == "-help" ]
 then
	echo "Usage: redistrb.d [control file]"
	exit 1
else
	"$aprgdir"/bin/get-input "$1" | "$aprgdir"/adiajobs/redistrb.c.d.x
fi
