#!/bin/bash
#  script to call  $aprgdir/adipls/adipls.c.d.x
if [ $# -eq 0 ] 
then
	"$aprgdir"/adipls/adipls.c.d.x
elif [ "$1" == "-help" ]
then
	echo "Usage: adipls.c.d [control file]"
	exit 1
else
	"$aprgdir"/bin/get-input "$1" | "$aprgdir"/adipls/adipls.c.d.x
fi 


