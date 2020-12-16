#!/bin/bash
#  script to call  $aprgdir/adiajobs/fdif.sin.d.x
if [ $# -eq 0 ] 
then
	"$aprgdir"/adiajobs/fdif.sin.d.x
elif [ "$1" == "-help" ] 
then
	echo "Usage: fdif.sin.d [control file]"
	exit 1
else
	"$aprgdir"/bin/get-input "$1" | "$aprgdir"/adiajobs/fdif.sin.d.x
fi
