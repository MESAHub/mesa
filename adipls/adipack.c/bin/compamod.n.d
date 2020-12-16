#!/bin/bash
#  script to call  $aprgdir/adiajobs/compamod.n.d.x
if [ $# -eq 0 ] 
then
	"$aprgdir"/adiajobs/compamod.n.d.x
elif [ "$1" == "-help" ] 
then
	echo "Usage: compamod.n.d [control file]"
	exit 1
else
	"$aprgdir"/bin/get-input "$1" | "$aprgdir"/adiajobs/compamod.n.d.x
fi
