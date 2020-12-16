#!/bin/bash
#  script to call  $aprgdir/adiajobs/prt-amdl.d.x

if [ $# -eq 0 ]
then
	#csh version had $aprgdir/adiajobs/prt-amdl.d.x $* but $* is $argv but thats not defined for 0 inputs?
	"$aprgdir"/adiajobs/prt-amdl.d.x 
elif [ "$1" == "-help" ] 
then
	echo "Usage: prt-amdl.d -ltinput file-gt"
	exit 1
else
	echo "$1" | "$aprgdir"/adiajobs/prt-amdl.d.x
fi
