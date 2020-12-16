#!/bin/bash
#  script to call  $aprgdir/adiajobs/scan-amdl.d.x
if [ $# -eq 0 ] || [ "$1" == "-help" ]
then
	echo "Usage: scan-amdl.d -ltinput file-gt"
	exit 1
else
	echo "$1" | "$aprgdir"/adiajobs/scan-amdl.d.x
fi
