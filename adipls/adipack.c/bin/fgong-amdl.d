#!/bin/bash
#  script to call fgong-amdl.d
if [ $# -lt 2 ] || [ "$1" == "-help" ] 
then
	echo "Usage: fgong-amdl.d -ltInput file-gt -ltoutput file-gt"
	exit 1
fi

echo -e "$1\n$2" | "$aprgdir"/adiajobs/fgong-amdl.d.x

