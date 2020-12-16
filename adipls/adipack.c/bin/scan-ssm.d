#!/bin/bash
#  script to call scan-ssm.d
if [ $# -eq 0 ]
then
	echo "usage:  scan-ssm -ltinput file-gt"
	exit 1
elif [ "$1" == "-help" ]
then
	echo "usage:  scan-ssm -ltinput file-gt" 
	exit 1
fi

echo "$1" | "$aprgdir"/adiajobs/scan-ssm.d.x 
