#!/bin/bash
#  script to call  $aprgdir/adiajobs/res-amde.d.x
#   usage:  res-amde -ltinput file-gt -ltoutput file-gt -ltcase-gt
if [ $# -eq 0 ]
then
	echo "usage:  res-amde -ltinput file-gt -ltoutput file-gt -ltcase-gt"
	echo "Input: full eigenfunction file"
	echo "case = 1: output y(1) and y(2)"
	echo "case = 2: output zhat(1) and zhat(2)"
	exit 1
elif [ "$1" == "-help" ] 
then
	echo "usage:  res-amde -ltinput file-gt -ltoutput file-gt -ltcase-gt"
	echo "Input: full eigenfunction file"
	echo "case = 1: output y(1) and y(2)"
	echo "case = 2: output zhat(1) and zhat(2)"
	exit 1
fi
echo -e "$1\n$2\n$3" | "$aprgdir"/adiajobs/res-amde.d.x

