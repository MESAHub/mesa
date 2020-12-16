#!/bin/bash
#  script to change between binary and ASCII grand summary files
if [ $# -eq 0 ] || [ "$1" == "-help" ] 
then
	echo "Usage: form-agsm.d -ltcase-gt -ltinput file-gt -ltoutput file-gt [ishort]"
	echo "If ishort is set and equal to 1, print only first mode"
	echo "for each degree in diagnostic output"
	exit 1
fi

ib=${$4:-0}

echo -e "$1\n$2\n$3\n$ib" | "$aprgdir"/adiajobs/form-agsm.d.x
