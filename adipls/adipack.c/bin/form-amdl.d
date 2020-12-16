#!/bin/bash
#  script to change between binary and ASCII pulsation model files
if [ $# -eq 0 ] || [ "$1" == "-help" ] 
then
	echo "Usage: form-amdl.d -ltcase-gt -ltinput file-gt -ltoutput file-gt"
	echo "case = 1: transform from binary to ASCII"
	echo "case = 2: transform from ASCII to binary"
	exit 1
else
	echo -e "$1\n$2 $3" | "$aprgdir"/adiajobs/form-amdl.d.x
fi
