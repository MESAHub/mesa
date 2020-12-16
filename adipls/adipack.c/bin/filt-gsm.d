#!/bin/bash
#  script to call  $aprgdir/adiajobs/filt-gsm.d.x
if [ $# -eq 0 ] 
then
	"$aprgdir"/adiajobs/filt-gsm.d.x
elif [ "$1" == "-help" ] 
then
	echo "Usage: filt-gsm.d -ltinput file-gt -ltoutput file-gt -ltcase-gt -ltparam.-gt"
	echo "Here the following cases are possible:"
	echo "case = 1: Eliminate modes with excessive difference"
	echo "          between variational and Richardson frequencies"
	echo "          -ltparam.-gt is maximum allowed difference (in microHz)"
	echo "case = 2: Eliminate modes with case number differing" \
	                "from -ltparam.-gt"
	echo "case = 3: Eliminate modes with case number equal" \
	                "to -ltparam.-gt"
	exit 1
elif [ $# -lt 4 ] 
then
	"$aprgdir"/adiajobs/filt-gsm.d.x
else
	echo -e "$1\n$2\n$3\n$4" | "$aprgdir"/adiajobs/filt-gsm.d.x
fi
