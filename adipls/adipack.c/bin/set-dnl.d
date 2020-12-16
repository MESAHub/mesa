#!/bin/bash
#  script to call  $aprgdir/adiajobs/set-dnl.d.x
if [ $# -eq 0 ]
then
	"$aprgdir"/adiajobs/set-dnl.d.x
elif [ $# -lt 8 ]
then
	echo "Usage: set-dnl.d -ltcase-gt -ltinput file-gt -ltoutput file-gt \\"
	echo "       -ltlmax-gt -ltn step-gt -ltl step-gt -ltscale-gt -ltierr-gt"
	echo "-ltcase-gt: type of modes used for input."
	echo "   case = 1: grand summary."
	echo "   case = 2: short summary."
	echo "   case = 3: observed frequencies."
	echo "   case = 4: grand summary, use (uncorrected) eigenfrequency"
	echo "   case = 5: grand summary, use Richardson extrapolated" \
			   "eigenfrequency"
	echo "-ltn step-gt: step in order (typically 1)"
	echo "-ltl step-gt: step in degree (typically 2)"
	echo "-ltscale-gt: determines scaling of difference"
	echo "   scale = 0: No scaling"
	echo "   scale = 1: scale by 3/(2*l+3)"
	echo "   scale = 2: scale by 2/(l+2)"
	echo "-ltierr-gt: if 1, read data with errors" \
		      "(only for observed frequencies)"
	exit 1
else
	echo -e "$1\n$2\n$3\n$4\n$5\n$6\n$7\n$8" | "$aprgdir"/adiajobs/set-dnl.d.x
fi
