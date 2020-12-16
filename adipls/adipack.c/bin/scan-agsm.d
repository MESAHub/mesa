#!/bin/bash
#  script to call scan-agsm.d
if [ $# -eq 0 ] || [ $1 == "-help" ]
then
	echo "Usage: scan-agsm.d -ltgrand summary file-gt [ishort]"
	echo "If ishort is set and equal to 1, print only first mode"
	echo "for each degree"
	exit 1
fi

ib=${$2:-0}

echo -e"$1\n$ib" | "$aprgdir"/adiajobs/scan-agsm.d.x
