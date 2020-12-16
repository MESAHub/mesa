#!/bin/bash
#  script to call  $aprgdir/adiajobs/sel-amdl.d.x
#  usage: sel-amdl -ltinput file-gt -ltoutput file-gt -ltmodel number-gt
if [ $# -eq 0 ]
then
	"$aprgdir"/adiajobs/sel-amdl.d.x
elif [ $# -lt 3 ] then
	echo "Usage: sel-amdl.d -ltinput file-gt -ltoutput file-gt -ltmodel number-gt"
	exit 1
else
	echo -e "$1\n$2\n$3" | "$aprgdir"/adiajobs/sel-amdl.d.x
fi
