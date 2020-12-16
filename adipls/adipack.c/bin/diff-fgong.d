#!/bin/bash
#  script to call  $aprgdir/adiajobs/diff-fgong.d.x
if [ $# -eq 0 ] 
then
	"$aprgdir"/adiajobs/diff-fgong.d.x
	exit 0
elif [ $# -lt 3 ] 
then
	echo 'Usage: diff-fgong.d -ltfirst file-gt -ltsecond file-gt -ltoutput file-gt \'
	echo "                    [case]" 
	echo "Here case is an optional case number (default 1):"
	echo "1: differences at fixed r/R "
	echo "2: differences at fixed q"
	echo "3: differences at fixed r/r(last point)"
	echo "4: differences at fixed p"
	echo "5: differences at fixed (R-r)/R"
	exit 1
fi

in1=$1
in2=$2
out=$3
case=${$4:-1}

echo -e "$in1\n$in2\n1 1\n$out\n$case\n0\n0" | "$aprgdir"/adiajobs/diff-fgong.d.x

