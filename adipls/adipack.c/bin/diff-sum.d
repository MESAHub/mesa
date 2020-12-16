#!/bin/bash
#  script to call  $aprgdir/adiajobs/diff-sum.d.x
if [ $# -eq 0 ] 
then
	"$aprgdir"/adiajobs/diff-sum.d.x
elif [ $1 == "-help" ] || [ $# -lt 6 ] 
then
	echo "Usage: diff-sum.d -ltcase 1-gt -ltcase 2-gt -ltinput 1-gt -ltinput 2-gt" \
			        "-ltexcess 1-gt -ltexcess 2-gt"
	echo "-ltcase-gt determines type of input file"
	echo "   case = 1: Grand summary"
	echo "   case = 2: Short summary"
	echo "   case = 3: ASCII file of "observed" frequencies"
	echo "-ltexcess 1-gt: Output file of modes in file 1 but not in file 2"
	echo "-ltexcess 2-gt: Output file of modes in file 2 but not in file 1"
	echo "Note: -ltexcess-gt files are output in same form as" \
	      "corresponding input file"
	exit 1
else
	echo -e "$1\n$2\n$3\n$4\n$5\n$6" | "$aprgdir"/adiajobs/diff-sum.d.x
fi
