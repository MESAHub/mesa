#!/bin/bash
#  script to call  $aprgdir/adiajobs/set-obs.d.x
#  usage (with arguments) set-obs.d -ltcase-gt -ltinput file-gt [-ltoutput file-gt]
#  If output file is not given, assumes that -ltinput file-gt contains 
#  trailer to be applied to agsm on input and obs on output.
#  with no arguments: prompts

if [ $# -eq 1 ] && [ "$1" == "-help" ]
then
  echo "usage (with arguments) set-obs.d -ltcase-gt -ltinput file-gt [-ltoutput file-gt]"
  echo "If output file is not given, assumes that -ltinput file-gt contains "
  echo "trailer to be applied to agsm on input and obs on output."
  echo "with no arguments: prompts"
  echo "case:"
  echo "1: grand summary, variational frequency."
  echo "2: short summary."
  echo "4: grand summary, from eigenfrequency in cs(20)."
  echo "   Note that this allows setting Cowling approximation frequency"
  echo "5: grand summary, from Richardson extrapolation frequency"
  echo "6: grand summary, from (possibly corrected) eigenfrequency in cs(21)"
  echo "If icasein gt 10, set according to icasein-10, including mode energy"
  exit 1
elif [ $# -lt 2 ] 
then
	"$aprgdir"/adiajobs/set-obs.d.x
elif [ $# -eq 2 ] 
then
	echo -e "$1""\n""agsm.$2""\n""obs.$2""\n""2" | "$aprgdir"/adiajobs/set-obs.d.x
else
	echo -e "$1\n$2\n$3\n2" | "$aprgdir"/adiajobs/set-obs.d.x
fi
