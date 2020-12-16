#!/bin/bash
#  script to call  $aprgdir/adiajobs/mer-sum.d.x
#  with no arguments, prompts
#  otherwise calling sequence is
#   mersum -ltcase-gt -ltfile1-gt -ltfile2-gt -ltoutput file-gt [-ltdiagnostics-gt]

if [ $# -eq 0 ]
then
	"$aprgdir"/adiajobs/mer-sum.d.x
	exit 0
fi

if [ $# -lt 4 ] 
then
	echo "Usage: mer-sum.d -ltcase-gt -ltfile1-gt -ltfile2-gt -ltoutput file-gt [-ltdiag.-gt]"
	echo "-ltcase-gt: type of modes used for input "
	echo "   case = 1: grand summary."
	echo "   case = 2: short summary."
	echo "   case = 3: observed frequencies, without errors."
	echo "   case = 4: observed frequencies, with errors."
	exit 1
fi

diag=${$5-0}

echo -e "$1\n$2\n$3\n$4\n$diag\n,,," | "$aprgdir"/adiajobs/mer-sum.d.x

