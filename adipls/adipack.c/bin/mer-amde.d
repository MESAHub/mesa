#!/bin/bash
#  script to call  $aprgdir/adiajobs/mer-amde.d.x
#  with no arguments, prompts
#  otherwise calling sequence is
#   mer-amde.d -ltcase-gt -ltfile1-gt -ltfile2-gt -ltoutput file-gt [-ltdiagnostics-gt]
#
#  case = 1: full eigenfunction file
#  case = 2: reduced eigenfunction file
#  case = 3: file of kernels for spherically symmetric rotation

if [ $# -eq 0 ]
then
	"$aprgdir"/adiajobs/mer-amde.d.x
	exit 0
fi

if [ $# == 1 ] 
then
	echo "calling sequence:"
	echo "mer-amde.d -ltcase-gt -ltfile1-gt -ltfile2-gt -ltoutput file-gt [-ltdiagnostics-gt]"
	echo "case = 1: full eigenfunction file"
	echo "case = 2: reduced eigenfunction file"
	echo "case = 3: file of kernels for spherically symmetric rotation"
	exit 1
fi

diag=${$5-0}

echo -e "$1\n$2\n$3\n$4\n$diag\n,,," | "$aprgdir"/adiajobs/mer-amde.d.x

