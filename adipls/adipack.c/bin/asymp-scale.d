#!/bin/bash
#  script to call  $aprgdir/auxprg/asymp-scale.d.x
if [ $# -eq 0 ]
then
	"$aprgdir"/auxprg/asymp-scale.d.x
else
	"$aprgdir"/bin/get-input "$1" | "$aprgdir"/auxprg/asymp-scale.d.x
fi	
