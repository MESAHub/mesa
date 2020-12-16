#!/bin/bash
#  script to call  $aprgdir/adiajobs/set-asscal.d.x

if [ $# -lt 2 ]
then
  echo "Usage: set-asscal.d -ltinput file-gt -ltoutput file-gt [-lttruncation-gt]"
  echo "       Here -lttruncation-gt defines the truncation radius"
  echo "       Default is 1 (i.e., the photosphere)"
  exit 1
fi

xtrscl=${$3:-1}

#Intentional space left between -1 and newline
echo -e "2 $1\n3 $2\n-1  \n$xtrscl" | "$aprgdir"/adiajobs/set-asscal.d.x

