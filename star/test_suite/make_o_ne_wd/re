#!/bin/bash

function check_okay {
	if [ $? -ne 0 ]
	then
		exit 1
	fi
}

echo $1
cp photos/$1 restart_photo
check_okay
date "+DATE: %Y-%m-%d%nTIME: %H:%M:%S"
if [[ -e star.exe ]];then
	./star.exe
else
	./star
fi
date "+DATE: %Y-%m-%d%nTIME: %H:%M:%S"
