#!/bin/bash

function check_okay {
	if [ $? -ne 0 ]
	then
		exit 1
	fi
}

svn up
check_okay
svn info > .svninfo
check_okay
VERSION=`grep "Revision" .svninfo`
check_okay
NUM=${VERSION#* }
NUM=$(($NUM+1))
echo $NUM
echo "set data/version_number"
echo $NUM > data/version_number
echo "set test_version"
echo $NUM > test_version
echo ""
echo "svn commit -m" "'"$1"'"
echo $1 > log_message
svn commit -F log_message
check_okay
echo ""
echo "data/version_number"
more data/version_number
echo ""
echo "cd star; cpv" $NUM
echo ""
date


#echo "to fix conflicts: svn resolve --accept 'working' <path to file with conflict>"


