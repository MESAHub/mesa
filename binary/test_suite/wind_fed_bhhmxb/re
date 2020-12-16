#!/bin/bash

photo_directory=photos

function most_recent_photo {
    (
        cd $photo_directory
        bp=$(ls -t b_* | head -1)
        echo ${bp#b_}
    )
}

if [ $# -eq 0 ]
then
    photo=$(most_recent_photo)
else
    photo=$1
fi

if [ -z "$photo" ] || ! [ -f "$photo_directory/b_$photo" ] || ! [ -f "$photo_directory/1_$photo" ]
then
    echo "Not all specified photos exist:" $photo
    exit 1
fi

echo "restart from $photo"
echo $photo > .restart

date "+DATE: %Y-%m-%d%nTIME: %H:%M:%S"
./binary
date "+DATE: %Y-%m-%d%nTIME: %H:%M:%S"
