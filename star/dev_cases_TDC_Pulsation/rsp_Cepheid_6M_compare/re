#!/bin/bash

shopt -u expand_aliases 

photo_directory=photos

function most_recent_photo {
    ls -tp "$photo_directory" | grep -v / | head -1
}

if [ $# -eq 0 ]
then
    photo=$(most_recent_photo)
else
    photo=$1
fi

if [ -z "$photo" ] || ! [ -f "$photo_directory/$photo" ]
then
    echo "specified photo ($photo) does not exist"
    exit 1
fi

echo "restart from $photo"
if ! cp "$photo_directory/$photo" restart_photo
then
    echo "failed to copy photo ($photo)"
    exit 1
fi

date "+DATE: %Y-%m-%d%nTIME: %H:%M:%S"
./star
date "+DATE: %Y-%m-%d%nTIME: %H:%M:%S"
