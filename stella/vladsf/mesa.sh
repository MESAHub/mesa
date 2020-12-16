#!/bin/bash

echo "terminal output in mesa_log.*"
echo "opacity table 1"
./xronfict.exe 1 > mesa_log.1 &
echo "opacity table 2"
./xronfict.exe 2  > mesa_log.2 &
echo "opacity table 3"
./xronfict.exe 3  > mesa_log.3 &
echo "opacity table 4"
./xronfict.exe 4  > mesa_log.4 &
echo "opacity table 5"
./xronfict.exe 5  > mesa_log.5 &
echo "opacity table 6"
./xronfict.exe 6  > mesa_log.6 &
date
wait
date
echo "all processes complete"


