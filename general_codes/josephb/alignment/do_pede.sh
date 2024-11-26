#! /bin/bash

steer="steer.txt"
rm -rf $steer
ls -1 $(pwd)/dat/*.bin >> $steer 2> /dev/null
echo "scaleerrors 5 5" >> $steer 2> /dev/null

pede $steer | tee millepede.out

echo "done"

