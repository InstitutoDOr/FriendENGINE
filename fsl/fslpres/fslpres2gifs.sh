#!/bin/sh

opera ${1}/index.html -geometry 1024x881 &

echo "before wait"
sleep 30
echo "after wait"

id=`xwininfo -root -tree | grep -m 1 1024x768 | awk '{print $1}'`

echo $id;

xwd -silent -nobdrs -id $id | convert - index.pnm
pnmscale -xsize=320 -ysize=240 index.pnm | convert - index.gif 

for i in ${1}/s_*htm ; do 
  j=`basename $i .htm | sed 's/\// /g' | sed 's/s_//'`

  opera ${i} &
  sleep 8
  xwd -silent -nobdrs -id $id | convert - ${j}.pnm
  pnmscale -xsize=320 -ysize=240 ${j}.pnm | convert - ${j}.gif
  echo $j.gif
done