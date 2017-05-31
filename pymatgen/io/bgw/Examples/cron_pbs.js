#!/bin/bash

#PBS -l select=1:ncpus=1
#PBS -l walltime=12:00:00
#PBS -q transfer
#PBS -A WPASC96170001
#PBS -N transfer
#PBS -o transfer.out
#PBS -e transfer.out
#PBS -j oe

file="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
echo file location: $file

# This can be changed by the User for customization
FW_DIR="FIREWORKS_JOBS"

cd $WORKDIR
if [ ! -d $FW_DIR ]; then
    mkdir $FW_DIR
fi

cd $FW_DIR

for i in `seq 120`; do
    echo "Checking for Fireworks and launching"
    date
    qlaunch -r -c ~/.fireworks_qe rapidfire -b 10 -m 5
    sleep 300
done

write_cron ()
{
   cat $file > ~/cron_pbs.js 
}

if [ ! -f ~/cron_pbs.js ]; then
    write_cron
fi  

qsub ~/cron_pbs.js

