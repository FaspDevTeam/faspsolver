#!/bin/sh

# backup server info
#SERVER="66.228.35.193"
SERVER="23.239.26.149"
#USERID="vcc"
USERID="root"

# optional parameters
OPT="--delete --progress --times --exclude *~"

# upload: Should be ran in doc directory with htdocs as a subdir
if [ ! -d htdocs/download ]; then
    mkdir htdocs/download
fi

# cp ../../download/.ht* htdocs/download
cp ../../download/*.zip htdocs/download
cp ../../download/*.pdf htdocs/download
rsync -avz $1 ${OPT} ./htdocs/* ${USERID}@${SERVER}:/opt/www/fasp
