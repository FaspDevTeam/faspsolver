#!/bin/sh

# backup server info
SERVER="45.79.93.25"
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
rsync -avz -e 'ssh -p 8522' $1 ${OPT} ./htdocs/* ${USERID}@${SERVER}:/www/wwwroot/11/fasp
