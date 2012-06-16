#!/bin/sh

# backup server info
SERVER="web.sourceforge.net"
USERID="faspsolver"

# optional parameters
OPT="--delete --progress --times --exclude *~"

# upload: Should be ran in doc directory with html as a subdir
rsync -avz $1 ${OPT} ./htdocs ${USERID}@${SERVER}:/home/project-web/fasp
