#!/bin/sh
# Backup FASP files as a redistributable zip file
rm -f faspsolver.zip
zip -r faspsolver.zip README INSTALL License Makefile VERSION  \
                      base data test tutorial modules util log \
                      vs19 *.txt doc/*.pdf doc/*.in doc/FAQ
