#!/bin/sh

make mrproper

DIR=`pwd | sed -e 's/^.*\/\(.*\)$/\1/'`
cd ..
tar cvzf "$DIR".tar.gz --exclude CVS --exclude \*.fit "$DIR"
cd "$DIR"
