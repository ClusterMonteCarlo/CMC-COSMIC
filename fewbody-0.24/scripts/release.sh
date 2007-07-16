#!/bin/sh

make mrproper
svn log > ChangeLog

FB_VERSION=`grep FB_VERSION fewbody.h | sed -e 's/#define FB_VERSION \"\(.*\)\"/\1/'`

DIR=`pwd | sed -e 's/^.*\/\(.*\)$/\1/'`
cd ..
mv "$DIR" "fewbody-$FB_VERSION"
tar -cvzf "fewbody-$FB_VERSION".tar.gz --exclude .svn "fewbody-$FB_VERSION"
mv "fewbody-$FB_VERSION" "$DIR"
cd "$DIR"
