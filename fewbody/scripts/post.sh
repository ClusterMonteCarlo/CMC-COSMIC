#!/bin/sh

FB_VERSION=`grep FB_VERSION fewbody.h | sed -e 's/#define FB_VERSION \"\(.*\)\"/\1/'`

scp ../"fewbody-$FB_VERSION".tar.gz x.dialup.mit.edu:www/code/fewbody/
