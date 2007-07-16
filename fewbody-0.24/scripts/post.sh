#!/bin/sh

FB_VERSION=`grep FB_VERSION fewbody.h | sed -e 's/#define FB_VERSION \"\(.*\)\"/\1/'`

scp ../"fewbody-$FB_VERSION".tar.gz leste.astro.northwestern.edu:public_html/code/fewbody/
