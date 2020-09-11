#!/bin/sh

cp fewbody.h fewbody.h~
sed -e 's/\(#define FB_VERSION \)\".*\"/\1\"'"$1"'\"/' -e 's/\(#define FB_NICK \)\".*\"/\1\"'"$2"'\"/' -e 's/\(#define FB_DATE \)\".*\"/\1\"'"`date`"'\"/' fewbody.h~ > fewbody.h
