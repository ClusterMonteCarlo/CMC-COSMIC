#!/bin/sh

cp cmc.h cmc.h~
sed -e 's/\(#define VERSION \)\".*\"/\1\"'"$1"'\"/' -e 's/\(#define NICK \)\".*\"/\1\"'"$2"'\"/' -e 's/\(#define DATE \)\".*\"/\1\"'"`date`"'\"/' cmc.h~ > cmc.h
