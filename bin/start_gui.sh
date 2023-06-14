#!/bin/sh

if [ -z "$CONTACTDIR" ]; then
   echo "ERROR: environment variable CONTACTDIR is not defined."
   echo "       Please refer to README.txt."
   exit
fi

java -jar $CONTACTDIR/bin/contact_gui.jar &
