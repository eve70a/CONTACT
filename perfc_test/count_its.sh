#!/bin/sh

grep 'Tang:.*It..=' $1 | \
awk 'BEGIN { itgs=0; ncase=0; }
     { itgs=itgs+$10; ncase=ncase+1; }
     END { print "ItGS: ", itgs, ", avg: " itgs/ncase }'
