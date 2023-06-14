#!/usr/bin/perl -w

use strict;

# cleanup: remove all generated files; note: this removes get_times.out as well

unlink(glob("*mat *.log *.out *.subs *.fx"));
unlink(glob("core.[0-9]*"));
