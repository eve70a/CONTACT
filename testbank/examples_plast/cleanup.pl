#!/usr/bin/perl -w

use strict;

# cleanup: remove all generated files

unlink(glob("*.mat *.out *.subs contact_addon.inp"));
unlink(glob("core.[0-9]* *.bak dump_gap.m"));
