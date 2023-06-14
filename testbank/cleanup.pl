#!/usr/bin/perl -w

use strict;

# cleanup: remove all generated files

unlink(glob("core.[0-9]* *.bak *.fx *.mat *.out *.subs"));
unlink(glob("dump_gap*.txt dump_gap*.m"));
unlink("spck_addon.inp");
unlink("um_addon.inp");
unlink("nucars_addon.inp");
