#!/usr/bin/perl -w

use strict;
use File::Copy qw(move);

foreach my $test ( 'm3_tempdep', 'conv_temp', 'm1_tempdep', 'euro_tractcurv' ) {
#foreach my $test ( 'cattaneo' ) {
   my $out = "$test.out";
   if (-e $out) {
      my $ref = "$test.ref_out";
      print "Moving $out to $ref\n";
      move $out, $ref;
   } else {
      print "Skipping $out (not found)\n";
   }
}
