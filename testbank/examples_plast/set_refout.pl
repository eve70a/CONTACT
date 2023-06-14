#!/usr/bin/perl -w

use strict;
use File::Copy qw(move);

foreach my $test ( 'carter_plast', 'catt2d_plast', 'catt3d_plast', 'tractcurve_plast',
                   'plastic_one_shift', 'mbench_a22_left' ) {
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
