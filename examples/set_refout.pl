#!/usr/bin/perl -w

use strict;
use File::Copy qw(move);

my @list=( 'cattaneo', 'carter2d', 'bentall', 'visc_cylindr', 'tractcurv', 'fastsim', 'subsurf', 
           'mbench_a22_left', 'mbench_a22_right', 'mbench_a22_varfric', 'mbench_a22_kpec', 
           'spence35', 'catt_to_cart', 'conformal', 'veldep_fric', 'ertz_temperature', 'plastic_3bl' );
foreach my $test ( @list ) {
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

foreach my $test ( 'subsurf' ) {
   my $subs = "$test.0002.subs";
   if (-e $subs) {
      my $ref = "$test.ref_subs";
      print "Moving $subs to $ref\n";
      move $subs, $ref;
   } else {
      print "Skipping $subs (not found)\n";
   }
}
