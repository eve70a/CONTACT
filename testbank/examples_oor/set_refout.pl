#!/usr/bin/perl -w

use strict;
use File::Copy qw(move);

my @list=( 'rounded_flat_d09',  'chalmers_flat_fz125' );  

foreach my $test ( @list ) {
   my $out = "$test.out";
   if (-e $out) {
      my $ref = "$test.ref_out";
      print "Moving $out to $ref\n";
      move $out, $ref;
   } else {
      print "Skipping $out (not found)\n";
   }
}
