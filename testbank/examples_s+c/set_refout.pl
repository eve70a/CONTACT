#!/usr/bin/perl -w

use strict;
use File::Copy qw(move);

my @list=( 'cross_brute',  'cross_locus',  'wing_brute',  'wing_locus',
           'cross+wing',   'cw_interrupt', 'two_patches',
           'mbench_brute', 'mbench_locus', 'mbench_intrup' );
#  @list=( 'cross_brute' );

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
