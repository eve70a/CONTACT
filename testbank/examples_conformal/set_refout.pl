#!/usr/bin/perl -w

use strict;
use File::Copy qw(move);

my @list=( 'test-rollers2', 'test-rollers3', 'mbench_conf', 'groove_li' );

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
