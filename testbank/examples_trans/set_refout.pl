#!/usr/bin/perl -w

use strict;
use File::Copy qw(move);

my @list=( 'force_site_b', 'iestim_site_b',
           't1_force_cp', 't1_merge_cp', 't1_npotmax', 't1_split_cp', 't1_zws',
           't3_t1_conformal', 't3_t1_force', 't3_t2_conformal', 't3_t2_force',
           'm1_trans_dqrel', 'm1_steady_roll', 'm1_catt_shift', 'm1_trans_roll', 'm1_totforce' );

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

