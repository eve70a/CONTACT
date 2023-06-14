#!/usr/bin/perl -w

use strict;
use File::Copy qw(move);

my @list=( 'stand_alone', 'contact_addon', 'test_mbench', 'test_table' );
#  @list=( 'test_table' );
my @files;

foreach my $test ( @list ) {

   # list reference outputs for each test:

   if ($test =~ /stand_alone/i) {
      @files = ( 'testbank.ref_out', 'dry_10mph_13750lbs.ref_out' );
   } elsif ($test =~ /contact_addon/i) {
      @files = ( 'test_caddon.ref_out', 'contact_addon.ref_out', 'contact_addon.ref_inp' );
   } elsif ($test =~ /test_mbench/i) {
      @files = ( 'test_mbench.ref_out', 'contact_mbench.ref_out' );
   } elsif ($test =~ /test_table/i) {
      @files = ( 'test_table_mx11.ref_out', 'test_table_mx11.ref_fx' );
   }

   # loop over all reference outputs per test:

   foreach my $ref_file ( @files ) {
      my $out_file = $ref_file;
      $out_file =~ s/\.ref_/./;
      if (-e $out_file) {
         print "Moving $out_file to $ref_file\n";
         move $out_file, $ref_file;
      } else {
         print "Skipping $out_file (not found)\n";
      }
   }
}

