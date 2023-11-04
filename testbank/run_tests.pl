#!/usr/bin/perl -w

use strict;
use File::Copy qw(move);

# specific settings for Dos/Windows versus Unix/Linux filesystems:
my $perl_platform = $^O;
my $dir_separator;
my $bits = '64';
my $contact_exe;
my $test_nonhz;
my $test_mbench;
my $test_table;
my $nc_addon;
my $tkdiff;
if ($perl_platform =~ /.*win.*/i) {
   $dir_separator = '\\';
   $contact_exe = "..\\bin\\contact_win$bits.exe";
   $test_nonhz = "..\\bin\\test_caddon_win$bits.exe";
   $test_mbench = "..\\bin\\test_mbench_win$bits.exe";
   $test_table = "..\\bin\\test_table_win$bits.exe"; 
#  $tkdiff = "\"c:\\program files\\WinMerge\\winmergeu.exe\"";
   $tkdiff = "winmergeu";
} else {
   $dir_separator = '/';
   $contact_exe = "../bin/contact_linux$bits";
   $test_nonhz = "../bin/test_caddon_linux$bits";
   $test_mbench = "../bin/test_mbench_linux$bits";
   $test_table = "../bin/test_table_linux$bits"; 
   $tkdiff = "tkdiff";
   $ENV{LD_LIBRARY_PATH} = "../bin:../bin/linux$bits";
}

my $status = 0;
my @list=( 'stand_alone', 'test_nonhz', 'test_mbench', 'test_table' );
#  @list=( 'test_nonhz', 'test_mbench', 'test_table' );
#  @list=( 'stand_alone' );
#  @list=( 'test_nonhz' );
#  @list=( 'test_mbench' );
#  @list=( 'test_table' );
my @files;

foreach my $test ( @list ) {
   if ( !$status ) {

      # run CONTACT, get exit code

      if ( $test =~ /stand_alone/i and -e "$contact_exe" ) {

         # Testing the stand-alone executable

         $ENV{OMP_NUM_THREADS} = 1;
         system($contact_exe . " 2 testbank");
         $status=$?;
         if ( !$status ) {
            system($contact_exe . " 2 dry_10mph_13750lbs");
            $status=$?;
         }
         @files = ( 'testbank.ref_out', 'dry_10mph_13750lbs.ref_out' );

      } elsif ( $test =~ /test_nonhz/i and -e "$test_nonhz" ) {

         # Testing the CONTACT library version for basic non-Hertzian contact, module 3

         $ENV{OMP_NUM_THREADS} = 1;
         system($test_nonhz . " > test_caddon.out");
         $status=$?;
         move "contact_addon.out", "contact_nonhz.out";
         move "contact_addon.inp", "contact_nonhz.inp";
         @files = ( 'test_caddon.ref_out', 'contact_nonhz.ref_out', 'contact_nonhz.ref_inp' );

      } elsif ( $test =~ /test_mbench/i and -e "$test_mbench" ) {

         # Testing the CONTACT library version for wheel/rail contact, module 1

         $ENV{OMP_NUM_THREADS} = 1;
         system($test_mbench . " > test_mbench.out");
         $status=$?;
         move "contact_addon.out", "contact_mbench.out";
         @files = ( 'test_mbench.ref_out', 'contact_mbench.ref_out' );

      } elsif ( $test =~ /test_table/i and -e "$test_table" ) {

         # Testing the Table program (using CONTACT library)

         $ENV{OMP_NUM_THREADS} = 4;
         system($test_table . " 11 test_table_mx11 > test_table_mx11.fx");
         $status=$?;
         @files = ( 'test_table_mx11.ref_out', 'test_table_mx11.ref_fx' );

      } else {

         print "\n ** Invalid test ($test) or executable not found.\n";
         @files = ();
      }

      if ( $status ) {
         print "An error occurred in CONTACT, skipping diff.\n";
         $status = 0;

      # compare output-files with reference results, get exit code

      } else {

         foreach my $ref_file ( @files ) {
            if ( !$status ) {

               my $out_file = $ref_file;
               $out_file =~ s/\.ref_/./;
               system($tkdiff . " $out_file $ref_file");
               $status = $?;

               # normal exit of tkdiff is ok:
               if ($status == 256) { $status = 0; }
            }
         } 
      }

      if ( $status ) {
         print "Diff aborted (status=$status), aborting loop.\n";
      }
   }
}

