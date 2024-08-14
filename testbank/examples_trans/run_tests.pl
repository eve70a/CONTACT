#!/usr/bin/perl -w

use strict;

# specific settings for Dos/Windows versus Unix/Linux filesystems:
my $perl_platform = $^O;
my $dir_separator;
if ($perl_platform =~ /.*win.*/i) {
   $dir_separator = '\\';
} else {
   $dir_separator = '/';
}
my $contactdir = "..".$dir_separator.".."; 
#  $contactdir = $ENV{CONTACTDIR};
my $bits = '64';
my $contact;
my $tkdiff;
if ($perl_platform =~ /.*win.*/i) {
   $contact = "$contactdir\\bin\\contact_win$bits.exe";
   $tkdiff = "winmergeu";
} else {
   $contact = "$contactdir/bin/contact";
   $tkdiff = "tkdiff";
   $ENV{LD_LIBRARY_PATH} = "$contactdir/bin/linux$bits";
}
$ENV{OMP_NUM_THREADS} = 2;

my $status = 0;
my $abort_on_error = 1;
my @list=( 'force_site_b', 'iestim_site_b',
           't1_force_cp', 't1_merge_cp', 't1_npotmax', 't1_split_cp', 't1_zws',
           't3_t1_conformal', 't3_t1_force', 't3_t2_conformal', 't3_t2_force',
           'm1_trans_dqrel', 'm1_steady_roll', 'm1_catt_shift', 'm1_trans_roll', 'm1_totforce' );
#  @list=( 't1_force_cp' );

foreach my $test ( @list ) {
   if ( !$status ) {

      # run CONTACT, get exit code
      system($contact . " 2 $test");
      $status=$?;

      if ( $status ) {
         print "An error occurred in CONTACT, skipping diff.\n";
         if ( not $abort_on_error ) {
            $status = 0;
         }
      } else {

         # run diff program, get status
         system($tkdiff . " $test.out $test.ref_out");
         $status = $?;
         if ($status == 256) { $status = 0; }
      }

      if ( $status ) {
         print "Diff aborted (status=$status), aborting loop.\n";
      }
   }
}

