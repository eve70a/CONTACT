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
#  $contactdir = "c:\\cmcc\\releases\\contact_v20.1";
my $bits = '64';
my $contact;
my $tkdiff;
if ($perl_platform =~ /.*win.*/i) {
   $contact = "$contactdir\\bin\\contact_win$bits.exe";
#  $tkdiff = "\"c:\\program files\\WinMerge\\winmergeu.exe\"";
   $tkdiff = "winmergeu";
} else {
   $contact = "$contactdir/bin/contact_linux$bits";
   $tkdiff = "tkdiff";
   $ENV{LD_LIBRARY_PATH} = "$contactdir/bin/linux$bits";
}
$ENV{OMP_NUM_THREADS} = 2;

my $status = 0;
my $abort_on_error = 1;
my @list=( 'carter_plast', 'catt2d_plast', 'catt3d_plast', 'plastic_one_shift', 
           'mbench_a22_left', 'tractcurve_plast' );
#  @list=( 'plastic_one_shift' );
#  @list=( 'catt3d_plast' );
#  @list=( 'carter_plast' );

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

