#!/usr/bin/perl -w

use strict;

# specific settings for Dos/Windows versus Unix/Linux filesystems:
my $perl_platform = $^O;
my $dir_separator;
my $bits = '64';
my $contact;
my $tkdiff;
if ($perl_platform =~ /.*win.*/i) {
   $dir_separator = '\\';
   $contact = "..\\bin\\contact_win$bits.exe";
#  $tkdiff = "\"c:\\program files\\WinMerge\\winmergeu.exe\"";
   $tkdiff = "winmergeu";
} else {
   $dir_separator = '/';
   $contact = "../bin/contact_linux$bits";
   $tkdiff = "tkdiff";
   $ENV{LD_LIBRARY_PATH} = "../bin/linux$bits";
}
$ENV{OMP_NUM_THREADS} = 2;

my $status = 0;
my @list=( 'cattaneo', 'carter2d', 'bentall', 'visc_cylindr', 'subsurf', 
           'mbench_a22_left', 'mbench_a22_right', 'mbench_a22_varfric', 'mbench_a22_kpec', 
           'fastsim', 'spence35', 'catt_to_cart', 'conformal', 'veldep_fric', 
           'ertz_temperature', 'plastic_3bl', 'tractcurv', 'wheelflat' );
#  @list=( 'mbench_a22_left', 'mbench_a22_right', 'mbench_a22_varfric', 'mbench_a22_kpec' );
#  @list=( 'mbench_a22_varfric' );
#  @list=( 'fastsim' );
#  @list=( 'wheelflat' );
#  @list=( 'veldep_fric' );
#  @list=( 'subsurf', 'spence35' );
#  @list=( 'ertz_temperature' );

foreach my $test ( @list ) {
   if ( !$status ) {

      # run CONTACT, get exit code
      system($contact . " 2 $test");
      $status=$?;

      if ( $status ) {
         print "An error occurred in CONTACT, skipping diff.\n";
         $status = 0;
      } else {

         # run diff program, get status
         system($tkdiff . " $test.out $test.ref_out");
         $status = $?;
         if ($status == 256) { $status = 0; }

         if ( !$status && $test =~ /subsurf/ ) {
            # run diff on subs-file, get status
            print " Compare subs...\n";
            system($tkdiff . " $test.0002.subs $test.ref_subs");
            $status = $?;
            if ($status == 256) { $status = 0; }
         }
      }

      if ( $status ) {
         print "Diff aborted (status=$status), aborting loop.\n";
      }
   }
}

