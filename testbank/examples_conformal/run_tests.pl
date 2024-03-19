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
my @list=( 'test-rollers2', 'test-rollers3', 'mbench_conf', 'groove_li' );
#  @list=( 'test-rollers3' );
#  @list=( 'groove_li' );

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
      }

      if ( $status ) {
         print "Diff aborted (status=$status), aborting loop.\n";
      }
   }
}

