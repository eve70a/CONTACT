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
my $clibrary;
my $tkdiff;
if ($perl_platform =~ /.*win.*/i) {
#  $contact = "$contactdir\\bin\\contact.exe";
   $contact = "$contactdir\\bin\\contact_win$bits.exe";
   $clibrary = "$contactdir\\bin\\test_varprof_win$bits.exe";
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
my @progrms=( $clibrary, $contact );
   @progrms=( $contact );
#  @progrms=( $clibrary );
my @list=( 'rounded_flat_d09',  'chalmers_flat_fz125' );  
#  @list=( 'rounded_flat_d09' );
#  @list=( 'chalmers_flat_fz125' );

foreach my $test ( @list ) {
   foreach my $prog ( @progrms ) {
      if ( !$status ) {

         if ( $prog eq $clibrary ) {
            # run CONTACT, get exit code
            print "$clibrary $test\n";
            system($clibrary . " $test");
            $status=$?;
         } else {
            # run CONTACT, get exit code
            print "$contact 2 $test\n";
            system($contact . " 2 $test");
            $status=$?;
         }
   
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
}

