#!/usr/bin/perl -w

use strict;

my $mode = $ARGV[0];   # choose "exe" or "clib"
if (not defined $mode) {
   $mode = "both";
}

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
   $contact = "$contactdir\\bin\\contact.exe";
#  $contact = "$contactdir\\bin\\contact_win$bits.exe";
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
my @progrms=( $clibrary, $contact );
#  @progrms=( $contact );
#  @progrms=( $clibrary );
my @list=( 'cross_brute',  'cross_locus',  'wing_brute',  'wing_locus',
           'cross+wing',   'cw_interrupt', 'two_patches',
           'mbench_brute', 'mbench_locus', 'mbench_intrup' );
#  @list=( 'cross_brute',  'cross_locus',  'wing_brute',  'wing_locus',
#          'cross+wing',   'cw_interrupt', 'two_patches' );
#  @list=( 'mbench_brute', 'mbench_intrup' );
#  @list=( 'mbench_intrup' );
   @list=( 'two_patches' );

foreach my $test ( @list ) {
   foreach my $prog ( @progrms ) {
      if ( ($prog eq $clibrary) and ($mode eq "exe") ) {

         print "Mode = $mode, skipping program $clibrary\n";

      } elsif ( ($prog eq $contact) and ($mode eq "clib") ) {

         print "Mode = $mode, skipping program $contact\n";

      } elsif ( !$status ) {

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
}

