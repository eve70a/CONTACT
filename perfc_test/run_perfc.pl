#!/usr/bin/perl -w
use strict;
use Sys::Hostname;
use POSIX;

my $HOST = hostname;
my $USER = getlogin();
#
# specific settings for Dos/Windows versus Unix/Linux filesystems:
my $perl_platform = $^O;
my $file_sep;
if ($perl_platform =~ /.*win.*/i) {
   $file_sep = '\\';
} else {
   $file_sep = '/';
}

# redirect stdout to log-file

open(STDOUT, ">run_perfc.log");
open(STDERR, ">run_perfc.log");

# prescribe which CONTACT version to use

my $VER="contact-open";
#  $VER="releases".$file_sep."contact_v23.2";

print " ----------------------------------------------------------------------\n";
my $time_str = strftime("%d-%m-%Y %H:%M:%S", localtime);
print " Starting performance test at $time_str\n";
print " Running at host $HOST ($perl_platform), started by $USER\n";
print " Using CONTACT version=$VER\n";
print " ----------------------------------------------------------------------\n";

# configuration: run all tests, or skip the finest grid tests?

my $test_norm=1;
my $test_tang=1;
my $test_mbench=1;
my $test_switch=1;
my $test_subsurf=1;
my $test_clib=1;
my $test_2dcases=1;
my $test_parall=1;

# shortest tests only: about  5 min
# including "1min":    about 20 min
# including "10min":   about  1 hour

# run tests that may take around one minute / 10 minutes to compute?
my $test_1min=1;
my $test_10min=1;

my $BASE;
if ($perl_platform =~ /.*win.*/i) {
   $BASE="c:\\CMCC";
} elsif ($HOST =~ /vortech/i || $HOST =~ /mail.cmcc/i ) {
   $BASE="/v3/CMCC";
   $ENV{'LD_LIBRARY_PATH'} = "$BASE/$VER/bin:$BASE/$VER/bin/linux64";
} else {
   $BASE="/home/wagm/edwin/repos";
}

my ($cntc, $clib);
if ($perl_platform =~ /.*win.*/i && $VER =~ /release/ ) {
   $cntc = "$BASE\\$VER\\bin\\contact.exe";
   $clib = "$BASE\\$VER\\bin\\test_table.exe";
} elsif ($perl_platform =~ /.*win.*/i ) {
   $cntc = "$BASE\\$VER\\bin\\contact_win64.exe";
   $clib = "$BASE\\$VER\\bin\\test_table_win64.exe";
} elsif ($VER =~ /release/ ) {
   $cntc = "$BASE/$VER/bin/contact";
   $clib = "$BASE/$VER/bin/test_table_linux64";
   # $clib = "$BASE/$VER/bin/test_table";
} else {
   $cntc = "$BASE/$VER/bin/contact_linux64";
   $clib = "$BASE/$VER/bin/test_table_linux64";
}


my $pfx="";
# my $pfx="taskset -c 1";

# run tests sequentially unless specified otherwise

$ENV{OMP_NUM_THREADS}=1;

# run tests for the normal contact problem

if ( $test_norm == 1 ) {
   my @list=("1p", "2p", "1f", "2f", "4p", "8p", "4f" );
   foreach my $tst ( @list ) {
      print " Starting problem norm_problm_$tst -------------------------\n";
      system("$pfx $cntc 2 norm_problm_$tst");
   }
   print " Starting problem norm_potcon_1p -------------------------------\n";
   system("$pfx $cntc 2 norm_potcon_1p");
}


# run tests for the tangential contact problem

if ( $test_tang == 1 ) {
   my @list=("1c", "1f", "1s", "2c", "2s", "4s");
   if ( $test_1min > 0 ) {
      @list=(@list, "8s", "4c", "2f");
   }
   if ( $test_10min > 0 ) {
      @list=(@list, "8c", "4f");
   }
   foreach my $tst ( @list ) {
      print " Starting problem tang_problm_$tst -------------------------\n";
      system("$pfx $cntc 2 tang_problm_$tst");
   }
   print " Starting problem tang_potcon_1c -------------------------------\n";
   system("$pfx $cntc 2 tang_potcon_1c");
   print " Starting problem tang_cnvxgs_1c -------------------------------\n";
   system("$pfx $cntc 2 tang_cnvxgs_1c");
}

# run tests for W/R contact module 1

if ( $test_mbench == 1 ) {
   my @list=("20kn_coarseprf", "20kn_dx05_zpos", "20kn_dx02_zpos",
             "100kn_coarseprf", "100kn_dx05_zpos" );
   if ( $test_1min > 0 ) {
      @list=(@list, "20kn_dx02_fz", "100kn_dx02_zpos", "100kn_dx02_fz" );
   }
   foreach my $tst ( @list ) {
      print " Starting problem mbench_$tst -------------------------\n";
      system("$pfx $cntc 2 mbench_$tst");
   }
}

# run tests for switches and crossings

if ( $test_switch == 1 ) {
   my @list=("mbench_intrup" );
   if ( $test_1min > 0 ) {
      @list=(@list, "mbench_locus", "mbench_brute" );
   }
   foreach my $tst ( @list ) {
      print " Starting problem $tst -------------------------\n";
      system("$pfx $cntc 2 $tst");
   }
}

# run 2D test-cases

if ( $test_2dcases == 1 ) {
   my @list=("carter2d_100", "carter2d_800", "visc_cylindr", "veldep_fric");
   foreach my $tst ( @list ) {
      print " Starting problem $tst --------------------------\n";
      system("$cntc 2 $tst");
   }
}

# run tests with CONTACT library for small grids, creating a table

if ( $test_clib == 1 ) {
   my @list=("11", "15", "19");
   if ( $test_1min > 0 ) {
      @list=(@list, "25");
   }
   if ( $test_10min > 0 ) {
      @list=(@list, "35");
   }
   foreach my $tst ( @list ) {
      print " Starting problem test_clib_seq$tst --------------------------\n";
      system("$pfx $clib $tst test_clib_seq$tst > test_clib_seq$tst.fx");
   }
}

# run parallel tests

delete $ENV{OMP_NUM_THREADS};

if ( $test_parall == 1 && $test_tang == 1 ) {

   my @list=("1c", "2c", "1s", "2s", "4s");
   if ( $test_1min > 0 ) {
      @list=(@list, "4c");
   }
   if ( $test_10min > 0 ) {
      @list=(@list, "8c");
   }
   foreach my $tst ( @list ) {
      print " Starting problem tang_parall_$tst"." --------------------------\n";
      system("$cntc 2 tang_parall_$tst");
   }
}

# run parallel tests with CONTACT library for small grids, creating a table

if ( $test_parall == 1 && $test_clib == 1 ) {
   my @list=("11", "15", "19");
   if ( $test_1min > 0 ) {
      @list=(@list, "25");
   }
   if ( $test_10min > 0 ) {
      @list=(@list, "35");
   }
   foreach my $tst ( @list ) {
      print " Starting problem test_clib_par$tst --------------------------\n";
      system("$pfx $clib $tst test_clib_par$tst > test_clib_par$tst.fx");
   }
}

# run tests for Panagiotopoulos and subsurface stresses (parallel runs)

if ( $test_subsurf == 1 ) {
   my @list=("spence35_nosubs", "spence35_135pt", "spence35_2025pt", "spence71_nosubs");
   if ( $test_1min > 0 ) {
      @list=(@list, "spence71_273pt", "spence71_8281pt");
   }
   foreach my $tst ( @list ) {
      print " Starting problem $tst --------------------------\n";
      system("$cntc 2 $tst");
   }
}

# cleanup: remove all .bac and .p files

unlink(glob("*.bac *.p *.subs"));

print " ----------------------------------------------------------------------\n";
$time_str = strftime("%d-%m-%Y %H:%M:%S", localtime);
print " Performance test complete at $time_str\n";
print " ----------------------------------------------------------------------\n";

# present overview of performance quantities

# ./get_times.sh > get_times.out

