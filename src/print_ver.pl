#!/usr/bin/perl
use strict;

my $strip = 0;  # strip (1) or add (0) version information?

# specific settings for Dos/Windows versus Unix/Linux filesystems:
my $perl_platform = $^O;

#-----------------------------------------------------------------------------
# check usage: print_ver.pl [ -strip ]
#-----------------------------------------------------------------------------
my $nargs = $#ARGV + 1;
if ( $nargs < 0 or $nargs > 1 ) { die "Usage: print_ver.pl [ -strip ]\n"; }
if ( $nargs == 0 ) {
   $strip = 0;
} else {
   if ( $ARGV[0] =~ m/-strip/ ) {
      $strip = 1;
   } else {
      # if ( $ARGV[0] =~ m/-h\>/ or $ARGV[0] =~ m/--help\>/ or $ARGV[0] =~ m/-\?\>/ ) {
      die "Usage: print_ver.pl [ -strip ]\n";
   }
}

#-----------------------------------------------------------------------------
# basic mode: add one line to the VERSION file, containing the svn version id
#-----------------------------------------------------------------------------

if ( $strip == 0 ) {
   my $ver;
   if ($perl_platform =~ /.*win.*/i) {
      system('subwcrev.exe . wcrev.txt wcrev.out -q');
      $ver = `type wcrev.out`; chomp($ver);
   } else {
      $ver = `svnversion ..`; chomp($ver);
   }

   open  VER, ">>VERSION";
   print VER "      version(7) = '(Subversion revision $ver)'\n";
   close VER;

#-----------------------------------------------------------------------------
# strip mode: remove all lines containing the svn version id from VERSION
#-----------------------------------------------------------------------------

} else {
   # get all lines from the existing file
   open VER, "<VERSION";
   my @lines = <VER>;
   close VER;
   # truncate file, output all lines that must be kept
   open VER, ">VERSION";
   foreach my $line ( @lines ) {
      if (not $line =~ m/Subversion revision/ ) {
         print VER "$line";
      }
   }
   close VER
}
