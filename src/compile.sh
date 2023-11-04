#!/bin/sh

DLL=../bin/contact_addon_linux64.so
EXE1=../bin/test_caddon_linux64
EXE2=../bin/test_mbench_linux64
EXE3=../bin/test_varprof_linux64
EXE4=../bin/test_table_linux64
EXE5=../bin/caddon_license_linux64

if [ 1 = 1 ]; then
   echo "Using Intel fortran..."
   FFOPTS="-fpp -qopenmp -names lowercase -assume nounderscore -convert big_endian"

   ifort $FFOPTS test_caddon.f90    $DLL -o $EXE1
   ifort $FFOPTS test_mbench.f90    $DLL -o $EXE2
   ifort $FFOPTS test_varprof.f90   $DLL -o $EXE3
   ifort $FFOPTS test_table.f90     $DLL -o $EXE4
   ifort $FFOPTS caddon_license.f90 $DLL -o $EXE5

else
   echo "Using gfortran..."
   FFOPTS="-cpp '-DCNAME_(x)="x"' -fopenmp -fno-underscoring"
   gfortran $FFOPTS test_caddon.f90    $DLL -o $EXE1
   gfortran $FFOPTS test_mbench.f90    $DLL -o $EXE2
   gfortran $FFOPTS test_varprof.f90   $DLL -o $EXE3
   gfortran $FFOPTS test_table.f90     $DLL -o $EXE4
   gfortran $FFOPTS caddon_license.f90 $DLL -o $EXE5
fi
