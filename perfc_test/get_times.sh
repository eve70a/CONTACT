#!/bin/sh

export show_norm=1
export show_tang=1
export show_mbench=1
export show_switch=1
export show_subsurf=1
export show_2dcases=1
export show_parall=1
export show_clib=1

export show_missing=1
export show_failed=1
export idebug=1

#
# Function for printing a separator between different tests
sub_header()
{
   echo " ------------------------------------------------------------------------"
}

#
# Function for retrieving various data-items from a CONTACT output-file
#     in:   $file   - full path-name of CONTACT output-file
#     out:  $ver, $itcg, $itgs, etc.
sub_getdata()
{
   if [ $idebug -ge 5 ]; then
      echo "Starting sub_getdata for file=$file."
   fi
   if [ ! -f "$file" ]; then
      test=`echo $file | sed 's/\.out//g'`
      if [ $show_missing -gt 0 -o $idebug -ge 5 ]; then
         echo "$test: No output for this experiment."
      fi
      is_ok=-1
   else
      is_ok=1
      # get version number of CONTACT executable
      ver=`grep 'Subversion revision' $file| sed 's/)//g' | awk '{print $4}'`
      if [ -z "$ver" ]; then
         ver=`grep Version: $file | awk '{print $5}'`
      fi
      # determine whether the run used module 1 (w/r contact) or not
      is_wr=`grep -c 'Wheel-rail cases' $file`
      if [ $idebug -ge 2 ]; then
         echo "is_wr=$is_wr."
      fi
      # get number of cases performed and total number of Panag. outer iterations
      if [ $is_wr -ge 1 ]; then
         ncase=`grep 'Wheel-rail cases' $file | awk '{print $4}'`
      else
         ncase=`grep 'Panag.process' $file | awk '{print $3}'`
      fi
      nout=`grep 'Algorithm Norm' $file | awk '{print $4}'`
      if [ -z "$ncase" -a -z "$nout" ]; then
         ncase=1
         if [ $show_failed -gt 0 -o $idebug -ge 5 ]; then
            echo "### $test  : No timings for this experiment, run failed?"
         fi
         is_ok=-1
      elif [ -z "$ncase" ]; then
         ncase=$nout
      fi
      # get total number of CG and GS iterations of all cases
      itcg=`grep 'Norm.*ItCG=' $file | awk 'BEGIN {sum=0} {sum += $NF} END {print sum}'`
      itgs=`grep -e 'Tang.*ItGS=' -e 'Tang.*ItCG=' -e 'Tang.*ItGD=' $file | \
                                awk 'BEGIN {sum=0} {sum += $NF} END {print sum}'`
      # get number of elements in contact area for final case (uses O=1)
      ncon_avg=`grep -A1 NCON $file | grep -v -e NCON -e -- | \
                awk 'BEGIN {tot=0; num=0} {tot += $2; num++} END {print tot/(num<1 ? 1 : num)}'`
      ncon=`grep -A1 NCON $file | tail -1 | awk '{print $2}'`
      if [ -z "$ncon" ]; then
         ncon=-
      fi
      nadh=`grep -A1 NADH $file | tail -1 | awk '{print $3}'`
      if [ -z "$nadh" ]; then
         nadh=-
      fi
      nslp=`grep -A1 NSLIP $file | tail -1 | awk '{print $4}'`
      if [ -z "$nslp" ]; then
         nslp=-
      fi
      # display overview of data collected so far
      if [ $idebug -ge 5 ]; then
         echo "ncase=$ncase, nout=$nout, itcg=$itcg, itgs=$itgs."
         echo "ncon=$ncon_avg, $ncon, nadh=$nadh, nslp=$nslp."
      fi
      # get wall-clock calculation times
      wall_tot=`grep '|Total' $file | awk '{print $9}'`
      wall_geom=`grep '|Geometric analysis' $file | awk '{print $10}'`
      wall_panag=`grep '|Panag.process' $file | awk '{print $9}'`
      wall_norm=`grep '|Algorithm Norm' $file | awk '{print $10}'`
      wall_tang=`grep '|Algorithm Tang' $file | awk '{print $10}'`
      wall_subs=`grep '|Subsurface points' $file | awk '{print $10}'`
      if [ -z "$wall_tot" ]; then
         wall_tot=-1
      fi
      if [ -z "$wall_geom" ]; then
         wall_geom=-1
      fi
      if [ -z "$wall_panag" ]; then
         wall_panag=-1
      fi
      if [ -z "$wall_norm" ]; then
         wall_norm=-1
      fi
      if [ -z "$wall_tang" ]; then
         wall_tang=-1
      fi
      if [ -z "$wall_subs" ]; then
         wall_subs=-1
      fi
      # display calculation times when requested
      if [ $idebug -ge 5 ]; then
         echo "wall_tot=$wall_tot, wall_geom=$wall_geom, wall_panag=$wall_panag,"
         echo "wall_norm=$wall_norm, wall_tang=$wall_tang, wall_subs=$wall_subs."
      fi
      # in sequential runs, copy wall-times to seq-times
      nthrd=`grep 'Parallel run using' $file | awk '{print $4}'`
      if [ -z "$nthrd" ]; then
         nthrd=1
         seq_tot=$wall_tot ; seq_norm=$wall_norm ; seq_tang=$wall_tang
         if [ $idebug -ge 5 ]; then
            echo "Sequential run copy tot/norm/tang to seq_..."
         fi
      fi
   fi
}

#
# Function for printing wall-clock time for the normal problem.
sub_print_norm()
{
   if [ $is_ok -gt 0 ]; then
      # note: how to round in awk? print %d rounds downwards, add 0.5
      echo "$test $ncon $itcg $wall_tot $ver" | \
      awk -F" " '{ if ($4 > 50) {
                      printf "%-14s: ncon=%6d, ItCG=%5d, tot.wall=%6d (ver=%5s)\n", $1, $2, $3, $4+0.5, $5 ;
                   } else {
                      printf "%-14s: ncon=%6d, ItCG=%5d, tot.wall=%6.1f (ver=%5s)\n", $1, $2, $3, $4, $5 ;
                   }
                 }' | \
      sed 's/\([0-9]\),\([0-9]\)/\1.\2/g'
   fi
}

#
# Function for printing parallel speedup for the normal problem.
# norm_parall_1s (8 tr): ItCG=   19, seq.norm=   0.1, par=   0.5, fac=  0.2 (1181)
sub_print_norm_parall()
{
   if [ $is_ok -gt 0 ]; then
      # note: how to round in awk? print %d rounds downwards, add 0.5
      echo "$test $nthrd $itcg $seq_norm $wall_norm $ver" | \
      awk -F" " '{ if ($5 > 50) {
                      printf "%-14s (%1d tr): ItCG=%5d, seq.norm=%6d, par=%6d, fac=%5.1f (%4s)\n", $1, $2, $3, $4+0.5, $5+0.5, $4/$5, $6 ;
                   } else {
                      printf "%-14s (%1d tr): ItCG=%5d, seq.norm=%6.1f, par=%6.1f, fac=%5.1f (%4s)\n", $1, $2, $3, $4, $5, $4/$5, $6 ;
                   }
                 }' | \
      sed 's/\([0-9]\),\([0-9]\)/\1.\2/g'
   fi
}

#
# Function for printing wall-clock time for the tangential problem.
# tang_problm_1c : nslp=  1872, ItGS=   48, wall.tang=   4.4 (ver= 1181)
sub_print_tang()
{
   if [ $is_ok -gt 0 ]; then
      # note: how to round in awk? print %d rounds downwards, add 0.5
      echo "$test $nslp $itgs $wall_tang $ver" | \
      awk -F" " '{ if ($4 > 50) {
                      printf "%-14s: nslp=%6d, ItGS=%5d, wall.tang=%6d (ver=%5s)\n", $1, $2, $3, $4+0.5, $5 ;
                   } else {
                      printf "%-14s: nslp=%6d, ItGS=%5d, wall.tang=%6.1f (ver=%5s)\n", $1, $2, $3, $4, $5 ;
                   }
                 }' | \
      sed 's/\([0-9]\),\([0-9]\)/\1.\2/g'
   fi
}

#
# Function for printing parallel speedup for the tangential problem.
# tang_parall_1s (8 tr): ItGS=   51, seq.tang=   0.3, par=   0.3, fac=  1.0 (1181)
sub_print_tang_parall()
{
   if [ $is_ok -gt 0 ]; then
      # note: how to round in awk? print %d rounds downwards, add 0.5
      echo "$test $nthrd $itgs $seq_tang $wall_tang $ver" | \
      awk -F" " '{ if ($5 > 50) {
                      printf "%-14s (%1d tr): ItGS=%5d, seq.tang=%6d, par=%6d, fac=%5.1f (%4s)\n", $1, $2, $3, $4+0.5, $5+0.5, $4/$5, $6 ;
                   } else {
                      printf "%-14s (%1d tr): ItGS=%5d, seq.tang=%6.1f, par=%6.1f, fac=%5.1f (%4s)\n", $1, $2, $3, $4, $5, $4/$5, $6 ;
                   }
                 }' | \
      sed 's/\([0-9]\),\([0-9]\)/\1.\2/g'
   fi
}

#
# Function for printing wall-clock time for making tables
# table_11x11 : ncase= 3221, ncon=    0, tot.wall=   7.9,   2.5 ms/cs (ver= 1181)
sub_print_table()
{
   if [ $is_ok -gt 0 ]; then
      # note: how to round in awk? print %d rounds downwards, add 0.5
      echo "$test $ncase $ncon $wall_tot $ver" | \
      awk -F" " '{ if ($4 > 50) {
                      printf "%-13s: ncase=%5d, ncon=%5d, tot.wall=%6d, %5d ms/cs (ver=%5s)\n", $1, $2, $3, $4+0.5, ($4*1000)/$2, $5 ;
                   } else {
                      printf "%-13s: ncase=%5d, ncon=%5d, tot.wall=%6.1f, %5.1f ms/cs (ver=%5s)\n", $1, $2, $3, $4, ($4*1000)/$2, $5 ;
                   }
                 }' | \
      sed 's/\([0-9]\),\([0-9]\)/\1.\2/g'
   fi
}

#
# Function for printing wall-clock time for making tables in parallel runs
# test_clib_par11: ncase= 3220, seq.wall=   7.7, par=  1.1, fac=  7.0 (ver= 1216)
sub_print_table_parall()
{
   if [ $is_ok -gt 0 ]; then
      # note: how to round in awk? print %d rounds downwards, add 0.5
      echo "$test $ncase $seq_tot $wall_tot $ver" | \
      awk -F" " '{ if ($4 > 50) {
                      printf "%-13s: ncase=%5d, seq.wall=%6d, par=%6d, fac=%4.1f (ver=%5s)\n", $1, $2, $3, $4+0.5, $3/$4, $5 ;
                   } else {
                      printf "%-13s: ncase=%5d, seq.wall=%6.1f, par=%6.1f, fac=%4.1f (ver=%5s)\n", $1, $2, $3, $4, $3/$4, $5 ;
                   }
                 }' | \
      sed 's/\([0-9]\),\([0-9]\)/\1.\2/g'
   fi
}

#
# Function for printing wall-clock time for wheel-rail contact calculations
sub_print_mbench()
{
   if [ $is_ok -gt 0 ]; then
      # note: how to round in awk? print %d rounds downwards, add 0.5
      echo "$test $ncon_avg $wall_geom $wall_panag $ver" | \
      awk -F" " '{ if ($4 > 50) {
                      printf "%-22s: ncon=%6d, t.geom=%6.1f  t.panag=%6d   (ver=%5s)\n", $1, $2, $3, $4+0.5, $5 ;
                   } else {                                                     
                      printf "%-22s: ncon=%6.1f, t.geom=%6.1f  t.panag=%6.1f   (ver=%5s)\n", $1, $2, $3, $4, $5 ;
                   }
                 }' | \
      sed 's/\([0-9]\),\([0-9]\)/\1.\2/g'
   fi
}

# Function for printing wall-clock time for subsurface stress calculations
sub_print_subsurf()
{
   if [ $is_ok -gt 0 ]; then
      # note: how to round in awk? print %d rounds downwards, add 0.5
      echo "$test $ncon $nout $wall_tot $wall_subs $ver" | \
      awk -F" " '{ if ($5 < 0) { # print overview
                      if ($4 > 50) {
                         printf "%-16s: ncon=%5d, nout=%5d,  tot.wall=%6d   (ver=%5s)\n", $1, $2, $3, $4+0.5, $6 ;
                      } else {
                         printf "%-16s: ncon=%5d, nout=%5d,  tot.wall=%6.1f   (ver=%5s)\n", $1, $2, $3, $4, $6 ;
                      }
                   } else { # print time subsurface
                      if ($5 > 50) {
                         printf "%-16s:                         wall.subs=%8d (ver=%5s)\n", $1, $5+0.5, $6 ;
                      } else {
                         printf "%-16s:                         wall.subs=%8.1f (ver=%5s)\n", $1, $5, $6 ;
                      }
                   }
                 }' | \
      sed 's/\([0-9]\),\([0-9]\)/\1.\2/g'
   fi
}

#
# Function for printing iterations and wall-clock time for 2D cases
sub_print_2d_norm()
{
   if [ $is_ok -gt 0 ]; then
      # note: how to round in awk? print %d rounds downwards, add 0.5
      echo "$test $ncase $ncon $itcg $wall_norm $nslp $itgs $wall_tang $ver" | \
      awk -F" " '{ if ($5 > 50) {
                      printf "%-14s: ncon=%6d, ItCG=%6.1f, wall.norm=%6d (ver=%5s)\n", $1, $3, $4/$2, $5+0.5, $9 ;
                   } else {
                      printf "%-14s: ncon=%6d, ItCG=%6.1f, wall.norm=%6.1f (ver=%5s)\n", $1, $3, $4/$2, $5, $9 ;
                   }
                 }' | \
      sed 's/\([0-9]\),\([0-9]\)/\1.\2/g'
   fi
}

#
# Function for printing iterations and wall-clock time for 2D cases
sub_print_2d_tang()
{
   if [ $is_ok -gt 0 ]; then
      # note: how to round in awk? print %d rounds downwards, add 0.5
      echo "$test $ncase $ncon $itcg $wall_norm $nslp $itgs $wall_tang $ver" | \
      awk -F" " '{ if ($8 > 50) {
                      printf "%-14s: nslp=%6d, ItGS=%6.1f, wall.tang=%6d (ver=%5s)\n", $1, $6, $7/$2, $8+0.5, $9 ;
                   } else {
                      printf "%-14s: nslp=%6d, ItGS=%6.1f, wall.tang=%6.1f (ver=%5s)\n", $1, $6, $7/$2, $8, $9 ;
                   }
                 }' | \
      sed 's/\([0-9]\),\([0-9]\)/\1.\2/g'
   fi
}

###########################################################################
# "Main program"
###########################################################################

#
# ensure that period "." is used as decimal split instead of comma ","
export LC_ALL=C

sub_header
grep -i -e 'performance test' -e 'Running at host' run_perfc.log | \
    awk 'BEGIN { days["jan"]=31; days["feb"]=28; days["mar"]=31;
                 days["apr"]=30; days["may"]=31; days["jun"]=30;
                 days["jul"]=31; days["aug"]=31; days["sep"]=30;
                 days["oct"]=31; days["nov"]=30; days["dec"]=31;
                 days[  1 ]=31; days[  2 ]=28; days[  3 ]=31;
                 days[  4 ]=30; days[  5 ]=31; days[  6 ]=30;
                 days[  7 ]=31; days[  8 ]=31; days[  9 ]=30;
                 days[ 10 ]=31; days[ 11 ]=30; days[ 12 ]=31;
                 found2=0;
               }
         {
            print $0;
            if (NF==6) {
                # Starting performance test at 01-08-2011 13:30:54
               if (/Starting/) { month1=substr($5,4,2); day1=substr($5,1,2);
                                 hr1=substr($6,1,2); mn1=substr($6,4,2); 
                                 sec1=substr($6,7,2); 
                               }
               if (/complete/) { month2=substr($5,4,2); day2=substr($5,1,2);
                                 hr2=substr($6,1,2); mn2=substr($6,4,2); 
                                 sec2=substr($6,7,2); found2=1; 
                               }
            } else {
               if (/Starting/) { month1=$6; day1=$7; hr1=substr($8,1,2);
                                 mn1=substr($8,4,2); sec1=substr($8,7,2); 
                                 printf "month1= %d, %d\n", month1, day1;
                               }
               if (/complete/) { found2=1; 
                                 month2=$6; day2=$7; hr2=substr($8,1,2);
                                 mn2=substr($8,4,2); sec2=substr($8,7,2); }
            }
         }
         END {
              if (found2==1) {
                 if (sec1-sec2>0) { mn2=mn2-1; sec2=sec2+60 }
                 if (mn1-mn2>0) { hr2=hr2-1; mn2=mn2+60 }
                 if (hr1-hr2>0) { day2=day2-1; hr2=hr2+24 }
                 if (day1-day2>0) { day2=day2+days[month1] }
                 printf " Total %d days, %d hr, %d min, %d sec\n", 
                             day2-day1, hr2-hr1, mn2-mn1, sec2-sec1;
              }
         }'
sub_header

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# retrieve and print times for the normal problem

if [ $show_norm -gt 0 ]; then
   for t in 1p 2p 4p 8p 1f 2f 4f ; do
      test=norm_problm_$t ; file=$test.out ; sub_getdata ; sub_print_norm
   done
   test=norm_potcon_1p ; file=$test.out ; sub_getdata ; sub_print_norm
   sub_header
fi

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# retrieve and print speedup for the normal problem

if [ $show_parall -gt 0 -a $show_norm -gt 0 ]; then
   for t in 1s 2s 4s 8s ; do
      file=tang_problm_$t.out ; sub_getdata
      if [ $is_ok -gt 0 ]; then
         file=tang_parall_$t.out ; sub_getdata
      fi
      test=norm_parall_$t ; sub_print_norm_parall
   done
   sub_header
fi

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# retrieve and print times for the tangential problem

if [ $show_tang -gt 0 ]; then
   for t in 1s 2s 4s 8s 1c 2c 4c 8c 1f 2f 4f ; do
      test=tang_problm_$t ; file=$test.out ; sub_getdata ; sub_print_tang
   done
   sub_header

   test=tang_potcon_1c ; file=$test.out ; sub_getdata ; sub_print_tang
   test=tang_cnvxgs_1c ; file=$test.out ; sub_getdata ; sub_print_tang
   sub_header
fi

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# retrieve and print speedup for the tangential problem

if [ $show_parall -gt 0 -a $show_tang -gt 0 ]; then
   for t in 1s 2s 4s 8s 1c 2c 4c ; do
      file=tang_problm_$t.out ; sub_getdata
      if [ $is_ok -gt 0 ]; then
         file=tang_parall_$t.out ; sub_getdata
      fi
      test=tang_parall_$t ; sub_print_tang_parall
   done
   sub_header
fi

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# retrieve and print times for making tables using the CONTACT library

if [ $show_clib -gt 0 ]; then

   # print data for sequential runs
   for m in 11 15 19 25 35; do
      test=test_clib_seq$m ; file=$test.out ; sub_getdata ; sub_print_table
   done
   sub_header

   # print data for parallel runs
   for m in 11 15 19 25 35; do
      test=test_clib_seq$m ; file=$test.out ; sub_getdata ;
      test=test_clib_par$m ; file=$test.out ; sub_getdata ; sub_print_table_parall
   done
   sub_header
fi

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# retrieve and print times for wheel-rail contact calculations

if [ $show_mbench -gt 0 ]; then
   test=mbench_20kn_coarseprf ; file=$test.out ; sub_getdata ; sub_print_mbench
   test=mbench_20kn_dx05_zpos ; file=$test.out ; sub_getdata ; sub_print_mbench
   test=mbench_20kn_dx02_zpos ; file=$test.out ; sub_getdata ; sub_print_mbench
   test=mbench_20kn_dx02_fz   ; file=$test.out ; sub_getdata ; sub_print_mbench
   sub_header
   test=mbench_100kn_coarseprf ; file=$test.out ; sub_getdata ; sub_print_mbench
   test=mbench_100kn_dx05_zpos ; file=$test.out ; sub_getdata ; sub_print_mbench
   test=mbench_100kn_dx02_zpos ; file=$test.out ; sub_getdata ; sub_print_mbench
   test=mbench_100kn_dx02_fz   ; file=$test.out ; sub_getdata ; sub_print_mbench
   sub_header
fi

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# retrieve and print times for switch-and-crossing calculations

if [ $show_switch -gt 0 ]; then
   test=mbench_intrup ; file=$test.out ; sub_getdata ; sub_print_mbench
   test=mbench_locus  ; file=$test.out ; sub_getdata ; sub_print_mbench
   test=mbench_brute  ; file=$test.out ; sub_getdata ; sub_print_mbench
   sub_header
fi

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# retrieve and print times for subsurface stress calculations

if [ $show_subsurf -gt 0 ]; then
   test=spence35_nosubs ; file=$test.out ; sub_getdata ; sub_print_subsurf
   test=spence35_135pt  ; file=$test.out ; sub_getdata ; sub_print_subsurf
   test=spence35_2025pt ; file=$test.out ; sub_getdata ; sub_print_subsurf
   test=spence71_nosubs ; file=$test.out ; sub_getdata ; sub_print_subsurf
   test=spence71_273pt  ; file=$test.out ; sub_getdata ; sub_print_subsurf
   test=spence71_8281pt ; file=$test.out ; sub_getdata ; sub_print_subsurf
   sub_header
fi

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# retrieve and print times for the 2D testcases

if [ $show_2dcases -gt 0 ]; then
   for test in carter2d_100 carter2d_800 veldep_fric visc_cylindr ; do
      file=$test.out ; sub_getdata ; sub_print_2d_norm
   done
   sub_header
   for test in carter2d_100 carter2d_800 veldep_fric visc_cylindr ; do
      file=$test.out ; sub_getdata ; sub_print_2d_tang
   done
   sub_header
fi

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
