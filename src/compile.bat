@echo off

rem Initialize environment variables outside of this script, e.g.:
rem "C:\Program files (x86)\Intel\Composer XE 2013 SP1\bin\compilervars.bat" intel64
rem "C:\Program Files (x86)\IntelSWTools\parallel_studio_xe_2020.1.086\bin\psxevars.bat" intel64

set CNTC_LIB=..\bin\contact_addon_win64.lib
set EXE1=..\bin\test_caddon_win64.exe
set EXE2=..\bin\test_mbench_win64.exe
set EXE3=..\bin\test_varprof_win64.exe
set EXE4=..\bin\test_table_win64.exe
set EXE5=..\bin\caddon_license_win64.exe

rem Since version 20.2, the use of OpenMP is disabled in the library version for Windows
set FORT=ifx
set FFOPTS=/nologo /I. /fpp /O3 /MT /nodebug /Qdiag-disable:10448 /names:lowercase /iface:nomixed_str_len_arg 

echo %FORT% %FFOPTS% test_caddon.f90 %CNTC_LIB% -o %EXE1%
%FORT% %FFOPTS% test_caddon.f90 %CNTC_LIB% -o %EXE1%
%FORT% %FFOPTS% test_mbench.f90 %CNTC_LIB% -o %EXE2%
%FORT% %FFOPTS% test_varprof.f90 %CNTC_LIB% -o %EXE3%
%FORT% %FFOPTS% test_table.f90 %CNTC_LIB% -o %EXE4%
%FORT% %FFOPTS% caddon_license.f90 %CNTC_LIB% -o %EXE5%
