# makefile - compilation of CONTACT software on Linux platform
#
# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# Usage: make PROD=clib_l64
#
# for 32-bits version, first use: 
# Csh, Tcsh: source /opt/intel/bin/ifortvars.csh ia32
# Sh, Bash:  . /opt/intel/bin/ifortvars.sh ia32

# PRODucts are composed as "TARGET"_"PLATF".
#
# Possible TARGET's are
#   "exe"  -  stand-alone version of CONTACT
#   "clib" -  CONTACT library version
#
# Possible PLATForms are
#   "l64"  -  64bit Linux system, using Intel fortran compiler
#   "l32"  -  32bit Linux system, using Intel fortran compiler
#   "w64"  -  64bit Windows system, using Intel fortran compiler
#   "w32"  -  32bit Windows system, using Intel fortran compiler
#
# Separate makefiles are used for Windows (nmake) vs. Linux. These are
# kept similar as much as possible.
#
# Include-files are used for settings for different platforms.

# set default product if no product is specified
ifeq ($(PROD),)
   PROD = exe_l64
endif

# set target and platform for the product that is requested
ifeq ($(PROD), exe_l64)
   TARGET = exe
   PLATF  = linux64
else ifeq ($(PROD), clib_l64)
   TARGET = clib
   PLATF  = linux64
else ifeq ($(PROD), exe_l32)
   TARGET = exe
   PLATF  = linux32
else ifeq ($(PROD), clib_l32)
   TARGET = clib
   PLATF  = linux32
endif

# include compiler settings, compiler flags etc. for the requested platform
include make_inc.linux

# generic parts of the makefile:
.SUFFIXES:
.SUFFIXES:	.f .f90 .c .cpp .o .obj .mod

# main targets of the makefile
ifeq ($(TARGET), exe)
all:	$(STND_EXE)
else ifeq ($(TARGET), clib)
all:	$(CLIB_DLL) $(CLIB_EXE1) $(CLIB_EXE2) $(CLIB_EXE3) $(CLIB_EXE4) $(CLIB_EXE5) $(CLIB_EXE6)
else ifeq ($(TARGET),)
all:
	@echo "ERROR: invalid product PROD=$(PROD)."
endif

# modules and object files needed for the main targets;
# $(P) == object file path, separate directories per target & platform
MODS_ALL = $(P)m_aijpj.$O           $(P)m_blas.$O            $(P)m_bspline.$O         \
           $(P)m_bspline_def.$O     $(P)m_bspline_get.$O     $(P)m_bspline_make.$O    \
           $(P)m_errormsg.$O        $(P)m_friclaw.$O         $(P)m_global_data.$O     \
           $(P)m_gridfunc.$O        $(P)m_gridfunc_sub.$O    $(P)m_grids.$O           \
           $(P)m_hertz.$O           $(P)m_hierarch_data.$O   $(P)m_inflcf.$O          \
           $(P)m_interp.$O          $(P)m_interp_gf.$O       $(P)m_interp_sub.$O      \
           $(P)m_leadedge.$O        $(P)m_licensing.$O                                \
           $(P)m_markers.$O         $(P)m_nonhz.$O           $(P)m_print_output.$O    \
           $(P)m_profile.$O         $(P)m_ptrarray.$O        $(P)m_readline.$O        \
           $(P)m_scontc.$O          $(P)m_sdis.$O            $(P)m_sinput.$O          \
           $(P)m_snorm.$O           $(P)m_solvpn.$O          $(P)m_solvpt.$O          \
           $(P)m_soutpt.$O          $(P)m_spline.$O          $(P)m_spline_def.$O      \
           $(P)m_spline_get.$O      $(P)m_spline_make.$O     $(P)m_stang.$O           \
           $(P)m_subsurf.$O         $(P)m_temperature.$O     $(P)m_timers_contact.$O  \
           $(P)m_timings.$O         $(P)m_varprof.$O         $(P)m_visc.$O            \
           $(P)m_wrprof.$O          $(P)m_wrprof_data.$O     $(P)m_wr_brentmeth.$O    \
           $(P)m_wr_input.$O        $(P)m_wr_locatecp.$O     $(P)m_wr_profiles.$O     \
           $(P)m_wr_output.$O       $(P)m_wr_rigslip.$O      $(P)m_wr_solvecp.$O      \
           $(P)m_wr_totforce.$O     $(P)m_wr_undefdist.$O    $(P)mkl_dfti.$O          
MODS_CLIB= $(P)m_caddon_data.$O
OBJS_ALL = 

# make target for the CONTACT stand-alone executable

$(STND_EXE):	contact.f90 VERSION $(MODS_ALL) $(OBJS_ALL)
	$(F90) $(FFOPTS) contact.f90 $(MODS_ALL) $(OBJS_ALL) $(MKLLIB) $(CRTLIB) -o $(STND_EXE) $(LICLIB)
#	strip --strip-all $(STND_EXE)

# make targets for CONTACT library dll and test programs, with OpenMP

$(CLIB_DLL):	contact_addon.f90 VERSION contact_addon.ifc $(MODS_ALL) $(MODS_CLIB) $(OBJS_ALL)
	$(F90) $(FFOPTS) $(DLLOPT) contact_addon.f90 $(MODS_ALL) $(MODS_CLIB) $(OBJS_ALL) $(MKLLIB) $(CRTLIB) -o $(CLIB_DLL) $(LICLIB)

$(CLIB_EXE1):	test_caddon.f90 $(CLIB_DLL)
	$(F90) $(FFOPTS) -assume nounderscore test_caddon.f90 $(CLIB_LIB) $(LIB32) -o $(CLIB_EXE1)

$(CLIB_EXE2):	test_mbench.f90 $(CLIB_DLL)
	$(F90) $(FFOPTS) -assume nounderscore test_mbench.f90 $(CLIB_LIB) $(LIB32) -o $(CLIB_EXE2)

$(CLIB_EXE3):	test_varprof.f90 $(CLIB_DLL)
	$(F90) $(FFOPTS) -assume nounderscore test_varprof.f90 $(CLIB_LIB) $(LIB32) -o $(CLIB_EXE3)

$(CLIB_EXE4):	test_table.f90 $(CLIB_DLL)
	$(F90) $(FFOPTS) -assume nounderscore test_table.f90 $(CLIB_LIB) $(LIB32) -o $(CLIB_EXE4)

$(CLIB_EXE5):	test_memory.f90 $(CLIB_DLL)
	$(F90) $(FFOPTS) -assume nounderscore test_memory.f90 $(CLIB_LIB) $(LIB32) -o $(CLIB_EXE5)

$(CLIB_EXE6):	usetab_table.f90 $(CLIB_DLL)
	$(F90) $(FFOPTS) -assume nounderscore usetab_table.f90 $(CLIB_LIB) $(LIB32) -o $(CLIB_EXE6)

# generic make targets for compiling f90, f and C source files into objects

$(P)%.mod $(P)%.o:	%.f90
	$(F90) -c $(FFOPTS) $< -o $(P)$*.o

$(P)%.o:	%.f
	$(F90) -c $(FFOPTS) $< -o $@

$(P)%.o:	%.c
	$(CC) -c $(CCOPTS) $< -o $@

# make targets for cleaning up

ifeq ($(TARGET), exe)
clean:
	$(RM) $(P)*.$O $(P)*.mod $(P)*.smod $(STND_EXE)
else ifeq ($(TARGET), clib)
clean:
	$(RM) $(P)*.$O $(P)*.mod $(P)*.smod $(CLIB_DLL) $(MANIF) $(EXPFILE)
	$(RM) $(CLIB_EXE1) $(CLIB_EXE2) $(CLIB_EXE3) $(CLIB_EXE4) $(CLIB_EXE5) $(CLIB_EXE6)
else ifeq ($(TARGET),)
clean:	all
endif

# Dependencies are listed automatically using modules_used.py:

include modules.incl

# additional dependencies/rules per file:

$(P)m_readline.mod      $(P)m_readline.$O:        m_readline.f90
	$(F90) -c $(FFOPTS) -assume bscc $< -o $*.o
