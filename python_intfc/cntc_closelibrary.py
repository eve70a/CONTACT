
#------------------------------------------------------------------------------------------------------------
# function [ ] = cntc_closelibrary();
#
# clean-up, close files and unload the library from Matlab
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 0: m=*, glob   - no icp needed

import ctypes
from platform              import system
from python_intfc          import cntc_dll

def cntc_closelibrary():

    cntc_dll.cntc_finalizelast()

    if (system()=='Linux'):
        if_body_is_empty=1
        # ctypes.cdll.dlclose(cntc_dll._handle)
    else:
        from ctypes          import WinDLL, wintypes
        kernel32 = ctypes.WinDLL('kernel32', use_last_error=True)
        kernel32.FreeLibrary.argtypes = [wintypes.HMODULE]
        kernel32.FreeLibrary(cntc_dll._handle)

# end function function cntc_closelibrary

#------------------------------------------------------------------------------------------------------------

