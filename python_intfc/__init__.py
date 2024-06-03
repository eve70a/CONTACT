#------------------------------------------------------------------------------------------------------------
# __init__.py  -  providing a Python interface to the CONTACT library version (trunk).
#
# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
#------------------------------------------------------------------------------------------------------------

import os
import sys
from   platform import system
from   ctypes   import cdll
from   .contact_addon_headers import contact_addon_headers

#------------------------------------------------------------------------------------------------------------

# TODO: unload if loaded before
# TODO: error checking, informative messages

print('Loading the package "python_intfc", system=', system())

my_path  = os.path.dirname( __file__ )
dll_path = os.path.normpath( os.path.join( my_path, '..', 'bin') )

print('Using CONTACT library from', dll_path)
# print('Current working dir is',os.getcwd())

if (system() == 'Linux'):
    dll_filename = os.path.join( dll_path, 'contact_addon_linux64.so' )
    # TODO: check for dir in non-empty LD_LIBRARY_PATH
    if (not os.environ['LD_LIBRARY_PATH']):
        print('You may need to add CONTACTDIR/bin/linux64 to LD_LIBRARY_PATH')
else:
    # extend path for location of additional dlls, cf. LD_LIBRARY_PATH on Linux
    # os.environ['PATH'] += dll_path
    dll_filename = os.path.join( dll_path, 'contact_addon_win64.dll' )

# load library, create ctypes object with function handles

# print('Loading library', dll_filename)
if (not os.path.exists(dll_filename)):
    print('Error: cannot find library file "%s"' % dll_filename)
    exit(1)
cntc_dll = cdll.LoadLibrary(dll_filename)

# define prototypes of the library functions

cntc_dll = contact_addon_headers( cntc_dll )

# list the functions of python_intfc to be exported to users

# from .file_cntc_somefunction          import func_cntc_somefunction
from .cntc_calculate                    import cntc_calculate                   as calculate
from .cntc_closelibrary                 import cntc_closelibrary                as closelibrary
from .cntc_finalize                     import cntc_finalize                    as finalize
from .cntc_finalizelast                 import cntc_finalizelast                as finalizelast
from .cntc_getcalculationtime           import cntc_getcalculationtime          as getcalculationtime
from .cntc_getcontactforces             import cntc_getcontactforces            as getcontactforces
from .cntc_getcontactlocation           import cntc_getcontactlocation          as getcontactlocation
from .cntc_getcontactpatchareas         import cntc_getcontactpatchareas        as getcontactpatchareas
from .cntc_getcreepages                 import cntc_getcreepages                as getcreepages
from .cntc_getdisplacements             import cntc_getdisplacements            as getdisplacements
from .cntc_getelementdivision           import cntc_getelementdivision          as getelementdivision
from .cntc_getfielddata                 import cntc_getfielddata                as getfielddata
from .cntc_getglobalforces              import cntc_getglobalforces             as getglobalforces
from .cntc_getgriddiscretization        import cntc_getgriddiscretization       as getgriddiscretization
from .cntc_gethertzcontact              import cntc_gethertzcontact             as gethertzcontact
from .cntc_getmagicnumbers              import cntc_getmagicnumbers             as getmagicnumbers
from .cntc_getmaximumpressure           import cntc_getmaximumpressure          as getmaximumpressure
from .cntc_getmaximumtemperature        import cntc_getmaximumtemperature       as getmaximumtemperature
from .cntc_getmaximumtraction           import cntc_getmaximumtraction          as getmaximumtraction
from .cntc_getmicroslip                 import cntc_getmicroslip                as getmicroslip
from .cntc_getnumcontactpatches         import cntc_getnumcontactpatches        as getnumcontactpatches
from .cntc_getnumelements               import cntc_getnumelements              as getnumelements
from .cntc_getpenetration               import cntc_getpenetration              as getpenetration
from .cntc_getpotcontact                import cntc_getpotcontact               as getpotcontact
from .cntc_getprofilevalues             import cntc_getprofilevalues            as getprofilevalues
from .cntc_getreferencevelocity         import cntc_getreferencevelocity        as getreferencevelocity
from .cntc_getsensitivities             import cntc_getsensitivities            as getsensitivities
from .cntc_gettractions                 import cntc_gettractions                as gettractions
from .cntc_getwheelsetposition          import cntc_getwheelsetposition         as getwheelsetposition
from .cntc_getwheelsetvelocity          import cntc_getwheelsetvelocity         as getwheelsetvelocity
from .cntc_initlibrary                  import cntc_initlibrary                 as initlibrary
from .cntc_initialize                   import cntc_initialize                  as initialize
from .cntc_readinpfile                  import cntc_readinpfile                 as readinpfile
from .cntc_resetcalculationtime         import cntc_resetcalculationtime        as resetcalculationtime
from .cntc_setcreepages                 import cntc_setcreepages                as setcreepages
from .cntc_setextrarigidslip            import cntc_setextrarigidslip           as setextrarigidslip
from .cntc_setflags                     import cntc_setflags                    as setflags
from .cntc_setfrictionmethod            import cntc_setfrictionmethod           as setfrictionmethod    
from .cntc_setglobalflags               import cntc_setglobalflags              as setglobalflags
from .cntc_sethertzcontact              import cntc_sethertzcontact             as sethertzcontact
from .cntc_setinterfaciallayer          import cntc_setinterfaciallayer         as setinterfaciallayer
from .cntc_setmaterialproperties        import cntc_setmaterialproperties       as setmaterialproperties
from .cntc_setmaterialparameters        import cntc_setmaterialparameters       as setmaterialparameters
from .cntc_setmetadata                  import cntc_setmetadata                 as setmetadata
from .cntc_setnormalforce               import cntc_setnormalforce              as setnormalforce
from .cntc_setpenetration               import cntc_setpenetration              as setpenetration
from .cntc_setpotcontact                import cntc_setpotcontact               as setpotcontact         
from .cntc_setprofileinputfname         import cntc_setprofileinputfname        as setprofileinputfname   
from .cntc_setprofileinputvalues        import cntc_setprofileinputvalues       as setprofileinputvalues
from .cntc_setreferencevelocity         import cntc_setreferencevelocity        as setreferencevelocity
from .cntc_setrollingstepsize           import cntc_setrollingstepsize          as setrollingstepsize    
from .cntc_setsolverflags               import cntc_setsolverflags              as setsolverflags
from .cntc_settangentialforces          import cntc_settangentialforces         as settangentialforces
from .cntc_settemperaturedata           import cntc_settemperaturedata          as settemperaturedata
from .cntc_settimestep                  import cntc_settimestep                 as settimestep
from .cntc_settrackdimensions           import cntc_settrackdimensions          as settrackdimensions    
from .cntc_setundeformeddistc           import cntc_setundeformeddistc          as setundeformeddistc
from .cntc_setverticalforce             import cntc_setverticalforce            as setverticalforce     
from .cntc_setwheelsetdimensions        import cntc_setwheelsetdimensions       as setwheelsetdimensions
from .cntc_setwheelsetflexibility       import cntc_setwheelsetflexibility      as setwheelsetflexibility
from .cntc_setwheelsetposition          import cntc_setwheelsetposition         as setwheelsetposition  
from .cntc_setwheelsetvelocity          import cntc_setwheelsetvelocity         as setwheelsetvelocity  
from .subs_addblock                     import subs_addblock                    as subs_addblock
from .subs_calculate                    import subs_calculate                   as subs_calculate
from .subs_getblocksize                 import subs_getblocksize                as subs_getblocksize
from .subs_getresults                   import subs_getresults                  as subs_getresults
