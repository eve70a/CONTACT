from ctypes import c_char_p, c_int, c_double, POINTER

############################################################################################################# 
# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

def contact_addon_headers( cntc_dll ):


    # void cntc_initializefirst                     ( int *ifcver,       int *ierror,       int *ioutput,
    #                                          const char* c_outpath, const char* c_expnam, int *len_outpath, 
    #                                                 int *len_expnam );

    cntc_dll.cntc_initializefirst.restype         = None
    cntc_dll.cntc_initializefirst.argtypes        = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int), 
                                                      c_char_p,          c_char_p,          POINTER(c_int),   
                                                      POINTER(c_int) ]

    # void cntc_initialize                          ( int *ire,          int *imodul,       int *ifcver,
    #                                                 int *ierror,   const char* c_outpath, int *len_outpath );

    cntc_dll.cntc_initialize.restype              = None
    cntc_dll.cntc_initialize.argtypes             = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int), 
                                                      POINTER(c_int),    c_char_p,          POINTER(c_int) ]
  
    # void cntc_setglobalflags                      ( int *lenflg,       int *params,       int *values );
  
    cntc_dll.cntc_setglobalflags.restype          = None
    cntc_dll.cntc_setglobalflags.argtypes         = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int) ]
  
    # void cntc_setflags                            ( int *ire,          int *icp,          int *lenflg, 
    #                                                 int *params,       int *values );
  
    cntc_dll.cntc_setflags.restype                = None
    cntc_dll.cntc_setflags.argtypes               = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_int),    POINTER(c_int) ]
  
    # void cntc_setmetadata                         ( int *ire,          int *icp,          int *lenmta, 
    #                                                 int *params,       double *values );
  
    cntc_dll.cntc_setmetadata.restype             = None
    cntc_dll.cntc_setmetadata.argtypes            = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int), 
                                                      POINTER(c_int),    POINTER(c_double) ]
  
    # void cntc_setsolverflags                      ( int *ire,          int *icp,          int *imeth, 
    #                                                 int *nints,        int *iparam,       int *nreals, 
    #                                                 double *rparam );
  
    cntc_dll.cntc_setsolverflags.restype          = None
    cntc_dll.cntc_setsolverflags.argtypes         = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int), 
                                                      POINTER(c_int),    POINTER(c_int),    POINTER(c_int), 
                                                      POINTER(c_double) ]
  
    # void cntc_setmaterialparameters               ( int *ire,          int *icp,          int *m_digit, 
    #                                                 int *nparam,       double *params );
  
    cntc_dll.cntc_setmaterialparameters.restype   = None
    cntc_dll.cntc_setmaterialparameters.argtypes  = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_int),    POINTER(c_double) ]
  
    # void cntc_setmaterialproperties               ( int *ire,          int *icp,          double *g1, 
    #                                                 double *nu1,        double *g2,        double *nu2 );
  
    cntc_dll.cntc_setmaterialproperties.restype   = None
    cntc_dll.cntc_setmaterialproperties.argtypes  = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_double),
                                                      POINTER(c_double), POINTER(c_double), POINTER(c_double) ]
  
    # void cntc_settemperaturedata                  ( int *ire,          int *icp,          int *imeth, 
    #                                                 int *nparam,       double *rparam );
  
    cntc_dll.cntc_settemperaturedata.restype      = None
    cntc_dll.cntc_settemperaturedata.argtypes     = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_int),    POINTER(c_double) ]
  
    # void cntc_settimestep                         ( int *ire,          int *icp,          double *dt );
  
    cntc_dll.cntc_settimestep.restype             = None
    cntc_dll.cntc_settimestep.argtypes            = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_double) ]
  
    # void cntc_setreferencevelocity                ( int *ire,          int *icp,          double *veloc );
  
    cntc_dll.cntc_setreferencevelocity.restype    = None
    cntc_dll.cntc_setreferencevelocity.argtypes   = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_double) ]
  
    # void cntc_setrollingstepsize                  ( int *ire,          int *icp,          double *chi, 
    #                                                 double *dq );
  
    cntc_dll.cntc_setrollingstepsize.restype      = None
    cntc_dll.cntc_setrollingstepsize.argtypes     = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_double),
                                                      POINTER(c_double) ]
  
    # void cntc_setfrictionmethod                   ( int *ire,          int *icp,          int *imeth, 
    #                                                 int *nparam,       double *params );
  
    cntc_dll.cntc_setfrictionmethod.restype       = None
    cntc_dll.cntc_setfrictionmethod.argtypes      = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_int),    POINTER(c_double) ]
  
    # void cntc_setinterfaciallayer                 ( int *ire,          int *icp,          int *imeth, 
    #                                                 int *nparam,       double *params );
  
    cntc_dll.cntc_setinterfaciallayer.restype     = None
    cntc_dll.cntc_setinterfaciallayer.argtypes    = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_int),    POINTER(c_double) ]
  
    # void cntc_sethertzcontact                     ( int *ire,          int *icp,          int *ipotcn, 
    #                                                 int *nparam,       double *rparam );
  
    cntc_dll.cntc_sethertzcontact.restype         = None
    cntc_dll.cntc_sethertzcontact.argtypes        = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_int),    POINTER(c_double) ]
  
    # void cntc_setpotcontact                       ( int *ire,          int *icp,          int *ipotcn, 
    #                                                 int *nparam,       double *rparam );
  
    cntc_dll.cntc_setpotcontact.restype           = None
    cntc_dll.cntc_setpotcontact.argtypes          = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_int),    POINTER(c_double) ]
  
    # void cntc_setverticalforce                    ( int *ire,          double *fz );
  
    cntc_dll.cntc_setverticalforce.restype        = None
    cntc_dll.cntc_setverticalforce.argtypes       = [ POINTER(c_int),    POINTER(c_double) ]
  
    # void cntc_setpenetration                      ( int *ire,          int *icp,          double *pen );
  
    cntc_dll.cntc_setpenetration.restype          = None
    cntc_dll.cntc_setpenetration.argtypes         = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_double) ]
  
    # void cntc_setnormalforce                      ( int *ire,          int *icp,          double *fn );
  
    cntc_dll.cntc_setnormalforce.restype          = None
    cntc_dll.cntc_setnormalforce.argtypes         = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_double) ]
  
    # void cntc_setundeformeddistc                  ( int *ire,          int *icp,          int *ibase, 
    #                                                 int *nparam,       double *prmudf );
  
    cntc_dll.cntc_setundeformeddistc.restype      = None
    cntc_dll.cntc_setundeformeddistc.argtypes     = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_int),    POINTER(c_double) ]
  
    # void cntc_setcreepages                        ( int *ire,          int *icp,          double *vx, 
    #                                                 double *vy,        double *phi );
  
    cntc_dll.cntc_setcreepages.restype            = None
    cntc_dll.cntc_setcreepages.argtypes           = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_double),
                                                      POINTER(c_double), POINTER(c_double) ]
  
    # void cntc_setextrarigidslip                   ( int *ire,          int *icp,          int *lenarr, 
    #                                                 double *wx,        double *wy );
  
    cntc_dll.cntc_setextrarigidslip.restype       = None
    cntc_dll.cntc_setextrarigidslip.argtypes      = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_double), POINTER(c_double) ]
  
    # void cntc_settangentialforces                 ( int *ire,          int *icp,          double *fx, 
    #                                                 double *fy );
  
    cntc_dll.cntc_settangentialforces.restype     = None
    cntc_dll.cntc_settangentialforces.argtypes    = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_double),
                                                      POINTER(c_double) ]
  
    # void cntc_settrackdimensions                  ( int *ire,          int *ztrack,       int *nparam, 
    #                                                 double *rparam );
  
    cntc_dll.cntc_settrackdimensions.restype      = None
    cntc_dll.cntc_settrackdimensions.argtypes     = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_double) ]
  
    # void cntc_settrackdimensions_new              ( int *ire,          int *ztrack,       int *nparam, 
    #                                                 double *rparam );
  
    cntc_dll.cntc_settrackdimensions_new.restype  = None
    cntc_dll.cntc_settrackdimensions_new.argtypes = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_double) ]
  
    # void cntc_settrackdimensions_old              ( int *ire,          int *ztrack,       int *nparam, 
    #                                                 double *rparam );
  
    cntc_dll.cntc_settrackdimensions_old.restype  = None
    cntc_dll.cntc_settrackdimensions_old.argtypes = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_double) ]
  
    # void cntc_setprofileinputfname                ( int *ire,       const char* c_fname,  int *len_fname, 
    #                                                 int *nints,        int *iparam,       int *nreals, 
    #                                                 double *rparam );
  
    cntc_dll.cntc_setprofileinputfname.restype    = None
    cntc_dll.cntc_setprofileinputfname.argtypes   = [ POINTER(c_int),    c_char_p,          POINTER(c_int),
                                                      POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_double) ]
  
    # void cntc_setprofileinputvalues               ( int *ire,          int *npoint,       double *values, 
    #                                                 int *nints,        int *iparam,       int *nreals, 
    #                                                 double *rparam );
  
    cntc_dll.cntc_setprofileinputvalues.restype   = None
    cntc_dll.cntc_setprofileinputvalues.argtypes  = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_double),
                                                      POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_double) ]
  
    # void cntc_setwheelsetdimensions               ( int *ire,          int *ewheel,       int *nparam, 
    #                                                 double *rparam );
  
    cntc_dll.cntc_setwheelsetdimensions.restype   = None
    cntc_dll.cntc_setwheelsetdimensions.argtypes  = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_double) ]
  
    # void cntc_setwheelsetposition                 ( int *ire,          int *imeth,        int *nparam, 
    #                                                 double *rparam );
  
    cntc_dll.cntc_setwheelsetposition.restype     = None
    cntc_dll.cntc_setwheelsetposition.argtypes    = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_double) ]
  
    # void cntc_setwheelsetvelocity                 ( int *ire,          int *imeth,        int *nparam, 
    #                                                 double *rparam );
  
    cntc_dll.cntc_setwheelsetvelocity.restype     = None
    cntc_dll.cntc_setwheelsetvelocity.argtypes    = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_double) ]
  
    # void cntc_setwheelsetflexibility              ( int *ire,          int *imeth,        int *nparam, 
    #                                                 double *rparam );
  
    cntc_dll.cntc_setwheelsetflexibility.restype  = None
    cntc_dll.cntc_setwheelsetflexibility.argtypes = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_double) ]
  
    # void subs_addblock                            ( int *ire,          int *icp,          int *iblk,
    #                                                 int *isubs,        int *nx,           int *ny,
    #                                                 int *nz,           double *xparam,    double *yparam,
    #                                                 double *zparam );
    
    cntc_dll.subs_addblock.restype                = None
    cntc_dll.subs_addblock.argtypes               = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_int),    POINTER(c_double), POINTER(c_double),
                                                      POINTER(c_double) ]
    
    # void cntc_calculate                           ( int *ire,          int *icp,          int *ierror );
    
    cntc_dll.cntc_calculate.restype               = None
    cntc_dll.cntc_calculate.argtypes              = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int) ]
    
    # void subs_calculate                           ( int *ire,          int *icp,          int *ierror );
    
    cntc_dll.subs_calculate.restype               = None
    cntc_dll.subs_calculate.argtypes              = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int) ]
    
    # void cntc_getprofilevalues                    ( int *ire,          int *itask,        int *nints, 
    #                                                 int *iparam,       int *lenarr,       double *values );
    
    cntc_dll.cntc_getprofilevalues.restype        = None
    cntc_dll.cntc_getprofilevalues.argtypes       = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_int),    POINTER(c_int),    POINTER(c_double) ]
    
    # void cntc_getprofilevalues_new                ( int *ire,          int *itask,        int *nints, 
    #                                                 int *iparam,       int *nreals,       double *rparam,
    #                                                 int *lenarr,       double *values );
    
    cntc_dll.cntc_getprofilevalues_new.restype    = None
    cntc_dll.cntc_getprofilevalues_new.argtypes   = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_int),    POINTER(c_int),    POINTER(c_double),
                                                      POINTER(c_int),    POINTER(c_double) ]
    
    # void cntc_getwheelsetposition                 ( int *ire,          int *lenarr,       double *values );
    
    cntc_dll.cntc_getwheelsetposition.restype     = None
    cntc_dll.cntc_getwheelsetposition.argtypes    = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_double) ]
    
    # void cntc_getwheelsetvelocity                 ( int *ire,          int *lenarr,       double *values );
    
    cntc_dll.cntc_getwheelsetvelocity.restype     = None
    cntc_dll.cntc_getwheelsetvelocity.argtypes    = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_double) ]
    
    # void cntc_getnumcontactpatches                ( int *ire,          int *npatch );
    
    cntc_dll.cntc_getnumcontactpatches.restype    = None
    cntc_dll.cntc_getnumcontactpatches.argtypes   = [ POINTER(c_int),    POINTER(c_int) ] 
    
    # void cntc_getcontactlocation                  ( int *ire,          int *icp,          int *lenarr,
    #                                                 double *values );
    
    cntc_dll.cntc_getcontactlocation.restype      = None
    cntc_dll.cntc_getcontactlocation.argtypes     = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_double) ]
    
    # void cntc_getreferencevelocity                ( int *ire,          int *icp,          double *veloc );
    
    cntc_dll.cntc_getreferencevelocity.restype    = None
    cntc_dll.cntc_getreferencevelocity.argtypes   = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_double) ]
    
    # void cntc_gethertzcontact                     ( int *ire,          int *icp,          int *lenarr, 
    #                                                 double *values );
    
    cntc_dll.cntc_gethertzcontact.restype         = None
    cntc_dll.cntc_gethertzcontact.argtypes        = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_double) ]
    
    # void cntc_getnumelements                      ( int *ire,          int *icp,          int *mx, 
    #                                                 int *my );
    
    cntc_dll.cntc_getnumelements.restype          = None
    cntc_dll.cntc_getnumelements.argtypes         = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_int) ]
    
    # void cntc_getgriddiscretization               ( int *ire,          int *icp,          double *dx, 
    #                                                 double *dy );
    
    cntc_dll.cntc_getgriddiscretization.restype   = None
    cntc_dll.cntc_getgriddiscretization.argtypes  = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_double),
                                                      POINTER(c_double) ]
    
    # void cntc_getpotcontact                       ( int *ire,          int *icp,          int *lenarr, 
    #                                                 double *values );
    
    cntc_dll.cntc_getpotcontact.restype           = None
    cntc_dll.cntc_getpotcontact.argtypes          = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_double) ]
    
    # void cntc_getpenetration                      ( int *ire,          int *icp,          double *pen );
    
    cntc_dll.cntc_getpenetration.restype          = None
    cntc_dll.cntc_getpenetration.argtypes         = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_double) ]
    
    # void cntc_getcreepages                        ( int *ire,          int *icp,          double *vx, 
    #                                                 double *vy,        double *phi );
    
    cntc_dll.cntc_getcreepages.restype            = None
    cntc_dll.cntc_getcreepages.argtypes           = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_double),
                                                      POINTER(c_double), POINTER(c_double) ]
    
    # void cntc_getcontactforces                    ( int *ire,          int *icp,          double *fn, 
    #                                                 double *tx,        double *ty,        double *mz );
    
    cntc_dll.cntc_getcontactforces.restype        = None
    cntc_dll.cntc_getcontactforces.argtypes       = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_double),
                                                      POINTER(c_double), POINTER(c_double), POINTER(c_double) ]
    
    # void cntc_getglobalforces                     ( int *ire,          int *icp,          int *lenarr, 
    #                                                 double *values );
    
    cntc_dll.cntc_getglobalforces.restype         = None
    cntc_dll.cntc_getglobalforces.argtypes        = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_double) ]
    
    # void cntc_getcontactpatchareas                ( int *ire,          int *icp,          double *carea, 
    #                                                 double *harea,     double *sarea );
    
    cntc_dll.cntc_getcontactpatchareas.restype    = None
    cntc_dll.cntc_getcontactpatchareas.argtypes   = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_double),
                                                      POINTER(c_double), POINTER(c_double) ]
    
    # void cntc_getelementdivision                  ( int *ire,          int *icp,          int *lenarr,
    #                                                 int *eldiv ); 
    
    cntc_dll.cntc_getelementdivision.restype      = None
    cntc_dll.cntc_getelementdivision.argtypes     = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_int) ]
    
    # void cntc_getmaximumpressure                  ( int *ire,          int *icp,          double *pnmax );
    
    cntc_dll.cntc_getmaximumpressure.restype      = None
    cntc_dll.cntc_getmaximumpressure.argtypes     = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_double) ]
    
    # void cntc_getmaximumtraction                  ( int *ire,          int *icp,          double *ptmax );
    
    cntc_dll.cntc_getmaximumtraction.restype      = None
    cntc_dll.cntc_getmaximumtraction.argtypes     = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_double) ]
    
    # void cntc_getmaximumtemperature               ( int *ire,          int *icp,          double *t1max,
    #                                                 double *t2max );
    
    cntc_dll.cntc_getmaximumtemperature.restype   = None
    cntc_dll.cntc_getmaximumtemperature.argtypes  = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_double),
                                                      POINTER(c_double) ]
    
    # void cntc_getfielddata                        ( int *ire,          int *icp,          int *ifld,
    #                                                 int *lenarr,       double *fld );
    
    cntc_dll.cntc_getfielddata.restype            = None
    cntc_dll.cntc_getfielddata.argtypes           = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_int),    POINTER(c_double) ]
    
    # void cntc_gettractions                        ( int *ire,          int *icp,          int *lenarr,
    #                                                 double *pn,        double *px,        double *py );
    
    cntc_dll.cntc_gettractions.restype            = None
    cntc_dll.cntc_gettractions.argtypes           = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_double), POINTER(c_double), POINTER(c_double) ]
    
    # void cntc_getmicroslip                        ( int *ire,          int *icp,          int *lenarr,
    #                                                 double *sx,        double *sy );
    
    cntc_dll.cntc_getmicroslip.restype            = None
    cntc_dll.cntc_getmicroslip.argtypes           = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_double), POINTER(c_double) ]
    
    # void cntc_getdisplacements                    ( int *ire,          int *icp,          int *lenarr,
    #                                                 double *un,        double *ux,        double *uy );
    
    cntc_dll.cntc_getdisplacements.restype        = None
    cntc_dll.cntc_getdisplacements.argtypes       = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_double), POINTER(c_double), POINTER(c_double) ]
    
    # void cntc_getsensitivities                    ( int *ire,          int *icp,          int *lenout,
    #                                                 int *lenin,        double *sens );
    
    cntc_dll.cntc_getsensitivities.restype        = None
    cntc_dll.cntc_getsensitivities.argtypes       = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_int),    POINTER(c_double) ]
    
    # void cntc_getcalculationtime                  ( int *ire,          int *icp,          double *tcpu,
    #                                                 double *twall );
    
    cntc_dll.cntc_getcalculationtime.restype      = None
    cntc_dll.cntc_getcalculationtime.argtypes     = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_double),
                                                      POINTER(c_double) ]
    
    # void cntc_resetcalculationtime                ( int *ire,          int *icp );
    
    cntc_dll.cntc_resetcalculationtime.restype    = None
    cntc_dll.cntc_resetcalculationtime.argtypes   = [ POINTER(c_int),    POINTER(c_int) ] 
    
    # void subs_getblocksize                        ( int *ire,          int *icp,          int *iblk,
    #                                                 int *nx,           int *ny,           int *nz );
    
    cntc_dll.subs_getblocksize.restype            = None
    cntc_dll.subs_getblocksize.argtypes           = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_int),    POINTER(c_int),    POINTER(c_int) ]
    
    # void subs_getresults                          ( int *ire,          int *icp,          int *iblk,
    #                                                 int *lenarr,       int *ncol,         int *icol,
    #                                                 double *values );
    
    cntc_dll.subs_getresults.restype              = None
    cntc_dll.subs_getresults.argtypes             = [ POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_int),    POINTER(c_int),    POINTER(c_int),
                                                      POINTER(c_double) ]
    
    # void cntc_finalize                            ( int *ire );
    
    cntc_dll.cntc_finalize.restype                = None
    cntc_dll.cntc_finalize.argtypes               = [ POINTER(c_int) ]
    
    # void cntc_finalizelast                        ( );
    
    cntc_dll.cntc_finalizelast.restype            = None
    cntc_dll.cntc_finalizelast.argtypes           = [ ]
  
    return cntc_dll

# end function contact_addon_headers

############################################################################################################## 
