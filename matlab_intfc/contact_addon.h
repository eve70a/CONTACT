/*
 Copyright 2008-2023 by Vtech CMCC.

 Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
*/

void cntc_initializefirst_new(int *ifcver, int *ierror, int *ioutput, const char *c_wrkdir,
                        const char *c_outdir, const char *c_expnam, int *len_wrkdir, 
                        int *len_outdir, int *len_expnam);

void cntc_initializefirst(int *ifcver, int *ierror, int *ioutput, const char *c_outdir, 
                        const char *c_expnam, int *len_outdir, int *len_expnam);

void cntc_initialize(int *ire, int *imodul, int *ifcver, int *ierror, const char *c_outdir,
                        int *len_outdir);

void cntc_setglobalflags(int *lenflg, int *params, int *values);

void cntc_readinpfile(int *ire, int *inp_type, const char *c_fname, int *len_fname, int *ierror);

void cntc_setflags(int *ire, int *icp, int *lenflg, int *params, int *values);

void cntc_setmetadata(int *ire, int *icp, int *lenmta, int *params, double *values);

void cntc_setsolverflags(int *ire, int *icp, int *gdigit, int *nints, int *iparam, int *nreals, 
                        double *rparam);

void cntc_setmaterialparameters(int *ire, int *icp, int *mdigit, int *nparam, double *rparam);

void cntc_settemperaturedata(int *ire, int *icp, int *imeth, int *nparam, double *rparam);

void cntc_settimestep(int *ire, int *icp, double *dt);

void cntc_setreferencevelocity(int *ire, int *icp, double *veloc);

void cntc_setrollingstepsize(int *ire, int *icp, double *chi, double *dq);

void cntc_setfrictionmethod(int *ire, int *icp, int *imeth, int *nparam, double *params);

void cntc_sethertzcontact(int *ire, int *icp, int *ipotcn, int *nparam, double *rparam);

void cntc_setpotcontact(int *ire, int *icp, int *ipotcn, int *nparam, double *rparam);

void cntc_setverticalforce(int *ire, double *fz);

void cntc_setpenetration(int *ire, int *icp, double *pen);

void cntc_setnormalforce(int *ire, int *icp, double *fn);

void cntc_setundeformeddistc(int *ire, int *icp, int *ibase, int *nparam, double *prmudf);

void cntc_setcreepages(int *ire, int *icp, double *vx, double *vy, double *phi);

void cntc_setextrarigidslip(int *ire, int *icp, int *lenarr, double *wx, double *wy);

void cntc_settangentialforces(int *ire, int *icp, double *fx, double *fy);

void cntc_setprofileinputfname(int *ire, const char *c_fname, int *len_fname, int *nints,
                               int *iparam, int *nreals, double *rparam );

void cntc_setprofileinputvalues(int *ire, int *npoint, double *values, int *nints,
                                int *iparam, int *nreals, double *rparam );

void cntc_settrackdimensions(int *ire, int *ztrack, int *nparam, double *rparam);

void cntc_setwheelsetdimensions(int *ire, int *ewheel, int *nparam, double *rparam);

void cntc_setwheelsetposition(int *ire, int *ewheel, int *nparam, double *rparam);

void cntc_setwheelsetvelocity(int *ire, int *ewheel, int *nparam, double *rparam);

void cntc_setwheelsetflexibility(int *ire, int *ewheel, int *nparam, double *rparam);

void subs_addblock(int *ire, int *icp, int *iblk, int *isubs, int *nx, int *ny, int *nz, 
                                          double *xparam, double *yparam, double *zparam);

void cntc_calculate(int *ire, int *icp, int *ierror);

void subs_calculate(int *ire, int *icp, int *ierror);

void cntc_getflags(int *ire, int *icp, int *nparam, int *iparam, int *values);

void cntc_getparameters(int *ire, int *icp, int *itask, int *lenarr, double *values);

void cntc_getprofilevalues(int *ire, int *itask, int *nints, int *iparam, int *lenarr, double *values);

void cntc_getprofilevalues_new(int *ire, int *itask, int *nints, int *iparam, int *nreals, double *rparam,
                                        int *lenarr, double *values);

void cntc_getwheelsetposition(int *ire, int *lenarr, double *values);

void cntc_getwheelsetvelocity(int *ire, int *lenarr, double *values);

void cntc_getnumcontactpatches(int *ire, int *npatch);

void cntc_getcontactlocation(int *ire, int *icp, int *lenarr, double *values);

void cntc_getreferencevelocity(int *ire, int *icp, double *veloc);

void cntc_gethertzcontact(int *ire, int *icp, int *lenarr, double *values);

void cntc_getnumelements(int *ire, int *icp, int *mx, int *my);

void cntc_getgriddiscretization(int *ire, int *icp, double *dx, double *dy);

void cntc_getpotcontact(int *ire, int *icp, int *lenarr, double *values);

void cntc_getpenetration(int *ire, int *icp, double *pen);

void cntc_getcreepages(int *ire, int *icp, double *vx, double *vy, double *phi);

void cntc_getcontactforces(int *ire, int *icp, double *fn, double *tx, double *ty, double *mz);

void cntc_getglobalforces(int *ire, int *icp, int *lenarr, double *values);

void cntc_getcontactpatchareas(int *ire, int *icp, double *carea, double *harea, double *sarea);

void cntc_getelementdivision(int *ire, int *icp, int *lenarr, int *eldiv); 

void cntc_getmaximumpressure(int *ire, int *icp, double *pnmax);

void cntc_getmaximumtraction(int *ire, int *icp, double *ptmax);

void cntc_getmaximumtemperature(int *ire, int *icp, double *t1max, double *t2max);

void cntc_getfielddata(int *ire, int *icp, int *ifld, int *lenarr, double *fld);

void cntc_gettractions(int *ire, int *icp, int *lenarr, double *pn, double *px, double *py);

void cntc_getmicroslip(int *ire, int *icp, int *lenarr, double *sx, double *sy);

void cntc_getdisplacements(int *ire, int *icp, int *lenarr, double *un, double *ux, double *uy);

void cntc_getsensitivities(int *ire, int *icp, int *lenout, int *lenin, double *sens);

void cntc_getcalculationtime(int *ire, int *icp, double *tcpu, double *twall);

void cntc_resetcalculationtime(int *ire, int *icp);

void subs_getblocksize(int *ire, int *icp, int *iblk, int *nx, int *ny, int *nz);

void subs_getresults(int *ire, int *icp, int *iblk, int *lenarr, int *ncol, int *icol, double *values);

void cntc_finalize(int *ire);

void cntc_finalizelast();

