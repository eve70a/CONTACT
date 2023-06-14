!------------------------------------------------------------------------------------------------------------
! m_hierarch_data - declare user-defined types that define a grouping of important variables into
!                   a hierarchical, tree-like data-structure.
!
! Copyright 1993-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_hierarch_data

! using "public", all generic modules used are exported to all using modules

use m_globals
use m_readline
use m_ptrarray
use m_markers
use m_interp_1d
use m_spline_def
use m_grids
use m_gridfunc
use m_interp_gf
use m_bspline_def
use m_bspline_make
use m_bspline_get
use m_spline_make
use m_spline_get
use m_friclaw
use m_inflcf

implicit none

public

   public  t_ic
   public  ic_init
   public  ic_unpack
   public  ic_pack

   public  t_material
   public  mater_init
   public  combin_mater
   public  inflcf_mater
   public  influe_mater

   public  t_potcon
   public  t_hertz
   public  t_geomet

   public  t_kincns
   public  t_solvers
   public  t_leadedge
   public  t_output

   public  t_subsblk
   public  t_subsurf
   public  subsurf_copy
   public  subsblk_copy
   public  subsblk_destroy

   public  t_metadata
   public  meta_init

   public  t_scaling

   public  t_probdata
   public  p_probdata
   public  setini
   public  gd_destroy

   !---------------------------------------------------------------------------------------------------------

   ! codes for the different solvers

   integer, parameter :: isolv_fastsm     = 0, isolv_normcg     = 1,                            &
                         isolv_tangcg     = 2, isolv_cnvxgs     = 3, isolv_stdygs = 4,          &
                         isolv_gdstdy     = 7

   !---------------------------------------------------------------------------------------------------------

   ! codes for computation of sensitivities

   integer, parameter :: nsens_out = 4, nsens_in  = 4
   integer, parameter :: iout_fn   = 1, iout_fx1  = 2, iout_fy1  = 3, iout_mz1  = 4
   integer, parameter :: iin_dpen  = 1, iin_dksi1 = 2, iin_deta1 = 3, iin_dphi1 = 4

   ! nsens_out    number of rows of sensitivities-matrix (1st index), output parameters,
   !              esp. forces and moments
   ! nsens_in     number of columns of sensitivities-matrix, input parameters that can be perturbed

   !---------------------------------------------------------------------------------------------------------

   ! integer control digits for type of problem to be solved:

   type :: t_ic
      integer :: config
      integer :: pvtime
      integer :: bound
      integer :: tang
      integer :: norm
      integer :: force
      integer :: heat
      integer :: stress
      integer :: sens
      integer :: varfrc
      integer :: frclaw_inp
      integer :: discns_inp
      integer :: discns_eff
      integer :: gapwgt
      integer :: gencr_inp
      integer :: mater
      integer :: ztrack
      integer :: ewheel
      integer :: rznorm
      integer :: rztang
      integer :: gausei_inp
      integer :: iestim
      integer :: matfil_surf
      integer :: matfil_subs
      integer :: output_surf
      integer :: output_subs
      integer :: flow
      integer :: nmdbg
      integer :: wrtinp
      integer :: return
      integer :: ilvout
      logical :: print_pmax

      ! c1, config  configuration or composition of wheel/rail problems
      !              0 = wheelset on track, left wheel,
      !              1 = wheelset on track, right wheel,
      !              4 = wheelset on rollers, left wheel,
      !              5 = wheelset on rollers, right wheel
      ! p, pvtime   the filling of the tractions Pv for the previous time instance,
      !              0 = full sequence; copy result Ps to Pv;
      !              1 = continuation for the normal part only;
      !              2 = no continuation, the entire Pv is cleared; T=1,2: initiation of contact;
      !             note: p has no effect for t=0 and 3.
      ! b, bound    selects the approach to be used for the normal problem, the traction bound:
      !              0 = full linearly elastic model and contact conditions;
      !              1 = (not used, reserved for semi-Hertzian approach);
      !              2 = elliptical traction bound derived from Hertzian problem data (requires IPOTCN -3:-1);
      !              3 = parabolical traction bound using Hertzian problem data (requires IPOTCN -3:-1).
      !              4 = elliptical traction bound derived from SDEC problem data (requires IPOTCN=-6).
      !              5 = non-Hertzian approximation of the normal problem using the KPEC approach.
      !              6 = non-Hertzian approximation of the normal problem using ANALYN method.
      ! t, tang     specifies the type of tangential problem to be solved
      !              0 = none,
      !              1 = shift,             using a world-fixed coordinate system,
      !              2 = transient rolling, using a contact-fixed moving coord.sys
      !              3 = steady rolling, direct method,
      !              5,6,7 = reserved
      !             note: options t=1 and 2 depend on the p-digit.
      !               initial shift:         t=1, p=2, continuation: t=1, p=0.
      !               initiation of rolling: t=2, p=2, continuation: t=2, p=0.
      ! n, norm     specifies the normal problem:
      !              0 = penetration/approach prescribed,
      !              1 = normal force prescribed
      ! f, force    specifies the tangential problem:
      !              0 = creepages cksi, ceta prescribed;
      !              1 = force fxrel1, creepage ceta prescribed;
      !              2 = forces fxrel1, fyrel1 prescribed
      ! h, heat     activates the temperature calculation:
      !              0 = no surface temperature calculation;
      !              1 = surface temperature calculation using parameters stored in memory
      !              3 = read new input data and perform temperature calculation for steady rolling
      ! s, stress   determines the operation of the subprogram STRESS
      !              0 = no subsurface stresses required for this case;
      !              1 = compute the stresses in the points already stored in XB;
      !              2 = read new flags, maintain points already stored in XB;
      !              3 = read new subsurface points and compute subsurface stresses for these points.
      !    sens     specifies whether sensitivities must be computed:
      !              0 = no sensitivities needed for this case
      !              1 = print sensitivities that are readily available
      !              2 = compute normal sensitivity
      !              3 = compute normal and tangential sensitivities
      ! v1, varfrc  concerns variation of friction across width/rail profile (module 1)
      !              0 = constant inputs, no variation
      !              1 = multiple inputs with linear interpolation
      !              2 = multiple sets of inputs are used, per row of the potential contact
      !                  used internally in module 3 for conformal contacts
      ! l, frclaw   concerns the friction law to be used:
      !              0 = Coulomb friction with static and kinetic friction coefficients;
      !              1 = use the friction law from parameters in storage;
      !              2 = friction law with linear/constant dependency on slip velocity
      !              3 = friction law with rational dependency on slip velocity
      !              4 = friction law with exponential dependency on slip velocity
      !              5 = used for Polach parameters in CONTACT library,
      !                  internally converted to exponential friction with L=4.
      !              6 = friction law with piecewise linear dependency on surface temperature
      !              7 = reserved
      !             Note: the variable ic%frclaw_inp contains the value obtained from (and written to) the
      !                   input-file, whereas fric%frclaw_eff is the effective law used in the program (<>1).
      ! d, discns   concerns the potential contact area and discretisation:
      !              0 = maintain the discretisation of the previous case;
      !              1 = form the discretisation from parameters in storage;
      !              2 = "wgt-wgt-2" planar contact with contact locus, weighed center and weighted angle,
      !                  read new input parameters and form the discretisation
      !              3 = "wx-surf",    same as 2, with extended rigid slip computation
      !              4 = "curved",     same as 2, with curved reference and conformal contact computation
      !              5 = "brute",      same as 2, using brute force instead of contact locus
      !              6 = "brute+steep" same as 5, using larger grid to accomodate steep slopes
      !              8 = "locus+view"  same as 2, using a hard-coded oblique view direction
      !              9 = "brute+view"  same as 5, using a hard-coded oblique view direction
      !             Note: the variable ic%discns_inp contains the value obtained from (and written to) the
      !                   input-file, whereas ic%discns_eff is the effective method used (2 -- 5),
      !                   with inputs 6,8,9 rewritten to 2,5.
      !     gapwgt  variations of weighted center method: 0 == default, 1 == linear weighting, 2 == quadratic
      ! c3, gencr   concerns the material parameters and influence coefficients
      !              0 = maintain inluence coefficients;
      !              1 = generate influence coefficients from parameters in storage;
      !              2 = analytical (half-space) influence coefficients, piecewise constant; read parameters
      !                  and compute.
      !              3 = analytical (half-space) influence coefficients, bilinear; read parameters and compute
      !              4 = reserved for Blanco's half-space + angle correction
      !              5 = reserved (layered half-space)
      !              9 = numerical influence coefficients; read parameters, read infl.coefficients from file
      !             Note: the variable ic%gencr_inp contains the value obtained from (and written to) the
      !                   input-file, whereas mater%gencr_eff is the effective type used (2, 3 or 9).
      ! m, mater    type of material model to be used:
      !              0 = (linearly) elastic material;
      !              1 = viscoelastic material (standard linear solid);
      !              2 = simplified theory with flexibility L;
      !              3 = simplified theory with flexibilities L1,L2,L3.
      !              4 = linearly elastic material + elasto-plastic third body layer
      !              5 = reserved (pseudo-viscous damping)
      !              6 = reserved (vertical slice/gap)
      !              Note: the m-digit is copied to t_material
      ! z1, ztrack  concerns the track geometry, profile and deviation (module 1)
      !              0 = maintain track dimensions, profile and deviations;
      !              1 = read new dimensions and profile;
      !              2 = read new track deviations;
      !              3 = read new dimensions, one profile for current side of track, and track deviations,
      ! e1, ewheel   concerns the wheel-set geometry, profile, position and velocity (module 1)
      !              0 = maintain wheel-set geometry, profile, position and velocity;
      !              1 = read new position data;
      !              2 = read new position and velocity data;
      !              3 = read new position and velocity data, geometry, profile for current side of wheelset
      !              4 = read new position and velocity data, including flexible wheelset deviations
      !                  for current side of wheelset
      !              5 = as 3, including flexible wheelset deviations for current side of wheelset
      ! z3, rznorm   concerns the right hand side of the normal problem (module 3)
      !              0 = maintain undeformed distance and planform;
      !              1 = form undeformed distance from parameters in storage;
      !              2 = read new parameters and compute undeformed distance.
      ! e3, rztang/exrhs  extra term of the rigid slip, rhs of tangential problem (module 3)
      !              0 = set the extra term equal to zero;
      !              1 = maintain the extra term of the previous case;
      !              9 = read a new extra term from the inputfile, element-by-element
      ! g, gausei   the choices w.r.t. iterative solvers
      !              0 = default solvers and settings, new max.its, eps
      !              1 = keep settings from previous case
      !              2 = use ConvexGS at all times, read omegah/s from input-file
      !              3 = use SteadyGS when possible (T=3), read omegah/s from input
      !              4 = use default solver, read inislp, omgslp from input
      !              5 = use GDsteady when possible (T=3, M=0), read fdecay, d_ifc etc from input
      !             Note: the variable ic%gausei_inp contains the value obtained from (and written to) the
      !                   input-file, whereas solv%gausei_eff is the effective value used in the program (<>1)
      ! i, iestim   governs the initial estimate:
      !              0 = no good initial estimate available;
      !              1 = use previous solution, regularize tractions;
      !              2 = use normal part of solution only;
      !              3 = previous case gives a good initial estimate.
      !              5 = keep initial estimate as set, obtained from coarse grid.
      ! a, matfil_surf   governs the use of the matlab-file <experim>.<case>.mat per case:
      !              0 = the mat-files is not created;
      !              1 = the detailed results of the case are written to a mat-file, for points inside
      !                  contact area.
      !              2 = results are evaluated/written to the mat-file for all points of the potential
      !                  contact area.
      ! as, matfil_subs  governs the writing of subsurface output to the subs-file:
      !              0 = the subs-file is not created;
      !              1 = the displacements and stress-invariants are written to the subs-file;
      !              2 = additionally, the full stress tensor is written to the subs-file.
      ! o, output_surf   governs the extent of the output to the output-file:
      !              0 = no results are printed;
      !              1 = minimum output, just the computed results;
      !              2 = overview of global quantities, inputs and outputs;
      !              3 = print picture of contact, adhesion and slip areas;
      !              4 = print (principal) profiles used in track coordinates;
      !              5 = detailed solution too, inside actual contact area;
      !              6 = detailed solution, full potential contact area.
      ! os, output_subs  governs the extent of the output on subsurf.stress to the out-file:
      !              0 = no results are printed;
      !              1 = minimum output, just the maximum values of primary stress invariants;
      !              2 = print additional maximum values: Tresca, principal stresses;
      !              3 = not used
      !              4 = print detailed results of subsurface stress (displacements, invariants).
      ! w, flow     governs the extent of the flow trace.
      !              0   = none,
      !              1-4 = flow of iteration processes Panag, Duvorol, Norm/Tang, Newton-Raphson,
      !              5,6 = flow of iterative solvers ConvexGS, ConjGrad,
      !              9   = full flow-trace.
      ! x, nmdbg    governs the extent of (debug) grid data printed to the output-file
      !              0-3 = none
      !              4   = global problem inputs Hs and outputs Igs,Ps,Us
      !              5   = intermediate results of the Panagiotopoulos' process
      !              6   = intermediate results of Duvorol and Veloc.dep iterations
      !              7   = intermediate results of NORM and TANG algorithms
      !              8   = intermediate results of Newton-Raphson loops
      !              9   = information per iteration of NormCG, ConvexGS, SteadyGS
      !    wrtinp   governs the writing of data to the input-file, used by the CONTACT add-on to SIMPACK
      !             Rail, allowing for off-line use of CONTACT after the SIMPACK run.
      !              0   = no writing of the input (default)
      !              1   = write input of each case to the .inp-file.
      ! r, return   return to main program.
      !              0 = calculate, stay in module, 1 = calculate and return,
      !              2 = skip calculation, stay, 3 = skip calculation and return
      ! print_pmax  .true.: compute/show pmax with aggregate results per patch
      !    ilvout   level of output to screen and .out-file, mainly for testing
      !              0 = minimum output, 1 = normal output

   end type t_ic

   !---------------------------------------------------------------------------------------------------------

   ! data with respect to the (elastic, viscoelastic) materials of the two bodies:

   type :: t_material
      real(kind=8) :: gg(2)
      real(kind=8) :: poiss(2)
      real(kind=8) :: ak
      real(kind=8) :: ga
      real(kind=8) :: nu
      real(kind=8) :: bktemp(2)
      real(kind=8) :: heatcp(2)
      real(kind=8) :: lambda(2)
      real(kind=8) :: dens(2)
      real(kind=8) :: betapl
      real(kind=8) :: fg(2)
      real(kind=8) :: tc(2)
      real(kind=8) :: vt(2)
      real(kind=8) :: akv
      real(kind=8) :: gav
      real(kind=8) :: nuv
      real(kind=8) :: flx_z
      real(kind=8) :: flx(3)
      real(kind=8) :: k0_mf
      real(kind=8) :: alfamf
      real(kind=8) :: betamf
      real(kind=8) :: k_eff
      real(kind=8) :: gg3
      real(kind=8) :: laythk
      real(kind=8) :: tau_c0
      real(kind=8) :: k_tau
      integer      :: bound_eff
      integer      :: mater_eff
      integer      :: gencr_eff
      integer      :: if_meth, if_ver, ninclin
      real(kind=8), dimension(:,:), allocatable :: surf_inclin
      character(len=256) :: fname_influe

      ! gg     [N/mm2] modulus of rigidity or shear modulus for bodies (1, 2) (steel: ~76 GPa = 76.000 N/mm2)
      ! poiss   [-]    Poisson ratio for bodies (1, 2) (steel: ~0.3, rubber: near 0.5, nearly incompressible)
      ! ak      [-]    difference parameter K, see Kalker (1.44)
      ! ga     [N/mm2] combined modulus of rigidity G, see Kalker (1.44)
      ! nu      [-]    combined Poisson's ratio "nu", see Kalker (1.44)

      ! bktemp  [C]     background / bulk / initial temperature for bodies (1, 2)
      ! heatcp [J/kg-C] specific heat capacity for bodies (1, 2)
      ! lambda [W/mm-C] thermal conductivity for bodies (1, 2)
      ! dens   [kg/mm3] density for bodies (1, 2)
      ! betapl  [-]     fraction of plastic work dissipated as heat

      ! fg      [-]    parameter for viscoelastic bodies (1, 2): Ratio of Cg to Cv spring constants, see
      !                user guide.
      ! tc      [s]    parameter for viscoelastic bodies (1, 2): creep relaxation time
      ! vt      [mm]   parameter for viscoelastic bodies (1, 2): relaxation distance t_c*V.
      ! akv     [-]    difference parameter K for viscoelastic materials
      ! gav    [N/mm2] combined modulus of rigidity G for viscoelastic materials
      ! nuv     [-]    combined Poisson's ratio "nu" for viscoelastic materials

      ! flx_z  [mm3/N] flexibility L_z used for the vertical compression of a thin sheet (B=1)
      ! flx    [mm3/N] flexibilities L1, L2 and L3 used in the simplified theory
      !                M = 2: L1 = L2 = L prescribed by user; M = 3: computed from linear theory
      !                M = 4: using L1 for flexibility of 3rd body layer
      ! k0_mf   [-]    slope reduction factor k_0 of modified Fastsim algorithm
      ! alfamf  [-]    slope reduction multiplier alpha_inf = k_inf/k_0 of modified Fastsim algorithm
      ! betamf  [-]    creep multiplier beta of modified Fastsim algorithm
      ! k_eff   [-]    effective slope reduction factor k of modified Fastsim algorithm (computed)

      ! gg3    [N/mm2] shear elastic modulus for body 3 (interface layer)
      ! laythk  [mm]   thickness of the interface layer
      ! tau_c0 [N/mm2] initial shear limit at which plasticity effects start to occur
      ! k_tau  [N/mm3] rate of increase of the shear limit with accumulated plastic deformation

      ! bound_eff      effective bound-digit used, see field "bound" in type t_ic.
      ! mater_eff      effective material model used, see field "mater" in type t_ic.
      ! gencr_eff      effective type of influence coefficients, see "gencr_inp" in type t_ic.
      ! if_meth        method number for Blanco-correction: 0 = fast (linear angle), 1 = detailed (full angle)
      ! if_ver         version number for Blanco-correction: v1 (orig), 2 (no xn), 3 (un) or 4 (split)
      ! ninclin        number of points in table of surface inclinations
      ! surf_inclin    table (y_i, alpha_i) of surface inclinations for Blanco-approach
      ! fname_influe   file-name for file with numerical influence coefficients as given by user

   end type t_material

   !---------------------------------------------------------------------------------------------------------

   ! data with respect to the discretisation, potential contact area:

   type :: t_potcon
      integer      :: ipotcn

      ! ipotcn  specification method for potential contact area:
      !          -1 == Hertzian ellipse by curvatures (a1,b1)
      !          -2 == Hertzian ellipse by curvature and ellipticity (a1,aob)
      !          -3 == Hertzian ellipse by semi-axes (aa,bb)
      !          -4 == Hertzian rectangular by curvature a1 and half-width bb
      !          -5 == Hertzian rectangular by half-length aa and half-width bb
      !          -6 == two half-ellipses by semi-axes (aa,bpos,bneg)
      !           1 == by lower-left corner (xl,yl) and step-sizes (dx,dy)
      !           2 == by lower-left (xl,yl) and upper-right corners (xh,yh)
      !           3 == by center of first element (xc1,yc1) and step-sizes (dx,dy)
      !           4 == by centers of first (xc1,yc1) and last elements (xcm,ycm)

      integer      :: mx, my, npot
      real(kind=8) :: dx, dy, dxdy

      ! mx             number of elements of potential contact in x-direction
      ! my             number of elements of potential contact in y-direction
      ! npot           total number of elements of potential contact area mx*my
      ! dx     [mm]    size of rectangular elements in x-direction
      ! dy     [mm]    size of rectangular elements in y-direction
      ! dxdy   [mm2]   surface dx*dy of each element

      real(kind=8) :: xl
      real(kind=8) :: yl
      real(kind=8) :: xh
      real(kind=8) :: yh
      real(kind=8) :: xc1
      real(kind=8) :: xcm
      real(kind=8) :: yc1
      real(kind=8) :: ycm

      ! xl     [mm]    x-coordinate of lower-left corner of potential contact area
      ! yl     [mm]    y-coordinate of lower-left corner of potential contact area
      ! xh     [mm]    x-coordinate of upper-right corner of pot. contact area
      ! yh     [mm]    y-coordinate of upper-right corner of pot. contact area
      ! xc1    [mm]    x-coordinate of center of first element (1,1) of pot.con.
      ! yc1    [mm]    y-coordinate of center of first element (1,1) of pot.con.
      ! xcm    [mm]    x-coordinate of center of last element (mx,my) of pot.con.
      ! ycm    [mm]    y-coordinate of center of last element (mx,my) of pot.con.

   end type t_potcon

   !---------------------------------------------------------------------------------------------------------

   ! data with respect to the Hertzian problem-description:

   type :: t_hertz
      real(kind=8) :: a1
      real(kind=8) :: aa
      real(kind=8) :: aob
      real(kind=8) :: b1
      real(kind=8) :: bb
      real(kind=8) :: bneg
      real(kind=8) :: bpos
      real(kind=8) :: scale
      real(kind=8) :: cp
      real(kind=8) :: rho

      ! a1, b1 [1/mm]  "curvature" in x-(rolling) and y-(transverse)-directions
      ! aa, bb [mm]    basic semi-axes of contact ellipse (rectangle) in x-/y-directions
      ! bneg,pos [mm]  semi-axes of two half ellipses in y-direction (SDEC method)
      ! aob    [-]     ratio of semi-axes aa/bb
      ! scale  [-]     multiplication factor for size of potential contact area in Hertzian problem

      ! cp     [mm]    characteristic length scale of contact patch, sqrt(aa*bb). See Kalker (1.58d)
      ! rho    [mm]    characteristic radius of curvature 2/(a1+b1). Equals 2/(d_1+d_2)?

      ! NOTE: cp and rho have not been introduced everywhere, e.g. are not set for non-Hertzian
      !       cases and not when F=1,2.

   end type t_hertz

   !---------------------------------------------------------------------------------------------------------

   ! data for one block of points for subsurface stress calculation

   type :: t_subsblk
      integer      :: isubs
      integer      :: ncolum
      integer      :: nx_inp, ixl_inp, ixinc, ixh_inp, nx_eff, ixl_eff, ixh_eff
      integer      :: ny_inp, iyl_inp, iyinc, iyh_inp, ny_eff, iyl_eff, iyh_eff
      integer      :: nz
      real(kind=8) :: zl, dz
      integer,      dimension(:),   pointer   :: ixlist_inp => NULL(), iylist_inp => NULL()
      integer,      dimension(:),   pointer   :: ixlist_eff => NULL(), iylist_eff => NULL()
      real(kind=8), dimension(:),   pointer   :: x      => NULL(), y      => NULL(), z      => NULL()
      real(kind=8), dimension(:,:), pointer   :: table  => NULL()

      ! isubs   -      type of input specification
      !                 1 - centers of all elements of potential contact, regular spacing dz
      !                 2 - centers of regular selection of elements,     regular spacing dz
      !                 3 - centers of irregular selection of elements,   regular spacing dz
      !                 5 - centers of all elements of potential contact, explicit selection of z
      !                 6 - centers of regular selection of elements,     explicit selection of z
      !                 7 - centers of irregular selection of elements,   explicit selection of z
      !                 9 - explicit specification of x, y and z
      ! ncolum    -    number of columns provided in table
      ! nx,ny,nz  -    number of points in resp. x, y and z-directions
      ! ixl, ixh  -    isubs=2,6: low and high element numbers in x-direction
      ! ixinc     -    isubs=2,6: increment in x-direction
      ! iyl, iyh  -    isubs=2,6: low and high element numbers in y-direction
      ! iyinc     -    isubs=2,6: increment in y-direction
      ! ixlist    -    isubs=3,7: explicit list of element numbers in x-direction
      ! iylist    -    isubs=3,7: explicit list of element numbers in y-direction
      ! zl, dz  [mm]   isubs=1-3: lowest z-coordinate, increment in z-coordinate
      ! x, y, z [mm]   coordinates in resp. x, y and z-directions
      ! table   [mm, N/mm^2]   computed values [xyz], u[xyz], sig{hyd,vm,tr}, sigma1-3, sigma

   end type t_subsblk

   !---------------------------------------------------------------------------------------------------------

   ! aggregate data for subsurface stress calculation

   type :: t_subsurf
      integer         :: nblock
      type(t_subsblk) :: blocks(MXBLCK)

      ! nblock         number of blocks of subsurface points
      ! blocks         input data for the blocks of subsurface points
   end type t_subsurf

   !---------------------------------------------------------------------------------------------------------

   ! data with respect to the geometry-description of the bodies in the contact zone:

   type :: t_geomet
      integer                                      :: ibase, iplan
      integer                                      :: nn, npatch
      real(kind=8),    dimension(:),   pointer     :: prmudf => NULL()
      real(kind=8),    dimension(:),   pointer     :: prmpln => NULL()
      real(kind=8),    dimension(:),   pointer     :: ysep   => NULL()
      real(kind=8),    dimension(:,:), pointer     :: facsep => NULL()
      type(t_gridfnc3)                             :: exrhs
      type(t_gridfnc3)                             :: hs1, hv1

      ! ibase          parameterisation method used for undeformed distance
      !                  1 == quadratic, 6 coeffs given in prmudf
      !                  2 == circular in x, pointwise in y
      !                  3 == quadratic + difference of two sines in x-direction
      !                  9 == h specified at centers of all elements
      ! iplan          parameterisation method used for planform
      !                  1 == unrestricted planform
      !                  2 == quadratic planform, 6 coeffs in prmpln
      !                  3 == union of two rectangles, 8 coeffs in prmpln
      !                  4 == weighted interaction between different (sub-) patches
      ! nn             number of heights for ibase==2
      ! npatch         number of (sub-) patches distinguished in iplan==4.
      !                data stored in ysep(1:np), facsep(1:np,1:np).
      ! prmudf [?]     parameters for computing the undeformed distance
      ! prmpln [?]     coefficients of the planform
      ! ysep   [mm]    y-positions in the middle between successive (sub-) patches
      ! facsep [-]     multiplication factors for effect of pressures on displacements between (sub-)patches
      ! exrhs [mm|-]   extra term of the rigid slip (T=1, [mm]) or creepage (T=2-3, [-])
      ! hs1    [mm]    1st column (vn): undeformed distance
      !                columns 2,3: rigid slip distance of body (1) w.r.t. (2) of the current time-step
      ! hv1    [mm]    undeformed distance/rigid slip of previous time, cf. hs1

   end type t_geomet

   !---------------------------------------------------------------------------------------------------------

   ! data with respect to the kinematic constants, relative movement of the two bodies with respect to
   ! each other (both input and output)

   type :: t_kincns
      real(kind=8) :: veloc
      real(kind=8) :: dq
      real(kind=8) :: dt
      real(kind=8) :: facphi
      real(kind=8) :: chi
      real(kind=8) :: pen
      real(kind=8) :: penv
      real(kind=8) :: cksi
      real(kind=8) :: ceta
      real(kind=8) :: cphi
      real(kind=8) :: fntrue
      real(kind=8) :: fnscal
      real(kind=8) :: fxrel1
      real(kind=8) :: fyrel1
      real(kind=8) :: fprev(3)
      real(kind=8) :: muscal
      logical      :: use_muscal

      ! veloc  [mm/s]  rolling velocity
      ! chi    [rad]   rolling direction
      ! dq     [mm]    distance traversed per time step in rolling problems
      ! dt     [s]     time step size
      ! facphi [-]     multiplication factor "1/6" for mysterious additional term for spin creepage

      ! pen    [mm]    approach/penetration of the two bodies in normal direction
      !                (prescribed when N=0, computed when N=1)
      ! penv   [mm]    approach/penetration of previous time instance
      ! cksi  [mm|-]   T=1: rigid shift in x-direction ([mm]),
      !                T=2-3: creepage in x-direction ([-]).
      !                (prescribed when F=0, computed when F>=1)
      ! ceta  [mm|-]   rigid shift (T=1) or creepage (T=2-3) in y-direction
      !                (prescribed when F<=1, computed when F=2)
      ! cphi [rad|rad/mm]  rigid spin shift (T=1, [rad]) or spin creepage (T=2-3, [rad/mm])
      ! TODO: it would be better to have separate variables for creep & shift

      ! fntrue [N]     total normal force, input when N=1, output when N=0
      ! fnscal [mm2]   total normal force relative to combined modulus of rigidity (output)
      ! fxrel1 [-]     total tangential force in x-direction _on_ body (1), relative to muscal*fntrue;
      !                input when F>=1, output when F=0.
      ! fyrel1 [-]     total tangential force in y-direction _on_ body (1), relative to muscal*fntrue;
      !                input when F=2, output when F<=1.
      ! fprev  [N,-]   total forces [fx1,fy1,fn1] of previous time instance

      ! use_muscal     flag indicating whether tangential forces are scaled by MU*FN (true) or by FN (false)
      ! muscal [-]     coefficient of friction used in scaling of tangential forces. Typically equal to FSTAT.

   end type t_kincns

   !---------------------------------------------------------------------------------------------------------

   ! variables related to the solution algorithms:

   type :: t_solvers
      integer      :: gausei_eff
      integer      :: solver_eff
      integer      :: gd_meth
      integer      :: maxout
      integer      :: maxin
      integer      :: maxnr
      integer      :: maxgs
      integer      :: mxsens
      integer      :: inislp
      real(kind=8) :: eps
      real(kind=8) :: epsens
      real(kind=8) :: omegah
      real(kind=8) :: omegas
      real(kind=8) :: omgslp
      real(kind=8) :: fdecay
      integer      :: kdown
      integer      :: kdowfb
      real(kind=8) :: betath
      real(kind=8) :: d_ifc
      real(kind=8) :: d_lin
      real(kind=8) :: d_cns
      real(kind=8) :: d_slp
      real(kind=8) :: pow_s
      integer      :: itnorm
      integer      :: ittang
      integer      :: itcg
      integer      :: itgs

      ! gausei_eff     effective solver configuration used in current case, see ic%gausei_inp above.
      ! solver_eff     actual solver to be used, see isolv_fastsm -- isolv_stdygs above
      ! gd_meth   primary variant of GDsteady: 1=E_trl, 2=E_down(k), 3=E_keep(f)
      ! maxout    maximum number of iterations for outer Panag. process
      ! maxin     maximum number of iterations for inner processes, algorithms Norm and Tang
      ! maxnr     maximum number of iterations for Newton and Newton-Raphson processes
      ! maxgs     maximum number of iterations for iterative solvers CG and CnvxGS
      ! mxsens    maximum number of iterations for CG and GS when computing sensitivities
      ! inislp    initial estimate for absolute slip velocity S in velocity dependent friction laws (L=2-4).
      !           when >0, start with S=0 and use increasing S,
      !           when =0, start with final S of previous case,
      !           when <0, start with S too large and go downwards
      ! eps       relative tolerance (epsilon) for solution processes
      ! epsens    relative tolerance for solution processes when computing sensitivities
      ! omegah    relaxation factor for CnvxGS and StdyGS for elements in the adhesion area
      ! omegas    relaxation factor for elements in the slip area
      ! omgslp    relaxation factor for velocity dependent friction iteration
      ! fdecay    reduction factor per grid point for sum(dp) in GDsteady variant 3, E_keep(f)
      ! kdown     number of grid points for negative dp in GDsteady variant 2, E_down(k)
      ! kdowfb    number of grid points for negative dp in E_down(k) when used as fallback method
      ! betath    threshold for coefficient beta in GDsteady for signalling stagnation
      ! d_ifc     diagonal scaling coefficient for GDsteady, value in adhesion at interface
      ! d_lin     diagonal scaling coefficient for GDsteady, change per grid point in adhesion area
      ! d_cns     diagonal scaling coefficient for GDsteady, value in adhesion far from interface
      ! d_slp     diagonal scaling coefficient for GDsteady, multiplication factor in slip area
      ! pow_s     diagonal scaling coefficient for GDsteady, exponent used in slip area

      ! itnorm    number of iterations in the NORM algorithm/subroutine, <0 means error
      ! ittang    number of iterations in the TANG algorithm/subroutine, <0 means error
      ! itcg      number of iterations used in NormCG, solving the pressures
      ! itgs      number of iterations used for solving the tangential tractions, by ConvexGS, SteadyGS or
      !           TangCG

   end type t_solvers

   !---------------------------------------------------------------------------------------------------------

   ! data with respect to the position of the leading edge of the contact area:

   type :: t_leadedge

      integer                               :: ixinc
      integer                               :: npos
      integer,      dimension(:),   pointer :: jbnd  => NULL()
      integer,      dimension(:),   pointer :: ixbnd => NULL()
      type(t_igrdfnc1)                      :: ii2j
      real(kind=8), dimension(:),   pointer :: xbnd  => NULL()
      real(kind=8), dimension(:),   pointer :: facdx => NULL()
      real(kind=8), dimension(:,:), pointer :: ubnd  => NULL()
      type(t_gridfnc3)                      :: facdt

      ! Note: we assume rolling in either positive (chi=0) or negative x-direction (chi=pi). In
      !       the former case, a leading edge position has x_interior <= x_bnd <= x_exterior, in the
      !       latter case int&ext are reversed. There may be 0, 1 or a few leading edge positions per
      !       grid row iy.

      ! ixinc   increment +1 or -1 for resp. rolling in positive (0) or negative x-direction (pi)
      ! npos    number of leading edge positions in whole grid, size of arrays ixbnd, xbnd, facdx, etc.
      ! jbnd    for each grid row iy, the first corresponding position j in arrays ixbnd, xbnd.
      !         Note: the number of positions occupied for row iy is jbnd(iy+1)-jbnd(iy).
      !         The array-size is my+1. Total number of positions is npos = jbnd(my+1)
      ! ixbnd   for each leading edge position, the element coordinate ix for the interior element
      !         adjacent to an exterior element at ix+1 (0) or ix-1 (pi). Array-size is given by npos.
      ! ii2j    for each element ii of the potential contact area, the corresponding index j in the leading
      !         edge arrays or 0 if ii is not a leading edge boundary.
      ! xbnd    for each leading edge position, the estimated x-coordinate. Array-size is given by npos.
      ! facdx   for each leading edge position, the relative position in terms of the step dx between ix
      !         and ix+1 (0) or ix-1 (pi). Array-size is given by npos.
      !         Note: facdx = 0 for center of ix, facdx = 1 at center of ix+1 or ix-1.
      ! facdt   for each element of the potential contact area, the fraction of the time-step that it is
      !         inside the actual contact area.
      !         facdt=0 for elements in the Exterior.
      !         facdt=1 for elements that were already in the contact area at the previous time instance.
      !         0 <= facdt <= 1 for elements entering the contact area in the current step.
      ! ubnd    for each leading edge position, the estimated displacement difference in tangential
      !         directions (2,3). Array-size is given by (npos,3).

   end type t_leadedge

   !---------------------------------------------------------------------------------------------------------

   ! data with respect to the solution and derived quantities:

   type :: t_output
      type(t_eldiv)    :: igs, igv
      type(t_leadedge) :: ledg
      type(t_gridfnc3) :: mus, muv, shft
      type(t_gridfnc3) :: ps, pv, us, uv, ss, sv
      type(t_gridfnc3) :: taucs, taucv, upls, uplv
      type(t_gridfnc3) :: temp1, temp2

      real(kind=8)     :: sens(nsens_out, nsens_in)
      real(kind=8)     :: mxtrue, mytrue, mztrue
      real(kind=8)     :: elen, frpow, pmax

      ! igs            element division for current solution, see t_eldiv.
      ! igv            element division for previous solution.
      ! ledg           administration w.r.t. position of leading edge of the contact area

      ! ps     [N/mm2] current surface tractions in normal and tangential directions
      ! mus    [-]     friction coefficient per element for current time step
      ! muv    [-]     friction coefficient per element for previous time step
      ! pv     [N/mm2] previous surface tractions, see ps
      ! us     [mm]    current surface displacement differences u = A*ps
      ! uv     [mm]    previous surface displacement differences u' = A'*pv
      ! taucs  [N/mm2] current shear limit tau_c for tangential plastic deformation
      ! taucv  [N/mm2] previous shear limit tau'_c for tangential plastic deformation
      ! upls   [mm]    current plastic deformation differences u_pl
      ! uplv   [mm]    previous plastic deformation differences u'_pl
      ! ss     [mm]    slip distance (shift) S_It at the current time-instance, integrated (micro-)slip
      !                velocity == s_abs * dt == s_rel * dq
      ! sv     [mm]    the slip distance S_It at the previous time-instance
      ! shft   [mm]    magnitude of the shift ss per element
      ! temp1  [*C]    current surface temperature of body 1
      ! temp2  [*C]    current surface temperature of body 2
      ! sens           sensitivities of output forces and moments
      !                  fntrue [N], fxrel1, fyrel1 [-], mztru1 [N.mm] etc.
      !                w.r.t. input parameters:
      !                  pen [mm], cksi1, ceta1 [mm], cphi1 [rad] etc. (T=1)
      !                  pen [mm], cksi1, ceta1 [-], cphi1 [rad/mm] etc. (T=2,3)
      !                Note: one array for whole contact-problem comprising kin1 & kin2 --> outpt1 & outpt2
      ! mxtrue [N.mm]  torsional moment around the x-axis
      ! mytrue [N.mm]  torsional moment around the y-axis
      ! mztrue [N.mm]  torsional moment around the z-axis
      ! elen   [N.mm]  elastic energy
      ! frpow  [N.m/s] frictional power dissipation
      ! pmax   [N.mm2] maximum pressure

   end type t_output

   !---------------------------------------------------------------------------------------------------------

   ! meta-data providing additional information about the calculation:

   type :: t_metadata
      integer            :: REid
      integer            :: CPid
      integer            :: actv_thrd
      integer            :: irun, iax, iside, ncase, itbrent
      integer            :: npatch, ipatch
      real(kind=8)       :: tim, s_ws, ynom_whl, rnom_whl, rnom_rol
      real(kind=8)       :: x_rw, y_rw, z_rw, rollrw, yawrw
      real(kind=8)       :: y_rr, z_rr, rollrr
      real(kind=8)       :: xcp_tr, ycp_tr, zcp_tr, deltcp_tr
      real(kind=8)       :: xcp_rr, ycp_rr, zcp_rr, scp_rr, deltcp_rr
      real(kind=8)       :: xcp_rw, ycp_rw, zcp_rw, scp_rw, deltcp_rw
      character(len=256) :: dirnam, expnam

      ! REid              result element ID, used by CONTACT add-on
      ! CPid              contact patch ID, used by CONTACT add-on
      ! actv_thrd         OpenMP thread number working on this gd, -1=none, used for locking
      ! irun              for Sentient: run number
      ! iax               for Sentient: axle number
      ! iside             for Sentient: side number
      ! ncase             case number for the (REid,CPid)
      ! itbrent           iteration number for vertical force iteration
      ! npatch            number of contact patches, used by module 1
      ! ipatch            contact patch number, used by module 1
      ! tim        [s]    (SIMPACK) simulation time
      ! s_ws       [mm]   location of wheel-set CM along the track curve
      ! ynom_whl   [mm]   lateral position of wheel origin in wheelset coordinates
      ! rnom_whl   [mm]   nominal wheel radius == vertical position of wheel origin in wheelset coords
      ! rnom_rol   [mm]   nominal roller radius
      ! [xyz]_rw   [mm]   location of right wheel origin m_rw in terms of track coordinates
      ! roll_rw    [rad]  roll angle of right wheel (marker m_rw) with respect to track coordinates
      ! yaw_rw     [rad]  yaw angle of right wheel (marker m_rw) with respect to track coordinates
      ! [yz]_rr    [mm]   location of right rail origin m_rr in terms of track coordinates
      ! roll_rr    [rad]  roll angle of right rail marker m_rr with respect to track coordinates
      ! [xyz]cp_tr [mm]   location of contact reference point m_ref (m_cp) in terms of track coordinates
      ! deltcp_tr  [rad]  contact angle: roll angle from track vertical to contact reference normal direction
      ! [xyz]cp_rr [mm]   location of contact reference point m_ref (m_cp) in terms of right rail coordinates
      ! scp_rr     [mm]   position of the contact reference point measured along the curved rail surface
      ! deltcp_rr  [rad]  roll angle from right rail vertical to contact reference normal direction
      ! [xyz]cp_rw [mm]   location of contact reference point m_ref (m_cp) in terms of right wheel coordinates
      ! scp_rw     [mm]   position of the contact reference point measured along the curved wheel surface
      ! deltcp_rw  [rad]  roll angle from right wheel vertical to contact reference normal direction
      ! dirnam            optional working folder for experiment relative to the program's working
      !                   folder (set in contact.f90)
      ! expnam            experiment name, stand-alone program: excluding the working directory

   end type t_metadata

   !---------------------------------------------------------------------------------------------------------

   ! scaling factors for unit conversions in the CONTACT library routines:

   type :: t_scaling
      integer            :: units
      real(kind=8)       :: len
      real(kind=8)       :: area
      real(kind=8)       :: forc
      real(kind=8)       :: veloc
      real(kind=8)       :: angle
      real(kind=8)       :: body

      ! units          code of the scheme used in the API of the library routines
      ! len            scaling factor for lengths:    internal [mm]   = len   * external
      ! area           scaling factor for areas:      internal [mm2]  = area  * external
      ! forc           scaling factor for forces:     internal [N]    = forc  * external
      ! veloc          scaling factor for velocities: internal [mm/s] = veloc * external
      ! angle          scaling factor for angles:     internal [rad]  = angle * external
      ! body           scaling factor for signs:      internal [on 1] = body  * external
      ! (temperature: offset + scaling:  internal [*C] = temp * (external - offset))

   end type t_scaling

   !---------------------------------------------------------------------------------------------------------

   ! aggregate of all data for a CONTACT calculation in module 3:

   type :: t_probdata
      type(t_metadata) :: meta   ! meta-data describing the calculation
      type(t_scaling)  :: scl    ! scaling factors for the CONTACT library
      type(t_ic)       :: ic     ! integer control digits
      type(t_material) :: mater  ! material-description of the bodies
      type(t_potcon)   :: potcon ! description of potential contact area
      type(t_grid)     :: cgrid  ! main discretisation grid for CONTACT
      type(t_hertz)    :: hertz  ! Hertzian problem-description
      type(t_geomet)   :: geom   ! geometry-description of the bodies
      type(t_friclaw)  :: fric   ! input parameters of friction law used
      type(t_kincns)   :: kin    ! kinematic description of the problem
      type(t_influe)   :: influ  ! combined influence coefficients
      type(t_solvers)  :: solv   ! variables related to solution algorithms
      type(t_output)   :: outpt1 ! solution and derived quantities
      type(t_subsurf)  :: subs   ! data of subsurface stress calculation
   end type t_probdata

   type :: p_probdata
      type(t_probdata), pointer :: gd => NULL() ! pointer to a gd hierarchical data-structure
   end type p_probdata

contains

!------------------------------------------------------------------------------------------------------------

   subroutine ic_init( ic )
!--purpose: Initialize control digits to sensible default values
      implicit none
!--subroutine arguments:
      type(t_ic) :: ic

      ic%config = 1
      ic%pvtime = 2
      ic%bound  = 0
      ic%tang   = 3
      ic%norm   = 1
      ic%force  = 0
      ic%heat   = 0
      ic%stress = 0
      ic%sens   = 1
      ic%discns_inp  = 2
      ic%discns_eff  = ic%discns_inp
      ic%gapwgt      = 0
      ic%frclaw_inp  = 0
      ic%gencr_inp   = 2
      ic%mater  = 0
      ic%rznorm = 2
      ic%rztang = 0
      ic%gausei_inp  = 1
      ic%iestim = 0
      ic%matfil_surf = 1
      ic%matfil_subs = 1
      ic%output_surf = 2
      ic%output_subs = 1
      ic%flow   = 4
      ic%nmdbg  = 0
      ic%wrtinp = 0
      ic%return = 1
      ic%print_pmax = .true.
      ic%ilvout = 1

   end subroutine ic_init

!------------------------------------------------------------------------------------------------------------

   subroutine ic_unpack (modul, cpbtnfs, vldcmze, xhgiaowr, ic)
!--purpose: unpack the control-words into array ic
      implicit none
!--subroutine parameters:
      type(t_ic) :: ic
      integer    :: modul, cpbtnfs, vldcmze, xhgiaowr
!--local variables:
      integer, parameter :: idebug = 0
      integer            :: ihulp, vdigit, zdigit, edigit

      ihulp         = 0
      ic%config     = (cpbtnfs - ihulp) / 1000000
         ihulp      = ihulp + 1000000 * ic%config
      ic%pvtime     = (cpbtnfs - ihulp) / 100000
         ihulp      = ihulp +  100000 * ic%pvtime
      ic%bound      = (cpbtnfs - ihulp) / 10000
         ihulp      = ihulp +   10000 * ic%bound
      ic%tang       = (cpbtnfs - ihulp) / 1000
         ihulp      = ihulp +    1000 * ic%tang
      ic%norm       = (cpbtnfs - ihulp) / 100
         ihulp      = ihulp +     100 * ic%norm
      ic%force      = (cpbtnfs - ihulp) / 10
         ihulp      = ihulp +      10 * ic%force
      ic%stress     = (cpbtnfs - ihulp) / 1
      if (idebug.ge.2) write(*,'(7(a,i2))') 'unpack: C=',ic%config,', P=',ic%pvtime,', B=',ic%bound,    &
         ', T=',ic%tang,', N=',ic%norm,', F=',ic%force,', S=',ic%stress


      ihulp         = 0
      ic%varfrc     = (vldcmze - ihulp) / 1000000
         ihulp      = ihulp + 1000000 * ic%varfrc
      ic%frclaw_inp = (vldcmze - ihulp) / 100000
         ihulp      = ihulp +  100000 * ic%frclaw_inp
      ic%discns_inp = (vldcmze - ihulp) / 10000
         ihulp      = ihulp +   10000 * ic%discns_inp
      ic%gencr_inp  = (vldcmze - ihulp) / 1000
         ihulp      = ihulp +    1000 * ic%gencr_inp
      ic%mater      = (vldcmze - ihulp) / 100
         ihulp      = ihulp +     100 * ic%mater
      zdigit        = (vldcmze - ihulp) / 10
         ihulp      = ihulp +      10 * zdigit
      edigit        = (vldcmze - ihulp) / 1
      if (idebug.ge.2) write(*,'(7(a,i2))') 'unpack: V=',vdigit,', L=',ic%frclaw_inp,                   &
         ', D=',ic%discns_inp, ', C=',ic%gencr_inp,', M=',ic%mater, ', Z=',zdigit,', E=',edigit

      if (modul.eq.1) then
         ic%ztrack = zdigit
         ic%ewheel = edigit
         ic%rztang = 0
         ic%rznorm = 0
      else
         ic%ztrack = 0
         ic%ewheel = 0
         ic%rznorm = zdigit
         ic%rztang = edigit
      endif

      ihulp          = 0
      ic%nmdbg       = (xhgiaowr  - ihulp) / 10000000
         ihulp       = ihulp + 10000000 * ic%nmdbg
      ic%heat        = (xhgiaowr  - ihulp) / 1000000
         ihulp       = ihulp +  1000000 * ic%heat
      ic%gausei_inp  = (xhgiaowr  - ihulp) / 100000
         ihulp       = ihulp +   100000 * ic%gausei_inp
      ic%iestim      = (xhgiaowr  - ihulp) / 10000
         ihulp       = ihulp +    10000 * ic%iestim
      ic%matfil_surf = (xhgiaowr  - ihulp) / 1000
         ihulp       = ihulp +     1000 * ic%matfil_surf
      ic%output_surf = (xhgiaowr  - ihulp) / 100
         ihulp       = ihulp +      100 * ic%output_surf
      ic%flow        = (xhgiaowr  - ihulp) / 10
         ihulp       = ihulp +       10 * ic%flow
      ic%return      = (xhgiaowr  - ihulp)

      if (idebug.ge.2) write(*,'(8(a,i2))') 'unpack: X=',ic%nmdbg,', H=',ic%heat, ', G=',ic%gausei_inp,  &
         ', I=',ic%iestim,', A=',ic%matfil_surf,', O=',ic%output_surf,', W=',ic%flow,  ', R=',ic%return

   end subroutine ic_unpack

!------------------------------------------------------------------------------------------------------------

   subroutine ic_pack (modul, cpbtnfs, vldcmze, xhgiaowr, ic)
!--purpose: pack the controls in Ic together in the control-words
      implicit none
!--subroutine parameters:
      type(t_ic) :: ic
      integer    :: modul, cpbtnfs, vldcmze, xhgiaowr

      cpbtnfs = 0
      if (modul.eq.1) cpbtnfs =  1000000*ic%config

      cpbtnfs =  cpbtnfs +                                                                              &
                    100000*ic%pvtime      +    10000*ic%bound       +    1000*ic%tang        +          &
                       100*ic%norm        +       10*ic%force       +       1*ic%stress

      if (modul.eq.1) then
         vldcmze =                                                    1000000*ic%varfrc      +          &
                    100000*ic%frclaw_inp  +    10000*ic%discns_inp  +    1000*ic%gencr_inp   +          &
                       100*ic%mater       +       10*ic%ztrack      +       1*ic%ewheel

      else
         vldcmze =  100000*ic%frclaw_inp  +    10000*ic%discns_inp  +    1000*ic%gencr_inp   +          &
                       100*ic%mater       +       10*ic%rznorm      +       1*ic%rztang
      endif

      xhgiaowr =                            10000000*ic%nmdbg       + 1000000*ic%heat        +          &
                    100000*ic%gausei_inp  +    10000*ic%iestim      +    1000*ic%matfil_surf +          &
                       100*ic%output_surf +       10*ic%flow        +       1*ic%return

   end subroutine ic_pack

!------------------------------------------------------------------------------------------------------------

   subroutine meta_init( m, ilevel )
!--purpose: Initialize meta-data to sensible default values
   implicit none
!--subroutine arguments:
   type(t_metadata) :: m
   integer          :: ilevel        ! 1 = time-varying parts; 2 = full init

   if (ilevel.ge.2) then
      m%REid      = 0
      m%CPid      = 0
      m%actv_thrd = -1
      m%irun      = 0
      m%iax       = 0
      m%iside     = 0
      m%ncase     = 0
      m%itbrent   = 0
      m%npatch    = 0
      m%ipatch    = 0
      m%tim       = 0d0
      m%dirnam    = ' '
      m%expnam    = ' '
   endif

   if (ilevel.ge.1) then
      m%s_ws      = 0d0
      m%ynom_whl  = 0d0
      m%rnom_whl  = 0d0
      m%rnom_rol  = 0d0

      m%x_rw      = 0d0
      m%y_rw      = 0d0
      m%z_rw      = 0d0
      m%rollrw    = 0d0
      m%yawrw     = 0d0

      m%y_rr      = 0d0
      m%z_rr      = 0d0
      m%rollrr    = 0d0

      m%xcp_tr    = 0d0
      m%ycp_tr    = 0d0
      m%zcp_tr    = 0d0
      m%deltcp_tr = 0d0

      m%xcp_rr    = 0d0
      m%ycp_rr    = 0d0
      m%zcp_rr    = 0d0
      m%scp_rr    = 0d0
      m%deltcp_rr = 0d0

      m%xcp_rw    = 0d0
      m%ycp_rw    = 0d0
      m%zcp_rw    = 0d0
      m%scp_rw    = 0d0
      m%deltcp_rw = 0d0
   endif

   end subroutine meta_init

!------------------------------------------------------------------------------------------------------------

   subroutine mater_init( m, ic )
!--purpose: Initialize material-data to sensible default values
   implicit none
!--subroutine arguments:
   type(t_material) :: m
   type(t_ic)       :: ic

   m%gg      = (/  82000.d0, 82000.d0 /)
   m%poiss   = (/    0.28d0,   0.28d0 /)
   m%flx_z   = 0d0
   m%flx     = (/ 0.00001d0, 0.00001d0, 0d0 /)
   m%bktemp  = (/       0d0,      0d0 /)
   m%heatcp  = (/     450d0,    450d0 /)
   m%lambda  = (/      50d3,     50d3 /)
   m%dens    = (/   7850d-9,  7850d-9 /)
   m%betapl  = 1d0
   m%k0_mf   = 1d0
   m%alfamf  = 1d0
   m%betamf  = 1d0
   m%gg3     = 82000.d0
   m%laythk  = 0d0
   m%tau_c0  = 1d20
   m%k_tau   = 0d0
   m%bound_eff = ic%bound
   m%mater_eff = ic%mater
   m%gencr_eff = ic%gencr_inp
   m%fname_influe  = ' '
   m%if_meth = 0
   m%if_ver  = 4
   m%ninclin = 0

   end subroutine mater_init

!------------------------------------------------------------------------------------------------------------

   subroutine combin_mater(m)
!--purpose: update the combined material parameters
      implicit none
!--subroutine arguments:
      type(t_material) :: m

      m%ga = 2d0 / (1d0/m%gg(1) + 1d0/m%gg(2))
      m%nu = m%ga * ( m%poiss(1)/m%gg(1) + m%poiss(2)/m%gg(2) ) / 2d0
      m%ak = (m%ga/4d0) * ( (1d0-2d0*m%poiss(1))/m%gg(1) - (1d0-2d0*m%poiss(2))/m%gg(2) )

   end subroutine combin_mater

!------------------------------------------------------------------------------------------------------------

   subroutine inflcf_mater(inflcf, mater)
!--purpose: update the material parameters for an influence coefficients matrix
      implicit none
!--subroutine arguments:
      type(t_inflcf)  , intent(inout) :: inflcf
      type(t_material), intent(in)    :: mater

      inflcf%ga       = mater%ga
      inflcf%ga_inv   = 1d0 / mater%ga
      inflcf%use_3bl  = (mater%mater_eff.eq.4 .and. mater%laythk.ge.1d-10)
      inflcf%flx_3bl  = mater%flx(1)
      inflcf%use_flxz = (mater%bound_eff.eq.1)
      inflcf%flx_z    = mater%flx_z
   end subroutine inflcf_mater

!------------------------------------------------------------------------------------------------------------

   subroutine influe_mater(influe, mater)
!--purpose: update the material parameters for the full influence coefficients structure
      implicit none
!--subroutine arguments:
      type(t_influe)  , intent(inout) :: influe
      type(t_material), intent(in)    :: mater

      call inflcf_mater(influe%cs,  mater)
      call inflcf_mater(influe%cv,  mater)
      call inflcf_mater(influe%csv, mater)
      call inflcf_mater(influe%ms,  mater)
   end subroutine influe_mater

!------------------------------------------------------------------------------------------------------------

   subroutine subsurf_copy(s_in, s_out)
!--purpose: copy input data for the subsurface stress calculation
      implicit none
!--subroutine arguments:
      type(t_subsurf)  , intent(in)    :: s_in
      type(t_subsurf)  , intent(inout) :: s_out
!--local variables:
      integer, parameter :: idebug = 0
      integer            :: iblk

      if (idebug.ge.1 .and. s_out%nblock.ge.1) then
         write(bufout,*) 's_out has',s_out%nblock,' blocks of points, discarding'
         call write_log(1, bufout)
      endif

      ! clean-up contents of s_out

      do iblk = 1, s_out%nblock
         call subsblk_destroy(s_out%blocks(iblk))
      enddo

      ! copy contents of s_in

      if (idebug.ge.1 .and. s_in%nblock.ge.1) then
         write(bufout,*) 's_in has',s_in%nblock,' blocks of points, copying'
         call write_log(1, bufout)
      endif

      do iblk = 1, s_in%nblock
         call subsblk_copy(s_in%blocks(iblk), s_out%blocks(iblk))
      enddo

      s_out%nblock = s_in%nblock

   end subroutine subsurf_copy

!------------------------------------------------------------------------------------------------------------

   subroutine subsblk_copy(b_in, b_out)
!--purpose: copy input data for a block of subsurface points
      implicit none
!--subroutine arguments:
      type(t_subsblk)  , intent(in)    :: b_in
      type(t_subsblk)  , intent(inout) :: b_out

      b_out%isubs   = b_in%isubs
      b_out%ncolum  = b_in%ncolum
      b_out%nx_inp  = b_in%nx_inp
      b_out%ixl_inp = b_in%ixl_inp
      b_out%ixinc   = b_in%ixinc
      b_out%ixh_inp = b_in%ixh_inp
      b_out%nx_eff  = b_in%nx_eff
      b_out%ixl_eff = b_in%ixl_eff
      b_out%ixh_eff = b_in%ixh_eff
      b_out%ny_inp  = b_in%ny_inp
      b_out%iyl_inp = b_in%iyl_inp
      b_out%iyinc   = b_in%iyinc
      b_out%iyh_inp = b_in%iyh_inp
      b_out%ny_eff  = b_in%ny_eff
      b_out%iyl_eff = b_in%iyl_eff
      b_out%iyh_eff = b_in%iyh_eff
      b_out%nz      = b_in%nz
      b_out%zl      = b_in%zl
      b_out%dz      = b_in%dz

      call copy_1d_int(b_in%ixlist_inp, b_out%ixlist_inp)
      call copy_1d_int(b_in%iylist_inp, b_out%iylist_inp)
      call copy_1d_int(b_in%ixlist_eff, b_out%ixlist_eff)
      call copy_1d_int(b_in%iylist_eff, b_out%iylist_eff)

      call copy_1d_real(b_in%x, b_out%x)
      call copy_1d_real(b_in%y, b_out%y)
      call copy_1d_real(b_in%z, b_out%z)

      if (associated(b_out%table))  deallocate(b_out%table)
      b_out%table  => NULL()

   end subroutine subsblk_copy

!------------------------------------------------------------------------------------------------------------

   subroutine subsblk_destroy(blk)
!--purpose: de-allocate space allocated for a block of subsurface data
      implicit none
!--subroutine arguments:
      type(t_subsblk)  , intent(inout) :: blk

      if (associated(blk%ixlist_inp)) deallocate(blk%ixlist_inp)
      if (associated(blk%iylist_inp)) deallocate(blk%iylist_inp)
      if (associated(blk%ixlist_eff)) deallocate(blk%ixlist_eff)
      if (associated(blk%iylist_eff)) deallocate(blk%iylist_eff)
      if (associated(blk%x))          deallocate(blk%x)
      if (associated(blk%y))          deallocate(blk%y)
      if (associated(blk%z))          deallocate(blk%z)
      if (associated(blk%table))      deallocate(blk%table)
      blk%ixlist_inp => NULL()
      blk%iylist_inp => NULL()
      blk%ixlist_eff => NULL()
      blk%iylist_eff => NULL()
      blk%x          => NULL()
      blk%y          => NULL()
      blk%z          => NULL()
      blk%table      => NULL()

   end subroutine subsblk_destroy

!------------------------------------------------------------------------------------------------------------

   subroutine solv_init ( solv )
!--purpose: Set appropriate initial values for the solver settings
      implicit none
!--subroutine arguments:
      type(t_solvers) :: solv

      solv%gausei_eff = 0
      solv%maxgs   =  999
      solv%maxin   =   20
      solv%maxnr   =   25
      solv%maxout  =    1
      solv%mxsens  =   30
      solv%eps     = 1d-5
      solv%epsens  = 1d-4
      solv%omegah  = 0.90
      solv%omegas  = 1.00
      solv%gd_meth =    1
      solv%kdown   =    3
      solv%fdecay  = 0.90
      solv%betath  = 0.10
      solv%kdowfb  =    1
      solv%d_ifc   = 1.00
      solv%d_lin   = 1.00
      solv%d_cns   = 1.00
      solv%d_slp   = 1.00
      solv%pow_s   = 0.90
      solv%inislp  =    0
      solv%omgslp  = 0.50

   end subroutine solv_init

!------------------------------------------------------------------------------------------------------------

   subroutine potcon_init ( potcon )
!--purpose: Set appropriate initial values for the potential contact area
      implicit none
!--subroutine arguments:
      type(t_potcon) :: potcon

      potcon%ipotcn = 1
      potcon%mx   = 22
      potcon%my   = 20
      potcon%npot = potcon%mx * potcon%my
      potcon%dx   = 0.5d0
      potcon%dy   = 0.5d0
      potcon%xl   = -5.0d0
      potcon%yl   = -5.0d0

   end subroutine potcon_init

!------------------------------------------------------------------------------------------------------------

   subroutine geom_init ( geom )
!--purpose: Set appropriate initial values for the geometry description
      implicit none
!--subroutine arguments:
      type(t_geomet) :: geom
!--local variables:
      integer j

      geom%ibase  = 1
      geom%iplan  = 1
      allocate(geom%prmudf(10))
      geom%prmudf(1:10)= (/ 0.0056566, 0.0, 0.0056566, (0.0, j=4,10) /)
      allocate(geom%prmpln(10))
      geom%prmpln(1:10)= (/ (0.0, j=1,10) /)

   end subroutine geom_init

!------------------------------------------------------------------------------------------------------------

   subroutine kincns_init ( kin, ic, fric, dx )
!--purpose: Set appropriate initial values for the kinematic data
      implicit none
!--subroutine arguments:
      type(t_kincns)  :: kin
      type(t_ic)      :: ic
      type(t_friclaw) :: fric
      real(kind=8)    :: dx

      kin%dt      =  1d0             ! for shifts
      kin%veloc   =  30000d0         ! 30 m/s
      kin%dq      =  dx
      kin%facphi  =  1d0 / 6d0
      kin%chi     =  0.0d0
      kin%fntrue  =  100000d0
      kin%fnscal  =  1.2195d0        ! 100kN / 82kN/mm2
      kin%pen     =  0.1035d0        ! Rx=500, Ry=300mm
      kin%penv    =  0.0d0
      kin%fxrel1  = -0.875d0
      kin%fyrel1  =  0.0d0
      kin%cksi    =  0.0025d0
      kin%ceta    =  0.0d0
      kin%cphi    =  0.0d0
      kin%fprev   = (/ kin%fxrel1, kin%fyrel1, kin%fntrue /)
      kin%use_muscal = ic%varfrc.eq.0
      kin%muscal  =  1.0d0
      if (kin%use_muscal) kin%muscal = fric%fstat()

   end subroutine kincns_init

!------------------------------------------------------------------------------------------------------------

   subroutine hertz_init ( hertz )
!--purpose: Set appropriate initial values for the Hertzian data
      implicit none
!--subroutine arguments:
      type(t_hertz) :: hertz

      hertz%aa    =  7.5d0
      hertz%bb    =  5.3d0
      hertz%bneg  =  hertz%bb
      hertz%bpos  =  hertz%bb
      hertz%a1    =  0.001           ! Rx=500mm
      hertz%b1    =  0.001667        ! Ry=300mm
      hertz%aob   =  1.4d0
      hertz%scale =  1.1d0

   end subroutine hertz_init

!------------------------------------------------------------------------------------------------------------

   subroutine setini ( gd )
!--purpose: Set appropriate initial values for the control/input-variables
      implicit none
!--subroutine arguments:
      type(t_probdata) :: gd
!--local variables:
      real(kind=8) :: xc1, yc1

      call meta_init  ( gd%meta, 2 )
      call ic_init    ( gd%ic )
      call solv_init  ( gd%solv )
      call mater_init ( gd%mater, gd%ic )
      call potcon_init( gd%potcon )

      xc1 = gd%potcon%xl + 0.5d0*gd%potcon%dx
      yc1 = gd%potcon%yl + 0.5d0*gd%potcon%dy
      call grid_create_uniform(gd%cgrid, nxarg=gd%potcon%mx, x0arg=xc1, dxarg=gd%potcon%dx,             &
                               nyarg=gd%potcon%my, y0arg=yc1, dyarg=gd%potcon%dy, zarg=0d0)

      call geom_init  ( gd%geom )
      call fric_init  ( gd%fric )
      call kincns_init( gd%kin, gd%ic, gd%fric, gd%cgrid%dx )
      call hertz_init ( gd%hertz )

      call gf3_new(gd%outpt1%ps, 'outpt1%ps', gd%cgrid)
      call gf3_new(gd%outpt1%pv, 'outpt1%pv', gd%cgrid)

      gd%subs%nblock = 0

   end subroutine setini

!------------------------------------------------------------------------------------------------------------

   subroutine gd_destroy(gd)
!--function: destroy the contact hierarchical data-structure gd
   implicit none
!--subroutine arguments:
   type(t_probdata) :: gd
!--local variables:
   integer          :: iblk

   ! cleanup of grid:

   call grid_destroy(gd%cgrid)

   ! allocatable & pointer arrays of t_geomet:

   call destroy_arr(gd%geom%prmudf)
   call destroy_arr(gd%geom%prmpln)
   call destroy_arr(gd%geom%ysep)
   call destroy_arr(gd%geom%facsep)
   call gf3_destroy(gd%geom%exrhs)
   call gf3_destroy(gd%geom%hs1)
   call gf3_destroy(gd%geom%hv1)

   ! TODO: cleanup subsurf.input
   !! call destroy_arr(gd%geom%xb)

   call fric_destroy( gd%fric )

   ! cleanup of influence coefficients influ:

   call inflcf_destroy(gd%influ%cs)
   call inflcf_destroy(gd%influ%cv)
   call inflcf_destroy(gd%influ%csv)
   call inflcf_destroy(gd%influ%ms)

   ! cleanup of output arrays outpt1:

   call eldiv_destroy(gd%outpt1%igs)
   call eldiv_destroy(gd%outpt1%igv)
   call destroy_arr(gd%outpt1%ledg%jbnd)
   call destroy_arr(gd%outpt1%ledg%ixbnd)
   call if1_destroy(gd%outpt1%ledg%ii2j)
   call destroy_arr(gd%outpt1%ledg%xbnd)
   call destroy_arr(gd%outpt1%ledg%facdx)
   call destroy_arr(gd%outpt1%ledg%ubnd)
   call gf3_destroy(gd%outpt1%ledg%facdt)
   call gf3_destroy(gd%outpt1%mus)
   call gf3_destroy(gd%outpt1%muv)
   call gf3_destroy(gd%outpt1%shft)
   call gf3_destroy(gd%outpt1%ps)
   call gf3_destroy(gd%outpt1%pv)
   call gf3_destroy(gd%outpt1%us)
   call gf3_destroy(gd%outpt1%uv)
   call gf3_destroy(gd%outpt1%ss)
   call gf3_destroy(gd%outpt1%sv)
   call gf3_destroy(gd%outpt1%taucs)
   call gf3_destroy(gd%outpt1%taucv)
   call gf3_destroy(gd%outpt1%upls)
   call gf3_destroy(gd%outpt1%uplv)
   call gf3_destroy(gd%outpt1%temp1)
   call gf3_destroy(gd%outpt1%temp2)

   ! cleanup subsurf data

   do iblk = 1, gd%subs%nblock
      call subsblk_destroy( gd%subs%blocks(iblk) )
   enddo

   end subroutine gd_destroy

!------------------------------------------------------------------------------------------------------------

end module m_hierarch_data
