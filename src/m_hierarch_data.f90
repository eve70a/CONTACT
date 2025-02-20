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
use m_interp
use m_spline
use m_grids
use m_gridfunc
use m_interp_gf
use m_spline_get
use m_bspline
use m_friclaw
use m_inflcf

implicit none

public

   public  t_ic
   public  ic_init
   public  ic_copy
   public  ic_unpack
   public  ic_pack
   public  ic_unpack_dbg
   public  ic_pack_dbg
   public  ic_destroy

   public  t_material
   public  mater_init
   public  mater_copy
   public  combin_mater
   public  inflcf_mater
   public  influe_mater
   public  mater_destroy

   public  t_potcon
   public  potcon_copy
   public  potcon_fill
   public  potcon_cgrid
   public  potcon_hertz
   public  potcon_get_overlap
   public  potcon_merge
   public  potcon_destroy

   public  t_hertz
   public  hertz_init
   public  hertz_copy
   public  hertz_destroy

   public  t_geomet
   public  geom_init
   public  geom_copy
   public  geom_destroy

   public  t_kincns
   public  kincns_init
   public  kincns_copy
   public  kincns_destroy

   public  t_solvers
   public  solv_init
   public  solv_copy
   public  solv_destroy

   public  t_output
   public  output_init
   public  output_copy
   public  output_destroy

   public  t_subsblk
   public  subsblk_copy
   public  subsblk_destroy

   public  t_subsurf
   public  subsurf_init
   public  subsurf_copy
   public  subsurf_destroy

   public  t_metadata
   public  meta_init
   public  meta_copy
   public  meta_destroy

   public  t_scaling
   public  scaling_copy
   public  scaling_destroy

   public  t_probdata
   public  p_probdata
   public  gd_init
   public  gd_copy
   public  gd_resize_gridfunc
   public  gd_merge_gridfunc
   public  gd_merge
   public  gd_destroy

   !---------------------------------------------------------------------------------------------------------

   ! codes for the different solvers

   integer, parameter :: isolv_fastrp = -1, isolv_fastsm =  0, isolv_normcg =  1,                       &
                         isolv_tangcg =  2, isolv_cnvxgs =  3, isolv_stdygs =  4,                       &
                         isolv_gdstdy =  7

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
      integer :: force1
      integer :: force3
      integer :: heat
      integer :: stress
      integer :: sens
      integer :: varfrc
      integer :: frclaw_inp
      integer :: discns1_inp
      integer :: discns1_eff
      integer :: discns3
      integer :: gapwgt
      integer :: gencr_inp
      integer :: mater
      integer :: mater2
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
      integer :: xflow
      integer :: x_profil
      integer :: x_smooth
      integer :: x_force
      integer :: x_locate
      integer :: x_cpatch
      integer :: x_inflcf
      integer :: x_nmdbg
      integer :: x_readln
      integer :: wrtinp
      integer :: return
      integer :: ilvout
      logical :: print_pmax
      logical :: print_uxavg

   contains
      procedure :: is_left_side     => ic_is_left_side
      procedure :: is_roller        => ic_is_roller
      procedure :: is_conformal     => ic_is_conformal
      procedure :: use_initial_cp   => ic_use_initial_cp
      procedure :: use_oblique      => ic_use_oblique
      procedure :: use_steep_slopes => ic_use_steep_slopes
      procedure :: use_supergrid    => ic_use_supergrid

      ! c1, config  configuration or composition of wheel/rail problems
      !              0 = wheelset on track, left wheel,
      !              1 = wheelset on track, right wheel,
      !              4 = wheelset on rollers, left wheel,
      !              5 = wheelset on rollers, right wheel
      ! p, pvtime   the filling of the tractions Pv for the previous time instance,
      !              0 = new timestep in a sequence; copy result Ps to Pv; Fcntc to Fprev
      !             (1 = new timestep for the normal part only, set tangential to zero)
      !              2 = initiation of new sequence, no continuation, the entire Pv is cleared
      !              3 = new iteration for a timestep; keep Pv & Fprev unmodified. 
      !             note: previous tractions have no effect for T=0 and 3.
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
      ! f1, force1, module 1: specifies the long/lat forces:
      !              0 = no rail deflection, pitch velocity omega_ws prescribed;
      !              1 = no rail deflection, total long. force fx_ws prescribed;
      !              2 = no rail deflection, total long. moment my_ws prescribed;
      !              3 = rail deflection with stiffness ky, kz and spring forces fy_rail, fz_rail given
      ! f3, force3, module 3: specifies the tangential problem:
      !              0 = creepages cksi, ceta prescribed;
      !              1 = relative force  fxrel, creepage ceta prescribed;
      !              2 = relative forces fxrel, fyrel prescribed
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
      !              1 = multiple inputs with linear interpolation for lateral direction across rail head
      !              2 = multiple inputs with linear interpolation for longitudinal direction along track
      !              3 = multiple sets of inputs are used, per row (iy) of the potential contact
      !                  used in module 3 for conformal contacts with varfrc=1 in module 1
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
      ! d1, discns1 concerns the potential contact area and discretisation for module 1:
      !              0, 1 = maintain the same approach and parameters as used in the previous case
      !              2 = "wgt-wgt-2" planar contact with contact locus, weighed center and weighted angle,
      !                  read new input parameters and form the discretisation
      !              3 = "wx-surf",    same as 2, with extended rigid slip computation
      !              4 = "curved",     same as 2, with curved reference and conformal contact computation
      !              5 = "brute",      same as 2, using brute force instead of contact locus
      !              6 = "brute+steep" same as 5, using larger grid to accomodate steep slopes
      !              7 = "icp-loc",    same as 2, using initial contact point instead of weighted center,
      !                                           using the local contact angle instead of a weighted value
      !              8 = "locus+view"  same as 2, using a hard-coded oblique view direction
      !              9 = "brute+view"  same as 5, using a hard-coded oblique view direction
      !             Note: the variable ic%discns1_inp contains the value obtained from (and written to) the
      !                   input-file, whereas ic%discns1_eff is the effective method used (2 -- 9).
      ! d3, discns3 concerns the potential contact area and discretisation for module 3:
      !              0 = maintain the discretisation grid of the previous case;
      !              1 = form the discretisation grid from parameters in storage;
      !              2 = read new parameters from input-file and form the discretisation
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
      !              5 = FaStrip method, combination of strip theory + simplified theory
      !              6 = reserved (pseudo-viscous damping)
      !              7 = reserved (vertical slice/gap)
      !              Note: the m-digit is copied to t_material
      ! m2, mater2  type of damping model to be used:
      !              0 = no damping
      !              1 = ad-hoc proportional damping with parameters from storage
      !              2 = ad-hoc proportional damping with new parameters from input
      ! z1, ztrack  concerns the track geometry, profile, deviation, and deflection parameters (module 1)
      !              0 = maintain track dimensions, profile, deviations, and deflection parameters;
      !              1 = read new dimensions and profile;
      !              2 = read new track deviations and deflection parameters (if F=3);
      !              3 = read new dimensions, profile, track deviations, and deflection parameters if F=3
      ! e1, ewheel   concerns the wheel-set geometry, profile, position and velocity (module 1)
      !              0 = maintain wheel-set geometry, profile, position and velocity;
      !              1 = read new position data;
      !              2 = read new position and velocity data;
      !              3 = read new position and velocity data, geometry, wheel profile
      !              4 = read new position and velocity data, including flexible wheelset deviations
      !              5 = as 3, including flexible wheelset deviations
      ! z3, rznorm   concerns the right hand side of the normal problem (module 3)
      !              0 = maintain undeformed distance and planform;
      !              1 = form undeformed distance from parameters in storage;
      !              2 = read new parameters and compute undeformed distance.
      ! e3, rztang/exrhs  extra term of the rigid slip, rhs of tangential problem (module 3)
      !              0 = set the extra term equal to zero;
      !              1 = maintain the extra term of the previous case;
      !              2 = add extra term w_y = -x\phi, removing spin contribution from y-direction;
      !              3 = add extra term w_x =  y\phi, removing spin contribution from x-direction;
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
      ! a, matfil_surf   governs the use of the matlab-file <experim>.<case>.mat per case:
      !             -2 = results are evaluated for all points of the potential contact area (soutpt),
      !                  no mat-files created  -- for internal use in module 1 calling module 3
      !             -1 = results are evaluated for points in the actual contact area (soutpt),
      !                  no mat-files created  -- for internal use in module 1 calling module 3
      !              0 = no mat-files created;
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
      ! x, xflow    governs additional debug outputs printed to the output-file
      !              0   = none,
      !              1   = additional debug output enabled
      !    x_profil governs the debug output on profile processing
      !    x_smooth governs the debug output on profile smoothing
      !    x_force  governs the debug output on total force iteration
      !    x_locate governs the debug output on the contact search
      !    x_cpatch governs the debug output on shuffling of contact patches
      !    x_inflcf governs the debug output on influence coefficients
      !              0   = none
      !              7   = print influence coefficients
      !    x_nmdbg  governs the debug output from the contact solvers
      !              0-3 = none
      !              4   = global problem inputs Hs and outputs Igs,Ps,Us
      !              5   = intermediate results of the Panagiotopoulos' process
      !              6   = intermediate results of Duvorol and Veloc.dep iterations
      !              7   = intermediate results of NORM and TANG algorithms
      !              8   = intermediate results of Newton-Raphson loops
      !              9   = information per iteration of NormCG, ConvexGS, SteadyGS
      !    x_readln governs the debug output on reading the input file
      !    wrtinp   governs the writing of data to the input-file, used by the CONTACT library,
      !             allowing for off-line analysis of cases using the stand-alone CONTACT program
      !              0   = no writing of the input (default)
      !              1   = write input to the .inp-file using actual control digits
      !              2   = write input to inp-file, full configuration
      ! r, return   return to main program.
      !              0 = calculate, stay in module, 1 = calculate and return,
      !              2 = skip calculation, stay, 3 = skip calculation and return
      ! print_pmax  .true.: compute/show pmax in aggregate results per patch
      ! print_uxavg .true.: compute/show ux_avg in aggregate results per patch
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
      real(kind=8) :: cdampn
      real(kind=8) :: cdampt
      real(kind=8) :: dfnmax
      real(kind=8) :: dftmax
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

      ! cdampn  [s]    damping coefficient for normal contact forces
      ! cdampt  [s]    damping coefficient for tangential contact forces
      ! dfnmax [N/s]   maximum value for dFn/dt used in damping of normal contact forces
      ! dftmax [N/s]   maximum value for dFt/dt used in damping of tangential contact forces

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
      real(kind=8),    dimension(:),   allocatable :: prmudf 
      real(kind=8),    dimension(:),   allocatable :: prmpln 
      real(kind=8),    dimension(:),   allocatable :: prmrig
      real(kind=8),    dimension(:,:), allocatable :: xylim     ! (npatch,4)
      real(kind=8),    dimension(:,:), allocatable :: facsep    ! (npatch,npatch)
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
      !                data stored in xylim(1:np,4), facsep(1:np,1:np).
      ! prmudf [?]     parameters for computing the undeformed distance
      ! prmpln [?]     coefficients of the planform
      ! prmrig [?]     parameters for computing the rigid slip
      ! xylim  [mm]    extent [xl,xh] x [yl,yh] of (sub-)patches in the potential contact area
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
      real(kind=8) :: spinxo
      real(kind=8) :: spinyo
      real(kind=8) :: fntrue
      real(kind=8) :: fxrel
      real(kind=8) :: fyrel
      real(kind=8) :: fcntc(3)
      real(kind=8) :: fprev(3)
      real(kind=8) :: fdamp(3)
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
      !                TODO: it would be better to have separate variables for creep & shift
      ! spinxo [mm]    x-component of 'spin center' (xo,yo), spin creepage linearization point
      ! spinyo [mm]    y-component of 'spin center' (xo,yo), spin creepage linearization point

      ! fntrue [N]     total normal force, input when N=1, output when N=0
      ! fxrel  [-]     total tangential force in x-direction _on_ body (1), relative to muscal*fntrue;
      !                input when F3>=1, output when F3=0.
      ! fyrel  [-]     total tangential force in y-direction _on_ body (1), relative to muscal*fntrue;
      !                input when F3=2, output when F3<=1.
      ! fcntc  [N]     total forces [fxtrue,fytrue,fntrue] of current time instance
      ! fprev  [N]     total forces [fxtrue,fytrue,fntrue] of previous time instance
      ! fdamp  [N]     damping forces of current time instance

      ! use_muscal     flag indicating whether friction parameters are constant within the contact
      !                patch, permitting scaling of tangential forces by MU*FN (true) or by FN (false)
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

   ! data with respect to the solution and derived quantities:

   type :: t_output
      type(t_eldiv)    :: igs, igv
      type(t_gridfnc3) :: mus, muv, shft
      type(t_gridfnc3) :: ps, pv, us, uv, ss, sv
      type(t_gridfnc3) :: taucs, taucv, upls, uplv
      type(t_gridfnc3) :: temp1, temp2

      real(kind=8)     :: sens(nsens_out, nsens_in)
      real(kind=8)     :: mxtrue, mytrue, mztrue
      real(kind=8)     :: elen, frpow, pmax

      ! igs            element division for current solution, see t_eldiv.
      ! igv            element division for previous solution.

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
      !                  fntrue [N], fxrel, fyrel [-], mztrue [N.mm] etc.
      !                w.r.t. input parameters:
      !                  pen [mm], cksi, ceta [mm], cphi [rad] etc. (T=1)
      !                  pen [mm], cksi, ceta [-], cphi [rad/mm] etc. (T=2,3)
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
      integer            :: whl_ver, rail_ver
      integer            :: actv_thrd
      integer            :: irun, iax, iside, ncase, itforce, itforc_out, itforc_inn
      integer            :: npatch, ipatch
      real(kind=8)       :: tim, s_ws, th_ws, ynom_whl, rnom_whl, rnom_rol
      real(kind=8)       :: x_w, y_w, z_w, roll_w, yaw_w
      real(kind=8)       :: y_r, z_r, roll_r
      real(kind=8)       :: xcp_tr, ycp_tr, zcp_tr, deltcp_tr
      real(kind=8)       :: xcp_r, ycp_r, zcp_r, scp_r, deltcp_r
      real(kind=8)       :: xcp_w, ycp_w, zcp_w, scp_w, deltcp_w
      character(len=256) :: wrkdir, outdir, expnam

      ! REid              result element ID, used by CONTACT add-on
      ! CPid              contact patch ID, used by CONTACT add-on
      ! whl/rail_ver      counter for wheel/rail profile versions written to file
      ! actv_thrd         OpenMP thread number working on this gd, -1=none, used for locking
      ! irun              for Sentient: run number
      ! iax               for Sentient: axle number
      ! iside             for Sentient: side number
      ! ncase             case number for the (REid,CPid)
      ! itforce           total iteration number: reset in wr_contact, incremented in wr_contact_pos
      ! itforc_out        iteration number for total force outer iteration
      ! itforc_inn        iteration number for total force inner iteration
      ! npatch            number of contact patches, used by module 1
      ! ipatch            contact patch number, used by module 1
      ! tim        [s]    (SIMPACK) simulation time
      ! s_ws       [mm]   location of wheel-set CM along the track curve
      ! th_ws      [rad]  pitch angle of wheel-set CM wrt track curve
      ! ynom_whl   [mm]   lateral position of wheel origin in wheelset coordinates
      ! rnom_whl   [mm]   nominal wheel radius == vertical position of wheel origin in wheelset coords
      ! rnom_rol   [mm]   nominal roller radius
      ! [xyz]_w    [mm]   location of right wheel origin m_w in terms of track coordinates
      ! roll_w     [rad]  roll angle of right wheel (marker m_w) with respect to track coordinates
      ! yaw_w      [rad]  yaw angle of right wheel (marker m_w) with respect to track coordinates
      ! [yz]_r     [mm]   location of right rail origin m_r in terms of track coordinates
      ! roll_r     [rad]  roll angle of right rail marker m_r with respect to track coordinates
      ! [xyz]cp_tr [mm]   location of contact reference point m_ref (m_cp) in terms of track coordinates
      ! deltcp_tr  [rad]  contact angle: roll angle from track vertical to contact reference normal direction
      ! [xyz]cp_r  [mm]   location of contact reference point m_ref (m_cp) in terms of right rail coordinates
      ! scp_r      [mm]   position of the contact reference point measured along the curved rail surface
      ! deltcp_r   [rad]  roll angle from right rail vertical to contact reference normal direction
      ! [xyz]cp_w  [mm]   location of contact reference point m_ref (m_cp) in terms of right wheel coordinates
      ! scp_w      [mm]   position of the contact reference point measured along the curved wheel surface
      ! deltcp_w   [rad]  roll angle from right wheel vertical to contact reference normal direction
      ! [xy]o_spin [mm]   offset from super-grid pot.contact origin to contact reference position
      ! wrkdir            optional 'effective working folder' for experiment, can be an absolute path or
      !                   relative to the program's working folder
      ! outdir            optional folder for output files, can be an absolute path or relative to wrkdir
      ! expnam            experiment name, excluding any folder names

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
      type(t_potcon)   :: potcon_inp   ! input values for new potential contact area
      type(t_potcon)   :: potcon_cur   ! description of potential contact area used at current time
      type(t_grid)     :: cgrid_inp    ! discretisation grid used during input
      type(t_grid)     :: cgrid_cur    ! main discretisation grid for CONTACT
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
      ic%force1 = 0
      ic%force3 = 0
      ic%heat   = 0
      ic%stress = 0
      ic%sens   = 1
      ic%discns1_inp = 2
      ic%discns1_eff = ic%discns1_inp
      ic%discns3     = 2
      ic%gapwgt      = 0
      ic%frclaw_inp  = 0
      ic%gencr_inp   = 2
      ic%mater  = 0
      ic%mater2 = 0
      ic%rznorm = 2
      ic%rztang = 0
      ic%gausei_inp  = 1
      ic%iestim = 0
      ic%matfil_surf = 1
      ic%matfil_subs = 1
      ic%output_surf = 2
      ic%output_subs = 1
      ic%flow     = 4
      ic%xflow    = 0
      ic%x_profil = 1
      ic%x_smooth = 0
      ic%x_force  = 0
      ic%x_locate = 0
      ic%x_cpatch = 0
      ic%x_readln = 0
      ic%x_inflcf = 0
      ic%x_nmdbg  = 0
      ic%wrtinp   = 0
      ic%return   = 1
      ic%print_pmax  = .true.
      ic%print_uxavg = .false.
      ic%ilvout   = 1

   end subroutine ic_init

!------------------------------------------------------------------------------------------------------------

   subroutine ic_copy( ic_in, ic_out )
!--purpose: Copy control digits 
      implicit none
!--subroutine arguments:
      type(t_ic) :: ic_in, ic_out

      ic_out = ic_in                 ! plain structure, intrinsic copying

   end subroutine ic_copy

!------------------------------------------------------------------------------------------------------------

   function ic_is_left_side(this)
!--function: determine from the ic-struct whether a left-side configuration is used
      implicit none
!--result value
      logical                   :: ic_is_left_side
!--subroutine arguments
      class(t_ic),   intent(in) :: this

      ic_is_left_side = (this%config.eq.0 .or. this%config.eq.4)

   end function ic_is_left_side

!------------------------------------------------------------------------------------------------------------

   function ic_is_roller(this)
!--function: determine from the ic-struct whether a wheel-on-roller configuration is used
      implicit none
!--result value
      logical                   :: ic_is_roller
!--subroutine arguments
      class(t_ic),   intent(in) :: this

      ic_is_roller = (this%config.ge.4)

   end function ic_is_roller

!------------------------------------------------------------------------------------------------------------

   function ic_is_conformal(this)
!--function: determine from the ic-struct whether a conformal computation is used
      implicit none
!--result value
      logical                   :: ic_is_conformal
!--subroutine arguments
      class(t_ic),   intent(in) :: this

      ic_is_conformal = (this%discns1_eff.eq.4)

   end function ic_is_conformal

!------------------------------------------------------------------------------------------------------------

   function ic_use_initial_cp(this)
!--function: determine from the ic-struct whether to use the initial contact point as reference position
      implicit none
!--result value
      logical                   :: ic_use_initial_cp
!--subroutine arguments
      class(t_ic),   intent(in) :: this

      ic_use_initial_cp = (this%discns1_eff.eq.7)

   end function ic_use_initial_cp

!------------------------------------------------------------------------------------------------------------

   function ic_use_oblique(this)
!--function: determine from the ic-struct whether the computation uses an oblique view direction
      implicit none
!--result value
      logical                   :: ic_use_oblique
!--subroutine arguments
      class(t_ic),   intent(in) :: this

      ic_use_oblique = (this%discns1_eff.eq.8 .or. this%discns1_eff.eq.9)

   end function ic_use_oblique

!------------------------------------------------------------------------------------------------------------

   function ic_use_steep_slopes(this)
!--function: determine from the ic-struct whether the computation needs the steep slopes extension
      implicit none
!--result value
      logical                   :: ic_use_steep_slopes
!--subroutine arguments
      class(t_ic),   intent(in) :: this

      ic_use_steep_slopes = (this%discns1_eff.eq.6)

   end function ic_use_steep_slopes

!------------------------------------------------------------------------------------------------------------

   function ic_use_supergrid(this)
!--function: determine whether to snap the grid to rail profile s_r-coordinates
      implicit none
!--result value
      logical                     :: ic_use_supergrid
!--subroutine arguments
      class(t_ic),     intent(in) :: this

      ic_use_supergrid = (this%tang.eq.1 .or. this%tang.eq.2)

   end function ic_use_supergrid

!------------------------------------------------------------------------------------------------------------

   subroutine ic_unpack (modul, cpbtnfs, vldcmze, xhgiaowr, ic)
!--purpose: unpack the control-words into array ic
      implicit none
!--subroutine parameters:
      type(t_ic) :: ic
      integer    :: modul, cpbtnfs, vldcmze, xhgiaowr
!--local variables:
      integer, parameter :: idebug = 0
      integer            :: ihulp, fdigit, ddigit, zdigit, edigit

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
      fdigit        = (cpbtnfs - ihulp) / 10
         ihulp      = ihulp +      10 * fdigit
      ic%stress     = (cpbtnfs - ihulp) / 1
      if (idebug.ge.2) write(*,'(7(a,i2))') 'unpack: C=',ic%config,', P=',ic%pvtime,', B=',ic%bound,    &
         ', T=',ic%tang,', N=',ic%norm,', F=',fdigit,', S=',ic%stress

      if (modul.eq.1) then
         ic%force1 = fdigit
         ic%force3 = -99
      else
         ic%force1 = -99
         ic%force3 = fdigit
      endif

      ihulp         = 0
      ic%varfrc     = (vldcmze - ihulp) / 1000000
         ihulp      = ihulp + 1000000 * ic%varfrc
      ic%frclaw_inp = (vldcmze - ihulp) / 100000
         ihulp      = ihulp +  100000 * ic%frclaw_inp
      ddigit        = (vldcmze - ihulp) / 10000
         ihulp      = ihulp +   10000 * ddigit
      ic%gencr_inp  = (vldcmze - ihulp) / 1000
         ihulp      = ihulp +    1000 * ic%gencr_inp
      ic%mater      = (vldcmze - ihulp) / 100
         ihulp      = ihulp +     100 * ic%mater
      zdigit        = (vldcmze - ihulp) / 10
         ihulp      = ihulp +      10 * zdigit
      edigit        = (vldcmze - ihulp) / 1
      if (idebug.ge.2) write(*,'(7(a,i2))') 'unpack: V=',ic%varfrc,', L=',ic%frclaw_inp,                &
         ', D=',ddigit, ', C=',ic%gencr_inp,', M=',ic%mater, ', Z=',zdigit,', E=',edigit

      if (modul.eq.1) then
         ic%discns1_inp = ddigit
         ic%discns3 = 0
         ic%ztrack  = zdigit
         ic%ewheel  = edigit
         ic%rztang  = 0
         ic%rznorm  = 0
      else
         ic%discns1_inp = 0
         ic%discns3 = ddigit
         ic%ztrack  = 0
         ic%ewheel  = 0
         ic%rznorm  = zdigit
         ic%rztang  = edigit
      endif

      ihulp          = 0
      ic%xflow       = (xhgiaowr  - ihulp) / 10000000
         ihulp       = ihulp + 10000000 * ic%xflow
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

      if (idebug.ge.2) write(*,'(8(a,i2))') 'unpack: X=',ic%xflow,', H=',ic%heat, ', G=',ic%gausei_inp,  &
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

      if (modul.eq.1) then
         cpbtnfs =                                                    1000000*ic%config      +          &
                    100000*ic%pvtime      +    10000*ic%bound       +    1000*ic%tang        +          &
                       100*ic%norm        +       10*ic%force1      +       1*ic%stress
      else
         cpbtnfs =  100000*ic%pvtime      +    10000*ic%bound       +    1000*ic%tang        +          &
                       100*ic%norm        +       10*ic%force3      +       1*ic%stress
      endif

      if (modul.eq.1) then
         vldcmze =                                                    1000000*ic%varfrc      +          &
                    100000*ic%frclaw_inp  +    10000*ic%discns1_inp +    1000*ic%gencr_inp   +          &
                       100*ic%mater       +       10*ic%ztrack      +       1*ic%ewheel

      else
         vldcmze =  100000*ic%frclaw_inp  +    10000*ic%discns3     +    1000*ic%gencr_inp   +          &
                       100*ic%mater       +       10*ic%rznorm      +       1*ic%rztang
      endif

      xhgiaowr =                            10000000*ic%xflow       + 1000000*ic%heat        +          &
                    100000*ic%gausei_inp  +    10000*ic%iestim      +    1000*ic%matfil_surf +          &
                       100*ic%output_surf +       10*ic%flow        +       1*ic%return

   end subroutine ic_pack

!------------------------------------------------------------------------------------------------------------

   subroutine ic_unpack_dbg (modul, psflcin, ic)
!--purpose: unpack the debug flags into array ic
      implicit none
!--subroutine parameters:
      type(t_ic) :: ic
      integer    :: modul, psflcin
!--local variables:
      integer, parameter :: idebug = 0
      integer            :: ihulp

      ihulp         = 0
      ic%x_profil   = (psflcin - ihulp) / 1000000
         ihulp      = ihulp + 1000000 * ic%x_profil
      ic%x_smooth   = (psflcin - ihulp) / 100000
         ihulp      = ihulp +  100000 * ic%x_smooth
      ic%x_force    = (psflcin - ihulp) / 10000
         ihulp      = ihulp +   10000 * ic%x_force
      ic%x_locate   = (psflcin - ihulp) / 1000
         ihulp      = ihulp +    1000 * ic%x_locate
      ic%x_cpatch   = (psflcin - ihulp) / 100
         ihulp      = ihulp +     100 * ic%x_cpatch
      ic%x_inflcf   = (psflcin - ihulp) / 10
         ihulp      = ihulp +      10 * ic%x_inflcf
      ic%x_nmdbg    = (psflcin - ihulp) / 1

      if (idebug.ge.2) then
         write(bufout,'(7(a,i2))') ' dbg_unpack: P=',ic%x_profil,', S=',ic%x_smooth,                    &
             ', F=',ic%x_force,', L=',ic%x_locate,', C=',ic%x_cpatch,', I=',ic%x_inflcf,                &
             ', N=',ic%x_nmdbg
         call write_log(1, bufout)
      endif

      if (modul.ne.1) then
         ic%x_profil = 0
         ic%x_smooth = 0
         ic%x_force  = 0
         ic%x_locate = 0
      endif

   end subroutine ic_unpack_dbg

!------------------------------------------------------------------------------------------------------------

   subroutine ic_pack_dbg(psflcin, ic)
!--purpose: pack the debug flags in ic together in a control-word
      implicit none
!--subroutine parameters:
      type(t_ic) :: ic
      integer    :: psflcin

      psflcin =                                                    1000000*ic%x_profil    +             &
                 100000*ic%x_smooth    +    10000*ic%x_force     +    1000*ic%x_locate    +             &
                    100*ic%x_cpatch    +       10*ic%x_inflcf    +       1*ic%x_nmdbg

   end subroutine ic_pack_dbg

!------------------------------------------------------------------------------------------------------------

   subroutine ic_destroy( ic )
!--purpose: Destroy control digits 
      implicit none
!--subroutine arguments:
      type(t_ic) :: ic

      if (.false.) ic%config = 0    ! plain structure, no action needed

   end subroutine ic_destroy

!------------------------------------------------------------------------------------------------------------

   subroutine meta_init( m, ilevel )
!--purpose: Initialize meta-data to sensible default values
   implicit none
!--subroutine arguments:
   type(t_metadata) :: m
   integer          :: ilevel        ! 1 = time-varying parts; 2 = full init

   if (ilevel.ge.2) then
      m%REid       = 0
      m%CPid       = 0
      m%whl_ver    = 0
      m%rail_ver   = 0
      m%actv_thrd  = -1
      m%irun       = 0
      m%iax        = 0
      m%iside      = 0
      m%ncase      = 0
      m%itforce    = 0
      m%itforc_out = 0
      m%itforc_inn = 0
      m%npatch     = 0
      m%ipatch     = 0
      m%tim        = 0d0
      m%wrkdir     = ' '
      m%outdir     = ' '
      m%expnam     = ' '
   endif

   if (ilevel.ge.1) then
      m%s_ws      = 0d0
      m%th_ws     = 0d0
      m%ynom_whl  = 0d0
      m%rnom_whl  = 0d0
      m%rnom_rol  = 0d0

      m%x_w       = 0d0
      m%y_w       = 0d0
      m%z_w       = 0d0
      m%roll_w    = 0d0
      m%yaw_w     = 0d0

      m%y_r       = 0d0
      m%z_r       = 0d0
      m%roll_r    = 0d0

      m%xcp_tr    = 0d0
      m%ycp_tr    = 0d0
      m%zcp_tr    = 0d0
      m%deltcp_tr = 0d0

      m%xcp_r     = 0d0
      m%ycp_r     = 0d0
      m%zcp_r     = 0d0
      m%scp_r     = 0d0
      m%deltcp_r  = 0d0

      m%xcp_w     = 0d0
      m%ycp_w     = 0d0
      m%zcp_w     = 0d0
      m%scp_w     = 0d0
      m%deltcp_w  = 0d0
   endif

   end subroutine meta_init

!------------------------------------------------------------------------------------------------------------

   subroutine meta_copy( m_in, m_out )
!--purpose: Copy meta-data data-structure
   implicit none
!--subroutine arguments:
   type(t_metadata) :: m_in, m_out

   m_out = m_in                     ! plain structure, intrinsic copying

   end subroutine meta_copy

!------------------------------------------------------------------------------------------------------------

   subroutine meta_destroy( meta )
!--purpose: Destroy meta-data data-structure
      implicit none
!--subroutine arguments:
      type(t_metadata) :: meta

      if (.false.) meta%REid = 0    ! plain structure, no action needed

   end subroutine meta_destroy

!------------------------------------------------------------------------------------------------------------

   subroutine scaling_copy( s_in, s_out )
!--purpose: Copy scaling parameter data-structure
   implicit none
!--subroutine arguments:
   type(t_scaling) :: s_in, s_out

   s_out = s_in                     ! plain structure, intrinsic copying

   end subroutine scaling_copy

!------------------------------------------------------------------------------------------------------------

   subroutine scaling_destroy( scl )
!--purpose: Destroy scaling parameter data-structure
   implicit none
!--subroutine arguments:
   type(t_scaling) :: scl

      if (.false.) scl%len = 1d0    ! plain structure, no action needed

   end subroutine scaling_destroy

!------------------------------------------------------------------------------------------------------------

   subroutine mater_init( m, ic )
!--purpose: Initialize material-data to sensible default values
   implicit none
!--subroutine arguments:
   type(t_material) :: m
   type(t_ic)       :: ic

   m%gg      = (/    1.00d0,   1.00d0 /)
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
   m%gg3     = 1d0
   m%laythk  = 0d0
   m%tau_c0  = 0d0
   m%k_tau   = 0d0
   m%cdampn  = 0d0
   m%cdampt  = 0d0
   m%dfnmax  = 0d0
   m%dftmax  = 0d0
   m%bound_eff = ic%bound
   m%mater_eff = ic%mater
   m%gencr_eff = ic%gencr_inp
   m%fname_influe  = ' '
   m%if_meth = 0
   m%if_ver  = 4
   m%ninclin = 0

   end subroutine mater_init

!------------------------------------------------------------------------------------------------------------

   subroutine mater_copy( m_in, m_out )
!--purpose: copy material-data 
   implicit none
!--subroutine arguments:
   type(t_material) :: m_in, m_out

   m_out%gg                 = m_in%gg              
   m_out%poiss              = m_in%poiss   
   m_out%ak                 = m_in%ak      
   m_out%ga                 = m_in%ga      
   m_out%nu                 = m_in%nu
   m_out%bktemp             = m_in%bktemp  
   m_out%heatcp             = m_in%heatcp  
   m_out%lambda             = m_in%lambda  
   m_out%dens               = m_in%dens    
   m_out%betapl             = m_in%betapl  
   m_out%fg                 = m_in%fg
   m_out%tc                 = m_in%tc
   m_out%vt                 = m_in%vt
   m_out%akv                = m_in%akv
   m_out%gav                = m_in%gav
   m_out%nuv                = m_in%nuv
   m_out%flx_z              = m_in%flx_z   
   m_out%flx                = m_in%flx     
   m_out%k0_mf              = m_in%k0_mf   
   m_out%alfamf             = m_in%alfamf  
   m_out%betamf             = m_in%betamf  
   m_out%k_eff              = m_in%k_eff
   m_out%gg3                = m_in%gg3     
   m_out%laythk             = m_in%laythk  
   m_out%tau_c0             = m_in%tau_c0  
   m_out%k_tau              = m_in%k_tau   
   m_out%cdampn             = m_in%cdampn  
   m_out%cdampt             = m_in%cdampt  
   m_out%dfnmax             = m_in%dfnmax  
   m_out%dftmax             = m_in%dftmax  
   m_out%bound_eff          = m_in%bound_eff 
   m_out%mater_eff          = m_in%mater_eff 
   m_out%gencr_eff          = m_in%gencr_eff 
   m_out%if_meth            = m_in%if_meth 
   m_out%if_ver             = m_in%if_ver  
   m_out%ninclin            = m_in%ninclin 
   if (allocated(m_in%surf_inclin)) then
      m_out%surf_inclin     = m_in%surf_inclin    ! Fortran2003: automatic allocation
   endif
   m_out%fname_influe       = m_in%fname_influe  

   end subroutine mater_copy

!------------------------------------------------------------------------------------------------------------

   subroutine mater_destroy( mater )
!--purpose: Destroy material parameter data-structure
   implicit none
!--subroutine arguments:
   type(t_material) :: mater

   if (allocated(mater%surf_inclin)) deallocate(mater%surf_inclin)

   end subroutine mater_destroy

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

   subroutine subsurf_init(subs)
!--purpose: initialize data structure for subsurface stress calculation
      implicit none
!--subroutine arguments:
      type(t_subsurf)  :: subs

      subs%nblock = 0

   end subroutine subsurf_init

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

      call subsurf_destroy(s_out)

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

   subroutine subsurf_destroy(subs)
!--purpose: initialize data structure for subsurface stress calculation
      implicit none
!--subroutine arguments:
      type(t_subsurf)  :: subs
!--local variables:
      integer          :: iblk

      do iblk = 1, subs%nblock
         call subsblk_destroy(subs%blocks(iblk))
      enddo
      subs%nblock = 0

   end subroutine subsurf_destroy

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

   subroutine solv_copy ( solv_in, solv_out )
!--purpose: Copy solver setting data-structure
      implicit none
!--subroutine arguments:
      type(t_solvers) :: solv_in, solv_out

      solv_out = solv_in                ! plain structure, intrinsic copying

   end subroutine solv_copy

!------------------------------------------------------------------------------------------------------------

   subroutine solv_destroy( solv )
!--purpose: Destroy solver setting data-structure
   implicit none
!--subroutine arguments:
   type(t_solvers) :: solv

      if (.false.) solv%maxgs = 999    ! plain structure, no action needed

   end subroutine solv_destroy

!------------------------------------------------------------------------------------------------------------

   subroutine potcon_init ( potcon )
!--purpose: Set appropriate initial values for the potential contact area
      implicit none
!--subroutine arguments:
      type(t_potcon) :: potcon

      potcon%ipotcn =  1
      potcon%mx     =  1
      potcon%my     =  1
      potcon%dx     =  2d0
      potcon%dy     =  2d0
      potcon%xl     = -1d0
      potcon%yl     = -1d0

      ! compute other parameters like npot, dxdy, xc1, and so on.

      call potcon_fill( potcon )

   end subroutine potcon_init

!------------------------------------------------------------------------------------------------------------

   subroutine potcon_fill ( p, ipot_arg )
!--purpose: Complete data for the potential contact area based on non-Hertzian inputs
      implicit none
!--subroutine arguments:
      type(t_potcon)          :: p
      integer,       optional :: ipot_arg
!--local variables
      integer                 :: mypotcn

      mypotcn = p%ipotcn
      if (present(ipot_arg)) mypotcn = ipot_arg

      if (mypotcn.lt.1 .or. mypotcn.gt.4) then
         call write_log('Internal error(potcon_fill): invalid ipotcn')
         call abort_run()
      endif

      ! ipotcn == 1: the variables xl,  yl,  dx,  dy,  mx and my are given
      ! ipotcn == 2: the variables xl,  yl,  xh,  yh,  mx and my are given
      ! ipotcn == 3: the variables xc1, yc1, dx,  dy,  mx and my are given
      ! ipotcn == 4: the variables xc1, yc1, xcm, ycm, mx and my are given

      ! - fill in dx, dy if not specified

      if (mypotcn.eq.2) then
         p%dx = (p%xh - p%xl) / p%mx
         p%dy = (p%yh - p%yl) / p%my
      elseif (mypotcn.eq.4) then
         p%dx = (p%xcm - p%xc1) / max(1,p%mx-1)
         p%dy = (p%ycm - p%yc1) / max(1,p%my-1)
      endif

      ! - fill in xl, yl or xc1, yc1 if not specified

      if (mypotcn.eq.3 .or. mypotcn.eq.4) then
         p%xl = p%xc1 - 0.5d0 * p%dx
         p%yl = p%yc1 - 0.5d0 * p%dy
      else
         p%xc1 = p%xl + 0.5d0 * p%dx
         p%yc1 = p%yl + 0.5d0 * p%dy
      endif

      ! - fill in xh, yh and/or xcm, ycm if not specified

      if (mypotcn.eq.1 .or. mypotcn.eq.3 .or. mypotcn.eq.4) then
         p%xh = p%xl + p%mx * p%dx
         p%yh = p%yl + p%my * p%dy
      endif
      if (mypotcn.ge.1 .and. mypotcn.le.3) then
         p%xcm = p%xc1 + (p%mx - 1) * p%dx
         p%ycm = p%yc1 + (p%my - 1) * p%dy
      endif

      ! - fill in npot and dxdy

      p%npot = p%mx * p%my
      p%dxdy = p%dx * p%dy

   end subroutine potcon_fill

!------------------------------------------------------------------------------------------------------------

   subroutine potcon_copy ( potcon_in, potcon_out )
!--purpose: Copy potential contact area description
      implicit none
!--subroutine arguments:
      type(t_potcon) :: potcon_in, potcon_out

      potcon_out = potcon_in            ! plain structure, intrinsic copying

   end subroutine potcon_copy

!------------------------------------------------------------------------------------------------------------

   subroutine potcon_cgrid ( p, cgrid )
!--purpose: determine contact grid based for given potential contact area description
      implicit none
!--subroutine arguments:
      type(t_potcon) :: p
      type(t_grid)   :: cgrid

      call grid_create_uniform(cgrid, nxarg=p%mx, x0arg=p%xc1, dxarg=p%dx,                              &
                                      nyarg=p%my, y0arg=p%yc1, dyarg=p%dy, zarg=0d0)

   end subroutine potcon_cgrid

!------------------------------------------------------------------------------------------------------------

   subroutine potcon_hertz(hz, pot)
!--purpose: Complete the potcon specification for a Hertzian problem (mx, my filled in)
      implicit none
!--subroutine arguments:
      type(t_hertz)  :: hz
      type(t_potcon) :: pot
!--local variables :
      integer mypotcn, n_int

      ! ipotcn = -6: the potential contact must be derived from the SDEC parameters

      if (pot%ipotcn.eq.-6) then

         hz%bb  = 0.5d0 * (hz%bneg + hz%bpos)
         pot%xl = -hz%scale * max(1d-9, hz%aa)
         pot%yl = -hz%scale * max(1d-9, hz%bb)
         pot%xh = -pot%xl
         pot%yh = -pot%yl
         mypotcn = 2

      ! ipotcn = -5..-4: potential contact derived from 2D Hertzian solution computed by hzsol.

      elseif (pot%ipotcn.le.-4) then

         pot%xl = -hz%scale * max(1d-9, hz%aa)
         pot%xh = -pot%xl
         pot%dx = (pot%xh - pot%xl) / pot%mx

         ! let my == (n_int + n_ext) \approx scale * n_int
         ! choose dy such that n_int * dy is precisely 2 * bb
         ! when scale<1 set n_ext = 0
         ! use even n_ext, avoid n_int<1 for large scale

         if (hz%scale.le.1d0) then
            n_int = pot%my
         else
            n_int = max(1, nint(pot%my / hz%scale))
            if (mod(pot%my-n_int,2).eq.1) n_int = n_int  + 1
         endif

         pot%dy = 2d0*max(1d-9,hz%bb) / n_int 
         pot%yl = -pot%my*pot%dy / 2d0
         pot%yh = -pot%yl
         mypotcn = 2

         if (.false.) then
            write(bufout,'(a,i6,a,f6.3,a,2i4)') ' my=',pot%my,', scale=',hz%scale,': n_int, n_ext=',    &
                        n_int, pot%my-n_int
            call write_log(1, bufout)
            write(bufout,'(a,f6.3,a,f7.1,a,f6.1)') ' dy=',pot%dy,', 2*bb=',2*hz%bb,', n_int*dy=',       &
                        n_int*pot%dy
            call write_log(1, bufout)
         endif

      ! ipotcn = -3..-1: potential contact derived from 3D Hertzian solution computed by hzsol.

      elseif (pot%ipotcn.le.-1) then

         pot%xl = -hz%scale * max(1d-9, hz%aa)
         pot%yl = -hz%scale * max(1d-9, hz%bb)
         pot%xh = -pot%xl
         pot%yh = -pot%yl
         mypotcn = 2

      endif

      ! complete entries in potcon description

      call potcon_fill ( pot, mypotcn )

   end subroutine potcon_hertz

!------------------------------------------------------------------------------------------------------------

subroutine potcon_get_overlap(p_old, p_new, kofs_x, kofs_y, ix0, ix1, iy0, iy1, is_ok, is_equal)
!--purpose: check if potential contact grids are matching (subset of same super-grid) and 
!           determine overlap [ix0:ix1] x [iy0:iy1] within p_new
   implicit none
!--subroutine arguments:
   type(t_potcon)            :: p_old, p_new
   integer,     intent(out)  :: kofs_x, kofs_y, ix0, ix1, iy0, iy1
   logical,     intent(out)  :: is_ok, is_equal
!--local variables:

   call grid_overlap(p_old%xl, p_old%dx, p_old%mx, p_old%yl, p_old%dy, p_old%my,                        &
                     p_new%xl, p_new%dx, p_new%mx, p_new%yl, p_new%dy, p_new%my,                        &
                     kofs_x, kofs_y, ix0, ix1, iy0, iy1, is_ok, is_equal)

end subroutine potcon_get_overlap

!------------------------------------------------------------------------------------------------------------

   subroutine potcon_merge( pot1, pot2, pot_tot, idebug, ierror )
!--purpose: Create new potential contact encompassing two existing pot.contact areas
   implicit none
!--subroutine arguments:
   type(t_potcon) :: pot1, pot2, pot_tot
   integer        :: idebug, ierror
!--local variables:
   logical        :: is_ok, is_equal
   integer        :: ix0, ix1, iy0, iy1, kofs_x, kofs_y

   if (idebug.ge.2) then
      write(bufout,'(2(a,f8.3),2(a,i4),2(a,f8.3))') ' pot1: xl,yl=',pot1%xl,',',pot1%yl,                &
                ', mx,my=',pot1%mx,',',pot1%my,', xh,yh=',pot1%xh,',',pot1%yh
      call write_log(1, bufout)
      write(bufout,'(2(a,f8.3),2(a,i4),2(a,f8.3))') ' pot2: xl,yl=',pot2%xl,',',pot2%yl,                &
                ', mx,my=',pot2%mx,',',pot2%my,', xh,yh=',pot2%xh,',',pot2%yh
      call write_log(1, bufout)
   endif

   ! check that the two grids have matching (dx,dy) and determine offsets (kx,ky)

   call potcon_get_overlap(pot1, pot2, kofs_x, kofs_y, ix0, ix1, iy0, iy1, is_ok, is_equal)

   ierror = 0
   if (.not.is_ok) ierror = 1

   if (idebug.ge.2) then
      write(bufout,'(2(a,i3))') ' offset k_x=',kofs_x,', k_y=',kofs_y
      call write_log(1, bufout)
   endif

   ! determine extent of new grid. kx,ky == #extra columns/rows at start of pot1 wrt pot2

   if (kofs_x.ge.0) then
      pot_tot%xl = pot1%xl
      pot_tot%mx = max(pot1%mx, pot2%mx+kofs_x)
   else
      pot_tot%xl = pot2%xl
      pot_tot%mx = max(pot2%mx, pot1%mx-kofs_x)
   endif
   if (kofs_y.gt.0) then
      pot_tot%yl = pot1%yl
      pot_tot%my = max(pot1%my, pot2%my+kofs_y)
   else
      pot_tot%yl = pot2%yl
      pot_tot%my = max(pot2%my, pot1%my-kofs_y)
   endif
   pot_tot%ipotcn = 1
   pot_tot%dx  = pot2%dx
   pot_tot%dy  = pot2%dy

   ! complete pot.con specification

   call potcon_fill(pot_tot)

   if (idebug.ge.2) then
      write(bufout,'(2(a,f8.3),2(a,i4),2(a,f8.3))') '  pot: xl,yl=',pot_tot%xl,',',pot_tot%yl,          &
                ', mx,my=',pot_tot%mx,',',pot_tot%my,', xh,yh=',pot_tot%xh,',',pot_tot%yh
      call write_log(1, bufout)
   endif

   end subroutine potcon_merge

!------------------------------------------------------------------------------------------------------------

   subroutine potcon_destroy( pot )
!--purpose: Destroy potential contact data-structure
   implicit none
!--subroutine arguments:
   type(t_potcon) :: pot

      if (.false.) pot%mx = 1    ! plain structure, no action needed

   end subroutine potcon_destroy

!------------------------------------------------------------------------------------------------------------

   subroutine geom_init ( geom, cgrid )
!--purpose: Set appropriate initial values for the geometry description
      implicit none
!--subroutine arguments:
      type(t_geomet) :: geom
      type(t_grid)   :: cgrid
!--local variables:
      integer j

      geom%ibase  = 1
      geom%iplan  = 1
      geom%nn     = 0
      geom%npatch = 0
      geom%prmudf = (/ 1d0, 0d0, 1d0, (0d0, j=4,10) /)  ! automatic allocation
      geom%prmpln = (/ (0d0, j=1,10) /)                 ! automatic allocation
      ! xylim, facsep

      call gf3_new( geom%exrhs, 'geom%exrhs', cgrid, nulify=.true. )
      call gf3_new( geom%hs1, 'geom%hs1', cgrid, nulify=.true. )
      call gf3_new( geom%hv1, 'geom%hv1', cgrid, nulify=.true. )

   end subroutine geom_init

!------------------------------------------------------------------------------------------------------------

   subroutine geom_copy ( geom_in, geom_out, cgrid )
!--purpose: Copy a geometry description
      implicit none
!--subroutine arguments:
      type(t_geomet) :: geom_in, geom_out
      type(t_grid)   :: cgrid

      geom_out%ibase  = geom_in%ibase
      geom_out%iplan  = geom_in%ibase
      geom_out%nn     = geom_in%nn   
      geom_out%npatch = geom_in%npatch

      geom_out%prmudf = geom_in%prmudf  ! automatic (re)allocation
      geom_out%prmpln = geom_in%prmpln  ! automatic (re)allocation
      if (allocated(geom_in%xylim))  geom_out%xylim  = geom_in%xylim   ! automatic (re)allocation
      if (allocated(geom_in%facsep)) geom_out%facsep = geom_in%facsep  ! automatic (re)allocation

      call gf3_new( geom_out%exrhs, 'geom%exrhs', cgrid, nulify=.true. )
      call gf3_new( geom_out%hs1,   'geom%hs1'  , cgrid, nulify=.true. )
      call gf3_new( geom_out%hv1,   'geom%hv1'  , cgrid, nulify=.true. )

      call gf3_copy( AllElm, geom_in%exrhs, geom_out%exrhs, ikALL )
      call gf3_copy( AllElm, geom_in%hs1,   geom_out%hs1,   ikALL )
      call gf3_copy( AllElm, geom_in%hv1,   geom_out%hv1,   ikALL )

   end subroutine geom_copy

!------------------------------------------------------------------------------------------------------------

   subroutine geom_destroy ( geom )
!--purpose: Destroy a geometry description
      implicit none
!--subroutine arguments:
      type(t_geomet) :: geom

      if (allocated(geom%prmudf)) deallocate(geom%prmudf)
      if (allocated(geom%prmpln)) deallocate(geom%prmpln)
      if (allocated(geom%xylim))  deallocate(geom%xylim)
      if (allocated(geom%facsep)) deallocate(geom%facsep)
      call gf3_destroy(geom%exrhs)
      call gf3_destroy(geom%hs1)
      call gf3_destroy(geom%hv1)

   end subroutine geom_destroy

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
      kin%veloc   =  1d3
      kin%dq      =  dx
      kin%facphi  =  1d0 / 6d0
      kin%chi     =  0d0
      kin%fntrue  =  0d0
      kin%pen     =  0d0
      kin%penv    =  0d0
      kin%fxrel   =  0d0
      kin%fyrel   =  0d0
      kin%cksi    =  0d0
      kin%ceta    =  0d0
      kin%cphi    =  0d0
      kin%spinxo  =  0d0
      kin%spinyo  =  0d0
      kin%fcntc   = (/ 0d0, 0d0, 0d0 /)
      kin%fprev   = (/ 0d0, 0d0, 0d0 /)
      kin%fdamp   = (/ 0d0, 0d0, 0d0 /)
      kin%use_muscal = (ic%varfrc.eq.0 .or. ic%varfrc.eq.2)
      kin%muscal  =  1d0
      if (kin%use_muscal) kin%muscal = fric%fstat()

   end subroutine kincns_init

!------------------------------------------------------------------------------------------------------------

   subroutine kincns_copy ( kin_in, kin_out )
!--purpose: Copy the kinematic data
      implicit none
!--subroutine arguments:
      type(t_kincns)  :: kin_in, kin_out

      kin_out = kin_in             ! plain structure, intrinsic copying

   end subroutine kincns_copy

!------------------------------------------------------------------------------------------------------------

   subroutine kincns_destroy( kin )
!--purpose: Destroy kinematic data data-structure
   implicit none
!--subroutine arguments:
   type(t_kincns) :: kin

      if (.false.) kin%dq = 1d0   ! plain structure, no action needed

   end subroutine kincns_destroy

!------------------------------------------------------------------------------------------------------------

   subroutine hertz_init ( hertz )
!--purpose: Set appropriate initial values for the Hertzian data
      implicit none
!--subroutine arguments:
      type(t_hertz) :: hertz

      hertz%aa    =  5d0
      hertz%bb    =  5d0
      hertz%bneg  =  hertz%bb
      hertz%bpos  =  hertz%bb
      hertz%a1    =  0.001d0   ! Rx=500mm
      hertz%b1    =  0.001d0
      hertz%aob   =  1d0
      hertz%scale =  1d0

   end subroutine hertz_init

!------------------------------------------------------------------------------------------------------------

   subroutine hertz_copy ( hz_in, hz_out )
!--purpose: Copy Hertzian data data-structure
      implicit none
!--subroutine arguments:
      type(t_hertz) :: hz_in, hz_out

      hz_out = hz_in               ! plain structure, intrinsic copying

   end subroutine hertz_copy

!------------------------------------------------------------------------------------------------------------

   subroutine hertz_destroy( hz )
!--purpose: Destroy Hertzian data data-structure
   implicit none
!--subroutine arguments:
   type(t_hertz) :: hz

      if (.false.) hz%aa = 1d0    ! plain structure, no action needed

   end subroutine hertz_destroy

!------------------------------------------------------------------------------------------------------------

   subroutine output_init ( outp, cgrid )
!--purpose: Initialize output-data data-structure
      implicit none
!--subroutine arguments:
      type(t_output) :: outp
      type(t_grid)   :: cgrid

      call eldiv_new( outp%igs, cgrid, nulify=.true. )
      call eldiv_new( outp%igv, cgrid, nulify=.true. )
      call eldiv_exter(outp%igs)

      call gf3_new(outp%mus,   'outpt1%mus',   cgrid, nulify=.true.)
      call gf3_new(outp%shft,  'outpt1%shft',  cgrid, nulify=.true.)
      call gf3_new(outp%ps,    'outpt1%ps',    cgrid, nulify=.true., lzero=.true.)
      call gf3_new(outp%us,    'outpt1%us',    cgrid, nulify=.true.)
      call gf3_new(outp%ss,    'outpt1%ss',    cgrid, nulify=.true.)
      call gf3_new(outp%taucs, 'outpt1%taucs', cgrid, nulify=.true.)
      call gf3_new(outp%upls,  'outpt1%upls',  cgrid, nulify=.true.)
      call gf3_new(outp%temp1, 'outpt1%temp1', cgrid, nulify=.true.)
      call gf3_new(outp%temp2, 'outpt1%temp2', cgrid, nulify=.true.)

      call gf3_new(outp%muv,   'outpt1%muv',   cgrid, nulify=.true.)
      call gf3_new(outp%pv,    'outpt1%pv',    cgrid, nulify=.true.)
      call gf3_new(outp%uv,    'outpt1%uv',    cgrid, nulify=.true.)
      call gf3_new(outp%sv,    'outpt1%sv',    cgrid, nulify=.true.)
      call gf3_new(outp%taucv, 'outpt1%taucv', cgrid, nulify=.true.)
      call gf3_new(outp%uplv,  'outpt1%uplv',  cgrid, nulify=.true.)

   end subroutine output_init

!------------------------------------------------------------------------------------------------------------

   subroutine output_copy ( outp_in, outp_out, cgrid )
!--purpose: Copy output-data data-structure
      implicit none
!--subroutine arguments:
      type(t_output)         :: outp_in, outp_out
      type(t_grid),  target  :: cgrid

      outp_out%sens   = outp_in%sens
      outp_out%mxtrue = outp_in%mxtrue
      outp_out%mytrue = outp_in%mytrue
      outp_out%mztrue = outp_in%mztrue
      outp_out%elen   = outp_in%elen  
      outp_out%frpow  = outp_in%frpow 
      outp_out%pmax   = outp_in%pmax  

      ! create new element division linking to new grid 'cgrid', copy values

      call eldiv_new( outp_out%igs, cgrid, nulify=.true. )
      call eldiv_new( outp_out%igv, cgrid, nulify=.true. )
      call eldiv_copy( outp_in%igs, outp_out%igs, ikALL )
      call eldiv_copy( outp_in%igv, outp_out%igv, ikALL )

      ! create new grid functions linking to new grid 'cgrid', new igs/igv, copy values

      call gf3_copy_full( outp_in%mus,   outp_out%mus,   cgrid, outp_out%igs )
      call gf3_copy_full( outp_in%shft,  outp_out%shft,  cgrid, outp_out%igs )
      call gf3_copy_full( outp_in%ps,    outp_out%ps,    cgrid, outp_out%igs )
      call gf3_copy_full( outp_in%us,    outp_out%us,    cgrid, outp_out%igs )
      call gf3_copy_full( outp_in%ss,    outp_out%ss,    cgrid, outp_out%igs )
      call gf3_copy_full( outp_in%taucs, outp_out%taucs, cgrid, outp_out%igs )
      call gf3_copy_full( outp_in%upls,  outp_out%upls,  cgrid, outp_out%igs )
      call gf3_copy_full( outp_in%temp1, outp_out%temp1, cgrid, outp_out%igs )
      call gf3_copy_full( outp_in%temp2, outp_out%temp2, cgrid, outp_out%igs )

      call gf3_copy_full( outp_in%muv,   outp_out%muv,   cgrid, outp_out%igv )
      call gf3_copy_full( outp_in%pv,    outp_out%pv,    cgrid, outp_out%igv )
      call gf3_copy_full( outp_in%uv,    outp_out%uv,    cgrid, outp_out%igv )
      call gf3_copy_full( outp_in%sv,    outp_out%sv,    cgrid, outp_out%igv )
      call gf3_copy_full( outp_in%taucv, outp_out%taucv, cgrid, outp_out%igv )
      call gf3_copy_full( outp_in%uplv,  outp_out%uplv,  cgrid, outp_out%igv )

   end subroutine output_copy

!------------------------------------------------------------------------------------------------------------

   subroutine output_destroy ( outp )
!--purpose: Destroy output-data data-structure
      implicit none
!--subroutine arguments:
      type(t_output) :: outp

      call eldiv_destroy(outp%igs)
      call eldiv_destroy(outp%igv)

      call gf3_destroy( outp%mus )
      call gf3_destroy( outp%shft )
      call gf3_destroy( outp%ps )
      call gf3_destroy( outp%us )
      call gf3_destroy( outp%ss )
      call gf3_destroy( outp%taucs )
      call gf3_destroy( outp%upls )
      call gf3_destroy( outp%temp1 )
      call gf3_destroy( outp%temp2 )

      call gf3_destroy( outp%muv )
      call gf3_destroy( outp%pv )
      call gf3_destroy( outp%uv )
      call gf3_destroy( outp%sv )
      call gf3_destroy( outp%taucv )
      call gf3_destroy( outp%uplv )

   end subroutine output_destroy

!------------------------------------------------------------------------------------------------------------

   subroutine gd_init ( gd )
!--purpose: Set appropriate initial values for the control/input-variables
      implicit none
!--subroutine arguments:
      type(t_probdata) :: gd

      call meta_init    ( gd%meta, 2 )
    ! call scaling_init ( gd%scl )
      call ic_init      ( gd%ic )
      call mater_init   ( gd%mater, gd%ic )
      call potcon_init  ( gd%potcon_inp )
      call potcon_init  ( gd%potcon_cur )
      call potcon_cgrid ( gd%potcon_inp, gd%cgrid_inp )
      call potcon_cgrid ( gd%potcon_cur, gd%cgrid_cur )
      call hertz_init   ( gd%hertz )
      call geom_init    ( gd%geom, gd%cgrid_cur )
      call fric_init    ( gd%fric )
      call kincns_init  ( gd%kin, gd%ic, gd%fric, gd%cgrid_cur%dx )
    ! call influe_init  ( gd%influ )
      call solv_init    ( gd%solv )
      call output_init  ( gd%outpt1, gd%cgrid_cur )
      call subsurf_init ( gd%subs )

   end subroutine gd_init

!------------------------------------------------------------------------------------------------------------

   subroutine gd_copy ( gd_in, gd_out )
!--purpose: Full copy of a gd data-structure
      implicit none
!--subroutine arguments:
      type(t_probdata) :: gd_in, gd_out

      call meta_copy    ( gd_in%meta,   gd_out%meta )
      call scaling_copy ( gd_in%scl,    gd_out%scl )
      call ic_copy      ( gd_in%ic,     gd_out%ic )
      call mater_copy   ( gd_in%mater,  gd_out%mater )
      call potcon_copy  ( gd_in%potcon_inp,  gd_out%potcon_inp )
      call potcon_copy  ( gd_in%potcon_cur,  gd_out%potcon_cur )
      call potcon_cgrid ( gd_out%potcon_inp, gd_out%cgrid_inp )
      call potcon_cgrid ( gd_out%potcon_cur, gd_out%cgrid_cur )
      call hertz_copy   ( gd_in%hertz,  gd_out%hertz )
      call geom_copy    ( gd_in%geom,   gd_out%geom, gd_out%cgrid_cur )
      call fric_copy    ( gd_in%fric,   gd_out%fric )
      call kincns_copy  ( gd_in%kin,    gd_out%kin )
      call influe_copy  ( gd_in%influ,  gd_out%influ )
      call solv_copy    ( gd_in%solv,   gd_out%solv )
      call output_copy  ( gd_in%outpt1, gd_out%outpt1, gd_out%cgrid_cur )
      call subsurf_copy ( gd_in%subs,   gd_out%subs )

   end subroutine gd_copy

!------------------------------------------------------------------------------------------------------------

   subroutine gd_resize_gridfunc(cgrid, cgrid_new, geom, outpt, keep_grid_ptr)
!--purpose: resize grid-functions used in hierarchical data-structure, shifting data
      implicit none
!--subroutine arguments:
      type(t_grid),  target   :: cgrid, cgrid_new
      type(t_geomet)          :: geom
      type(t_output)          :: outpt
      logical,       optional :: keep_grid_ptr
!--local variables:
      logical                 :: keep_ptr

      keep_ptr = .false.
      if (present(keep_grid_ptr)) keep_ptr = keep_grid_ptr

      ! resize arrays for element division

      call eldiv_resize(outpt%igs, cgrid_new)
      call eldiv_resize(outpt%igv, cgrid_new)
      if (keep_ptr) outpt%igs%grid => cgrid
      if (keep_ptr) outpt%igv%grid => cgrid

      ! resize hs1 for right-hand side

      call gf3_resize(geom%hs1,   cgrid_new, 1d0)
      call gf3_resize(geom%hv1,   cgrid_new, 1d0)
      call gf3_resize(geom%exrhs, cgrid_new, 1d0)
      if (keep_ptr) geom%hs1%grid   => cgrid
      if (keep_ptr) geom%hv1%grid   => cgrid
      if (keep_ptr) geom%exrhs%grid => cgrid

      ! resize arrays for friction coefficients, tractions, displacements and shift

      call gf3_resize(outpt%mus, cgrid_new)
      call gf3_resize(outpt%ps,  cgrid_new, 0d0)
      call gf3_resize(outpt%us,  cgrid_new)
      call gf3_resize(outpt%ss,  cgrid_new)
      if (keep_ptr) outpt%mus%grid => cgrid
      if (keep_ptr) outpt%ps%grid  => cgrid
      if (keep_ptr) outpt%us%grid  => cgrid
      if (keep_ptr) outpt%ss%grid  => cgrid

      call gf3_resize(outpt%muv, cgrid_new)
      call gf3_resize(outpt%pv,  cgrid_new, 0d0)
      call gf3_resize(outpt%uv,  cgrid_new)
      call gf3_resize(outpt%sv,  cgrid_new)
      if (keep_ptr) outpt%muv%grid => cgrid
      if (keep_ptr) outpt%pv%grid  => cgrid
      if (keep_ptr) outpt%uv%grid  => cgrid
      if (keep_ptr) outpt%sv%grid  => cgrid

      ! resize arrays for plastic deformation, yield point, temperature

      call gf3_resize(outpt%taucs, cgrid_new)
      call gf3_resize(outpt%upls,  cgrid_new)
      call gf3_resize(outpt%taucv, cgrid_new)
      call gf3_resize(outpt%uplv,  cgrid_new)
      call gf3_resize(outpt%temp1, cgrid_new)
      call gf3_resize(outpt%temp2, cgrid_new)
      if (keep_ptr) outpt%taucs%grid => cgrid
      if (keep_ptr) outpt%upls%grid  => cgrid
      if (keep_ptr) outpt%taucv%grid => cgrid
      if (keep_ptr) outpt%uplv%grid  => cgrid
      if (keep_ptr) outpt%temp1%grid => cgrid
      if (keep_ptr) outpt%temp2%grid => cgrid

   end subroutine gd_resize_gridfunc

!------------------------------------------------------------------------------------------------------------

   subroutine gd_merge_gridfunc(geom_tot, geom_add, outpt_tot, outpt_add, idebug)
!--purpose: merge grid-functions from gd data-structure gd_add into gd data-structure gd_tot
!           requires matching grids; contact regions must not overlap in x- or y-direction.
      implicit none
!--subroutine arguments:
      type(t_geomet)          :: geom_tot, geom_add
      type(t_output)          :: outpt_tot, outpt_add
      integer                 :: idebug
!--local variables:
      integer, parameter :: npatch = 2
      integer        :: mx, ii, ipatch
      logical        :: is_ok
      real(kind=8)   :: xylim(npatch,4)
      type(t_eldiv)  :: mask

      associate(cgrid_add => outpt_add%ps%grid, cgrid_tot => outpt_tot%ps%grid,                         &
                ix0_add => outpt_add%igs%ixmin, ix1_add => outpt_add%igs%ixmax,                         &
                iy0_add => outpt_add%igs%iymin, iy1_add => outpt_add%igs%iymax,                         &
                ix0_tot => outpt_tot%igs%ixmin, ix1_tot => outpt_tot%igs%ixmax,                         &
                iy0_tot => outpt_tot%igs%iymin, iy1_tot => outpt_tot%igs%iymax)

      ! check that grids are equal (TODO: check xl,yl,dx,dy?)

      is_ok  = equal_grid_sizes(cgrid_tot%dx, cgrid_add%dx, cgrid_tot%dy, cgrid_add%dy)
      is_ok  = is_ok .and. abs(cgrid_tot%x(1)-cgrid_add%x(1)).lt.1d-3*cgrid_tot%dx
      is_ok  = is_ok .and. abs(cgrid_tot%y(1)-cgrid_add%y(1)).lt.1d-3*cgrid_tot%dy
      is_ok  = is_ok .and. cgrid_add%ntot.eq.cgrid_tot%ntot
      if (.not.is_ok) then
         call write_log(' Internal error(gd_merge_gridfunc): grids not equal.')
         call abort_run()
      endif

      ! create mask-array for sub-patches within pot.contact

      call eldiv_new(mask, cgrid_tot)
      call eldiv_exter(mask)

      ! determine extent of actual contact areas

      ipatch = 1
      mx     = cgrid_tot%nx
      call areas(outpt_tot%igs)
      xylim(ipatch,1) = cgrid_tot%x(ix0_tot) - 0.1d0 * cgrid_tot%dx
      xylim(ipatch,2) = cgrid_tot%x(ix1_tot) + 0.1d0 * cgrid_tot%dx
      xylim(ipatch,3) = cgrid_tot%y(iy0_tot*mx) - 0.1d0 * cgrid_tot%dy
      xylim(ipatch,4) = cgrid_tot%y(iy1_tot*mx) + 0.1d0 * cgrid_tot%dy

      ipatch = 2
      mx     = cgrid_add%nx
      call areas(outpt_add%igs)
      xylim(ipatch,1) = cgrid_add%x(ix0_add) - 0.1d0 * cgrid_add%dx
      xylim(ipatch,2) = cgrid_add%x(ix1_add) + 0.1d0 * cgrid_add%dx
      xylim(ipatch,3) = cgrid_add%y(iy0_add*mx) - 0.1d0 * cgrid_add%dy
      xylim(ipatch,4) = cgrid_add%y(iy1_add*mx) + 0.1d0 * cgrid_add%dy

      ! warn in case actual contact areas are overlapping
      ! ok:  [y0_add,y1_add] < [y0_tot,y1_tot]  or  [y0_add,y1_add] > [y0_tot,y1_tot]

      is_ok = (ix1_add.lt.ix0_tot .or. ix0_add.gt.ix1_tot .or. iy1_add.lt.iy0_tot .or. iy0_add.gt.iy1_tot)
      if (.not.is_ok) then
         call write_log(' ERROR(gd_merge_gridfunc): contact areas (boxes) overlap with each other.')
      endif

      if (idebug.ge.2 .or. .not.is_ok) then
         write(bufout,'(4(a,i3),a)') ' gd_tot: actual contact on ix=[',ix0_tot,',',ix1_tot,             &
                '], iy=[',iy0_tot,',',iy1_tot,']'
         call write_log(1, bufout)
         call wrigs (outpt_tot%igs, .true., 0d0, .false.)
         write(bufout,'(4(a,i3),a)') ' gd_add: actual contact on ix=[',ix0_add,',',ix1_add,             &
                '], iy=[',iy0_add,',',iy1_add,']'
         call write_log(1, bufout)
         call wrigs (outpt_add%igs, .true., 0d0, .false.)
      endif

      call eldiv_cpatches(mask, npatch, xylim, idebug)

      if (idebug.ge.2 .or. .not.is_ok) then
         call write_log(' mask array, *=tot, S=add:')
         call wrigs (mask, .true., 0d0, .false.)
      endif

      ! merge hs1, hv1 for right-hand side, adding contents of gd_add (ip=2 in mask) to gd_tot (ip=1)

      call gf3_msk_copy(2, mask, geom_add%hs1, geom_tot%hs1, ikALL)
      call gf3_msk_copy(2, mask, geom_add%hv1, geom_tot%hv1, ikALL)

      ! merge arrays for friction coefficients, tractions, displacements and shift

      call gf3_msk_copy(2, mask, outpt_add%mus, outpt_tot%mus, ikALL)
      call gf3_msk_copy(2, mask, outpt_add%ps,  outpt_tot%ps,  ikALL)
      call gf3_msk_copy(2, mask, outpt_add%us,  outpt_tot%us,  ikALL)
      call gf3_msk_copy(2, mask, outpt_add%ss,  outpt_tot%ss,  ikALL)

      call gf3_msk_copy(2, mask, outpt_add%muv, outpt_tot%muv, ikALL)
      call gf3_msk_copy(2, mask, outpt_add%pv,  outpt_tot%pv,  ikALL)
      call gf3_msk_copy(2, mask, outpt_add%uv,  outpt_tot%uv,  ikALL)
      call gf3_msk_copy(2, mask, outpt_add%sv,  outpt_tot%sv,  ikALL)

      ! merge arrays for plastic deformation, yield point, temperature

      call gf3_msk_copy(2, mask, outpt_add%taucs, outpt_tot%taucs, ikALL)
      call gf3_msk_copy(2, mask, outpt_add%upls,  outpt_tot%upls,  ikALL)
      call gf3_msk_copy(2, mask, outpt_add%taucv, outpt_tot%taucv, ikALL)
      call gf3_msk_copy(2, mask, outpt_add%uplv,  outpt_tot%uplv,  ikALL)
      call gf3_msk_copy(2, mask, outpt_add%temp1, outpt_tot%temp1, ikALL)
      call gf3_msk_copy(2, mask, outpt_add%temp2, outpt_tot%temp2, ikALL)

      ! merge arrays for element division

      do ii = 1, cgrid_tot%ntot
         outpt_tot%igs%el(ii) = max(outpt_tot%igs%el(ii), outpt_add%igs%el(ii))
         outpt_tot%igv%el(ii) = max(outpt_tot%igv%el(ii), outpt_add%igv%el(ii))
      enddo

      call eldiv_destroy(mask)

      end associate

   end subroutine gd_merge_gridfunc

!------------------------------------------------------------------------------------------------------------

   subroutine gd_merge ( gd_tot, gd_add, idebug )
!--purpose: merge gd data-structure gd_add into gd data-structure gd_tot
!           requires matching grids; contact regions must not overlap in y-direction.
!           gd_add must lie to left of gd_tot: y-range [y0_add,y1_add] < [y0_tot,y1_tot].
      implicit none
!--subroutine arguments:
      type(t_probdata) :: gd_tot, gd_add
      integer          :: idebug
!--local variables:
      type(t_potcon)   :: pot_new
      type(t_grid)     :: cgrid_new
      integer          :: ierror

      call potcon_merge( gd_tot%potcon_cur, gd_add%potcon_cur, pot_new, idebug, ierror )
      call potcon_cgrid(pot_new, cgrid_new)

      if (ierror.ne.0) then
         call write_log(' Internal error(gd_merge): grids cannot be merged')
         call abort_run()
      endif

      ! enlarge all grid-functions in gd_tot and gd_add to the size of cgrid_new

      call potcon_copy(pot_new, gd_tot%potcon_cur)
      call potcon_copy(pot_new, gd_add%potcon_cur)

      call gd_resize_gridfunc(gd_tot%cgrid_cur, cgrid_new, gd_tot%geom, gd_tot%outpt1, .true.)
      call gd_resize_gridfunc(gd_add%cgrid_cur, cgrid_new, gd_add%geom, gd_add%outpt1, .true.)

      call grid_copy(cgrid_new, gd_tot%cgrid_cur)
      call grid_copy(cgrid_new, gd_add%cgrid_cur)
      call grid_destroy(cgrid_new)

      ! merge grid-functions in geom, outpt using gd_add on left side into gd_tot on right side

      call gd_merge_gridfunc(gd_tot%geom, gd_add%geom, gd_tot%outpt1, gd_add%outpt1, idebug)

   end subroutine gd_merge

!------------------------------------------------------------------------------------------------------------

   subroutine gd_destroy(gd)
!--function: destroy the contact hierarchical data-structure gd
   implicit none
!--subroutine arguments:
   type(t_probdata) :: gd

   call meta_destroy    ( gd%meta )
   call scaling_destroy ( gd%scl )
   call ic_destroy      ( gd%ic )
   call mater_destroy   ( gd%mater )
   call potcon_destroy  ( gd%potcon_inp )
   call potcon_destroy  ( gd%potcon_cur )
   call grid_destroy    ( gd%cgrid_inp )
   call grid_destroy    ( gd%cgrid_cur )
   call hertz_destroy   ( gd%hertz )
   call geom_destroy    ( gd%geom )
   call fric_destroy    ( gd%fric )
   call kincns_destroy  ( gd%kin )
   call influe_destroy  ( gd%influ )
   call solv_destroy    ( gd%solv )
   call output_destroy  ( gd%outpt1 )
   call subsurf_destroy ( gd%subs )

   end subroutine gd_destroy

!------------------------------------------------------------------------------------------------------------

end module m_hierarch_data
