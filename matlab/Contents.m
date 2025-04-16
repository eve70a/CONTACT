% 
% Matlab scripts for visualisation of output of the CONTACT program (trunk).
%
% Loading results into Matlab.
%
%   loadcase   - load the results for one case for a given experiment name.
%   diffcase   - compute the difference of results for two cases
%   loadstrs   - load the results for a subsurface stress calculation.
%   diffstrs   - compute the difference of results for two cases w.r.t.
%                subsurface stress calculations.
%
%   parse_out1 - sample script to read output of the contact patches for
%                wheel-rail contact cases (module 1) from an .out-file.
%   parse_out3 - sample script to read the creepages and forces of generic
%                contact cases (module 3) from an .out-file.
%
% Visualizing the results.
%
%   plot2d     - 2D plots of tractions for rows or columns of the contact area.
%   plot3d     - 3D plots of various quantities for the entire contact area.
%   plotstrs   - plot sub-surface displacements and stresses.
%
% Working with wheel/rail profiles.
%
%   read_profile    - generic routine for reading profiles
%   read_slices     - lower-level routine for reading variable profile slcs file
%   read_simpack    - lower-level routine for reading Simpack prr/prw files
%   read_miniprof   - lower-level routine for reading Miniprof ban/whl files
%   modify_profile  - lower-level routine for making some profile adjustments
%   plot_update     - show difference of two profiles at magnification factor
%   resample_slices - lower-level routine for 2D interpolation of variable profile
%   smooth_profile  - compute smoothed profile using parametric smoothing spline approximation
%   write_simpack   - routine for writing Simpack prr/prw files
%   write_miniprof  - routine for writing Miniprof ban/whl files
%
% Spline curves and spline surfaces.
%
%   make_spline         - compute parametric smoothing spline curve (profile)
%   eval_spline         - evaluate parametric smoothing spline curve (profile)
%   eval_spline_deriv   - evaluate derivatives of parametric smoothing spline curve
%   make_2dspline       - compute parametric smoothing spline surface (variable profile)
%   eval_2dspline       - evaluate parametric smoothing spline surface (variable profile)
%   eval_2dspline_deriv - evaluate derivatives of parametric smoothing spline surface
%   plot_2dspline       - create 3d surface plot of a parametric spline surface (variable profile)
%   helpers: make_arclength, make_ppform, make_reflec, eval_bspline_basisfnc, solve_cubic_eq

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
