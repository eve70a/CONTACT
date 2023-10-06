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
%   resample_slices - helper routine for 2D interpolation of variable profile
%   write_simpack   - routine for writing Simpack prr/prw files
%   write_miniprof  - routine for writing Miniprof ban/whl files

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
