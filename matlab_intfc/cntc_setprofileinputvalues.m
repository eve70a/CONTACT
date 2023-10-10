
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_setprofileinputvalues(ire, values, iparam, rparam)
%
% Set a wheel or rail profile for a wheel-rail contact problem using a table of values
%
%  values         - lists of profile points. Can be (npoint,2) or (2,npoint), or
%                        1d array of length (2*npoint), ordered [y1,z1, ... yi,zi, ... yn,zn]
%  iparam         - integer configuration parameters
%                     1: itype     0 = rail, 1 = wheel profile, -1 = taken from file extension (default)
%                     2:  -        not used
%                     3: mirrory   0 or -1 = no mirroring (default), 1 = mirror y coordinate values
%                     4: mirrorz   0 = autodetect (default), -1 = no mirroring, 1 = mirror z values
%                     5: errhndl   configuration of error handling. 
%                                   -2 = continue as much as possible, suppress error messages; 
%                                   -1 = suppress warnings; 0: warn and continue (default);
%                                    1 = signal errors and abort
%                     6: ismooth   selection of smoothing method. 0 = original smoothing spline (default),
%                                    1 = weighted PP smoothing spline, 2 = weighted smoothing B-spline (best)
%  rparam         - real configuration parameters
%                     1: sclfac    scaling factor for conversion to [mm], e.g. 1e3 for data given in [m]
%                                  default (sclfac<=0): using the active unit convention
%                     2: smooth    smoothing parameter lambda for non-weighted spline or l_filt for
%                                  weighted spline smoothing
%                     3: maxomit   fraction: signal error if more than maxomit of profile points are
%                                  discarded after cleanup of profile. Default 0.5, use 1 to disable check.
%                     4: zigthrs   angle threshold for zig-zag detection. Default 5/6*pi, >=pi to disable.
%                     5: kinkhigh  angle threshold for kink detection. Default pi/6, >=pi to disable.
%                     6: kinklow   angle threshold for neighbouring points in kink detection. 
%                                  default kinkhigh/5.
%                     7: kinkwid   half-width of window used for kink detection, [len], default 2 mm
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 2: m=1, wtd    - no icp needed

function [ ] = cntc_setprofileinputvalues(ire, values, iparam, rparam)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(values))
      disp('ERROR in cntc_setprofileinputvalues: values is mandatory.');
      return
   end
   if (nargin<3 | isempty(iparam))
      iparam = [  0 ];
   end
   if (nargin<4 | isempty(rparam))
      rparam = [ -1, 0. ];
   end

   % check sizes, reshape input values to 1-d array of 2*npoint elements

   arrsiz = size(values);
   if (min(arrsiz)==2) % 2-d array
      npoint = max(arrsiz);
      % transpose to (2,npoint) array
      if (arrsiz(1)>2)
         values = values';
      end
      % reshape to 1D array
      values = reshape(values, 2*npoint, 1);
   elseif (min(arrsiz)==1)
      npoint = max(arrsiz) / 2;
   else
      disp('ERROR in cntc_setprofileinputvalues: values must have 2*npoint values.');
      return
   end

   nints     = length(iparam);
   nreals    = length(rparam);
   calllib(libname,'cntc_setprofileinputvalues', ire, npoint, values, nints, iparam, nreals, rparam);

end % cntc_setprofileinputvalues

%------------------------------------------------------------------------------------------------------------

