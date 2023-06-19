
%------------------------------------------------------------------------------------------------------------
% function [ values ] = cntc_getprofilevalues(ire, itask, iparam, rparam)
%
% Get a wheel or rail profile for a wheel-rail contact problem as a table of values
%
%  itask          - select type of outputs:
%                     0: npnt   number of points used in requested sampling method
%                     1: r/w    get (yr,zr) value for rail or (yw,zw) for wheel profile
%                     2: trk    get (ytr,ztr) values for rail or wheel profile (principal profile)
%                     3: gaug   get left-most point within gauge height, offset, point at gauge height
%                               [ygauge1, zgauge1, yoffs, zoffs, ygauge2, zgauge2]  [length]
%                     4: arc    get arc-length parameter s along profile
%                     5: angl   get surface inclination atan2(dz, dy) [angle]
%  iparam         - integer configuration parameters
%                     1: itype     0 = rail, 1 = wheel profile
%                     2: isampl   -1 = sampling cf. original input data;
%                                  0 = sampling cf. spline representation (default);
%                                  1 = sampling cf. spline representation at spacing ds_out
%                            kchk>=2 = sampling cf. spline representation with integer refinement factor
%  rparam         - real configuration parameters
%                     1: ds_out  step-size ds used with sampling method isampl=1, default 1mm
%
%  tasks 1,2,4: no unit conversion or scaling are applied for profile values
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 2: m=1, wtd    - no icp needed

function [ values ] = cntc_getprofilevalues(ire, itask, iparam, rparam)
   global libname;
   CNTC_err_profil = -32;

   values = [];
   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(itask))
      disp('ERROR in cntc_getprofilevalues: itask is mandatory.');
      return
   end
   if (itask<0 | itask>5)
      disp(sprintf('ERROR in cntc_getprofilevalues: itask=%d not available.', itask));
      return
   end
   if (nargin<3 | isempty(iparam))
      itype = 0; isampl = 0;
      iparam = [ itype, isampl ];
   end
   if (nargin<4 | isempty(rparam))
      ds_out = 1.0;
      rparam = [ ds_out ];
   end
   nints  = length(iparam);
   nreals = length(rparam);

   % always get the number of points in the profile at requested s sampling positions

   p_npnt = libpointer('doublePtr',-1);
   calllib(libname,'cntc_getprofilevalues_new', ire, 0, nints, iparam, nreals, rparam, 1, p_npnt);

   npnt = round(p_npnt.value);

   if (npnt==CNTC_err_profil | npnt<=0)
      disp(sprintf('cntc_getprofilevalues: the rail and/or wheel profile could not be found or processed (%d).',npnt));
      return;
   end

   % create output-array dependent on task to be performed

   if (itask==0)

      values = [ npnt ];

   elseif (itask>=1 & itask<=5)

      lenarr = [ npnt*2, npnt*2, 6, npnt, npnt ];
      lenarr = lenarr(itask);
      p_val  = libpointer('doublePtr',zeros(lenarr,1));

      calllib(libname,'cntc_getprofilevalues_new', ire, itask, nints, iparam, nreals, rparam, lenarr, p_val);
      values = p_val.value;

      if (itask==1 | itask==2)
         values = reshape(values, npnt, 2);
      end

   else

      values = [];

   end

end % cntc_getprofilevalues

%------------------------------------------------------------------------------------------------------------

