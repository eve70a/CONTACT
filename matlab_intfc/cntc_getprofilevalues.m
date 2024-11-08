
%------------------------------------------------------------------------------------------------------------
% function [ values ] = cntc_getprofilevalues(ire, itask, iparam, rparam)
%
% Get a wheel or rail profile for a wheel-rail contact problem as a table of values
%
%  itask          - select type of outputs:
%                    -1: ierr   <0: error codes, 0=profile loaded ok, 1=not set
%                     0: npnt   number of points used in requested sampling method
%                     1: r/w    get (yr,zr) values for rail or (yw,zw) for wheel profile
%                     2: trk    get (ytr,ztr) values for rail or wheel profile (principal profile)
%                     3: gaug   get left-most point within gauge height, offset, point at gauge height
%                               [ygauge1, zgauge1, yoffs, zoffs, ygauge2, zgauge2]  [length]
%                     4: arc    get arc-length parameter s along profile
%                     5: angl   get surface inclination atan2(dz, dy) [angle]
%                     6: 2d-u   evaluate varprof [xyz]_r or [xyz]_w for cross-section at constant u_i = x_i
%                     7: 2d-v   evaluate varprof [xyz]_r or [xyz]_w for interpolation path at constant v_j
%                     8: 2d-y   evaluate varprof [xyz]_r or [xyz]_w for longitudinal slice at constant y_j
%                     9: spl2d  dump 2D spline data-structure for inspection in Matlab
%  iparam         - integer configuration parameters
%                     1: itype     0 = rail, 1 = wheel profile
%                     2: isampl   -1 = sampling cf. original input data;
%                                  0 = sampling cf. spline representation (default for tasks 0--5);
%                                  1 = sampling cf. spline representation at spacing ds_out (required for
%                                      tasks 6--8)
%                            kchk>=2 = sampling cf. spline representation with integer refinement factor
%  rparam         - real configuration parameters
%                     1: ds_out  step-size ds/du/dv used with sampling method isampl=1, default 1mm
%                     2: c1_out  sample position on variable profiles (task 1, 2: x_out)
%                     3: c2_sta  start of sample range on variable profiles
%                     4: c2_end  end of sample range on variable profiles
%                   task 6: c1 = u, c2 = v; task 7: c1 = v, c2 = u; task 8: c1 = y, c2 = u.
%
%  units: s,x,y,z [mm], th [rad], u,v [-]
%  tasks 1--9: no unit conversion or scaling are applied for profile values
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
   if (itask<-1 | itask>9)
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

   if (itask==-1)       % task -1: get error code

      lenarr = 2;
      p_val  = libpointer('doublePtr',zeros(lenarr,1));

      calllib(libname,'cntc_getprofilevalues_new', ire, itask, nints, iparam, nreals, rparam, lenarr, p_val);
      values = round(p_val.value);

   elseif (itask>=0 & itask<=8)

      % task 0--8: first get the number of points in the profile at requested s sampling positions

      itask_npnt = 0;
      if (itask>=6), itask_npnt = itask; end

      p_npnt = libpointer('doublePtr',-1);
      calllib(libname,'cntc_getprofilevalues_new', ire, itask_npnt, nints, iparam, nreals, rparam, 1, p_npnt);

      npnt = round(p_npnt.value);

      if (npnt==CNTC_err_profil | npnt<=0)
         disp(sprintf(['cntc_getprofilevalues: the rail and/or wheel profile could not be found or ',...
                                                                                'processed (%d).'], npnt));
         return;
      end

      % create output-array dependent on task to be performed

      if (itask==0)

         values = [ npnt ];

      elseif (itask>=1 & itask<=8)

         lenarr = [ npnt*2, npnt*2, 6, npnt, npnt, npnt*3, npnt*3, npnt*3, 1 ];
         lenarr = lenarr(itask);
         p_val  = libpointer('doublePtr',zeros(lenarr,1));

         calllib(libname,'cntc_getprofilevalues_new', ire, itask, nints, iparam, nreals, rparam, ...
                                                                                         lenarr, p_val);
         values = p_val.value;

         if (itask==1 | itask==2)
            values = reshape(values, npnt, 2);
         elseif (itask==6 | itask==7 | itask==8)
            values = reshape(values, npnt, 3);
            ix = find(values(:,1)> 98d9); values(ix,1) = NaN;
            ix = find(values(:,2)>998d0); values(ix,2) = NaN;
            ix = find(values(:,3)>998d0); values(ix,3) = NaN;
         end

      end % (itask==0 | 1--8)

   elseif (itask==9)    % dump spline structure in m-file

      lenarr = 1;
      p_val  = libpointer('doublePtr',zeros(lenarr,1));

      calllib(libname,'cntc_getprofilevalues_new', ire, itask, nints, iparam, nreals, rparam, lenarr, p_val);
      values = [];

   else

      values = [];

   end % if (itask==-1)

end % cntc_getprofilevalues

%------------------------------------------------------------------------------------------------------------

