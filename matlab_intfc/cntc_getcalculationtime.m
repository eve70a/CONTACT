
%------------------------------------------------------------------------------------------------------------
% function [ tcpu, twall ] = cntc_getcalculationtime(ire, icp)
%
% return accumulated cpu-time and wall-clock-time used since last timer reset for a contact problem
%  tcpu, twall   - cpu- and wall-clock times used [time]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 6: m=*, cp     - require icp>0, default 1

function [ tcpu, twall ] = cntc_getcalculationtime(ire, icp)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = 1;
   end
   if (icp<=0)
      disp(sprintf('ERROR in cntc_getcalculationtime: not available for icp=%d',icp));
      return
   end

   p_tcpu  = libpointer('doublePtr',-1);
   p_twall = libpointer('doublePtr',-1);

   calllib(libname,'cntc_getcalculationtime', ire, icp, p_tcpu, p_twall);

   tcpu  = p_tcpu.value;
   twall = p_twall.value;

end % cntc_getcalculationtime

%------------------------------------------------------------------------------------------------------------

