
%------------------------------------------------------------------------------------------------------------
% function [ sens ] = cntc_getsensitivities(ire, icp, lenout, lenin)
%
% return the sensitivities of the total forces for a contact problem
%
%  lenout, lenin      - requested number of outputs (forces) and inputs (creepages or shifts)
%  sens(lenout,lenin) - matrix of sensitivities
%
%  the inputs are ordered   1: pen, 2: cksi, 3: ceta, 4: cphi
%  in rolling, T=2,3, the units are pen [length], cksi, ceta [-],      cphi [angle/length]
%  in shifts,  T=1,   the units are pen [length], cksi, ceta [length], cphi [angle]
%
%  the outputs are ordered  1: fn,  2: fx,   3: fy,   4: mz
%  the units are fn [force], fx, fy [-], mz [force.length]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 6: m=*, cp     - require icp>0, default 1

function [ sens ] = cntc_getsensitivities(ire, icp, lenout, lenin)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = 1;
   end
   if (icp<=0)
      disp(sprintf('ERROR in cntc_getsensitivities: not available for icp=%d',icp));
      return
   end
   if (nargin<4 | isempty(lenout) | isempty(lenin))
      disp('ERROR in cntc_getsensitivities: lenout, lenin are mandatory.');
      return
   end

   lenarr = lenout * lenin;
   p_sens = libpointer('doublePtr',zeros(lenarr,1));

   calllib(libname,'cntc_getsensitivities', ire, icp, lenout, lenin, p_sens);

   sens = reshape(p_sens.value, lenout, lenin);

end % cntc_getsensitivities

%------------------------------------------------------------------------------------------------------------

