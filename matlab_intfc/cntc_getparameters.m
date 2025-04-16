
%------------------------------------------------------------------------------------------------------------
% function [ values ] = cntc_getparameters(ire, icp, itask)
%
% used for retrieving various parameters from a contact problem
%
%  itask          - selected group of parameters (1: cntc_getcpresults, 2: material, 3: friction)
%  values         - values of the parameters obtained from CONTACT
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 6: m=*, cp     - require icp>0, default 1

function [ values ] = cntc_getparameters(ire, icp, itask)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = 1;
   end
   if (nargin<3 | isempty(itask))
      itask = 1; % values used in plot3d
   end

   lenarr = 30;
   p_values = libpointer('doublePtr',zeros(lenarr,1));

   calllib(libname,'cntc_getparameters', ire, icp, itask, lenarr, p_values);

   tmp    = p_values.value;
   values = struct();
   if (itask==1)
      values.veloc  = tmp( 1);
      values.chi    = tmp( 2);
      values.dq     = tmp( 3);
      values.spinxo = tmp( 4);
      values.spinyo = tmp( 5);
      values.tau_c0 = tmp( 6);
   elseif (itask==2)
      values.gg1    = tmp( 1);
      values.gg2    = tmp( 2);
      values.gg     = tmp( 3);
      values.poiss1 = tmp( 4);
      values.poiss2 = tmp( 5);
      values.poiss  = tmp( 6);
      values.ak     = tmp( 7);
      values.flx1   = tmp( 8);
      values.flx2   = tmp( 9);
      values.flx3   = tmp(10);
      values.k0_mf  = tmp(11);
      values.alfamf = tmp(12);
      values.betamf = tmp(13);
      values.k_eff  = tmp(14);
      values.gg3    = tmp(15);
      values.laythk = tmp(16);
      values.tau_c0 = tmp(17);
      values.k_tau  = tmp(18);
      values.cdampn = tmp(19);
      values.cdampt = tmp(20);
      values.dfnmax = tmp(21);
      values.dftmax = tmp(22);
   elseif (itask==3)
      disp(tmp(1:10)')
      values.lmeth  = tmp( 1);
      values.nvf    = tmp( 2);
      values.memdst = tmp( 3);
      values.mem_s0 = tmp( 4);
      values.alpha  = tmp( 5+3*[0:values.nvf-1]);
      values.fstat  = tmp( 6+3*[0:values.nvf-1]);
      values.fkin   = tmp( 7+3*[0:values.nvf-1]);
   end

end % cntc_getparameters

%------------------------------------------------------------------------------------------------------------

