
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_setmaterialparameters(ire, icp, m_digit, rparam)
%
% set the material parameters for a contact problem
%
% M-digit  = 0: rparam = [ poiss1, poiss2, g1,  g2                                 ];
%            1: rparam = [ poiss1, poiss2, g1,  g2,  fg1,   fg2,    vt1,    vt2    ];
%            2: rparam = [ poiss1, poiss2, g1,  g2,  flx,   k0_mf,  alfamf, betamf ];
%            3: rparam = [ poiss1, poiss2, g1,  g2,  k0_mf, alfamf, betamf         ];
%            4: rparam = [ poiss1, poiss2, g1,  g2,  g3,    laythk, tau_c0, k_tau  ];
% values >= 10 are used to configure the M2-digit, M2 = m_digit - 10
% M2-digit = 2: rparam = [ cdampn, cdampt, dfnmax, dftmax ];
%
% dimensions: poiss1, poiss2  [-],        g1, g2    [force/area],   
%             fg1, fg2        [-],        vt1, vt2  [time],
%             flx         [volume/force], k0_mf, alfamf, betamf [-], 
%             g3, tau_c0  [force/area],   laythk    [length],        k_tau [force/volume]
%             cdampn, cdampt  [time],     dfnmax, dftmax   [force/time]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 5: m=*, wtd    - default icp=-1

function [ ] = cntc_setmaterialparameters(ire, icp, m_digit, rparam)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = -1; % default: W/R contact, all patches
   end
   if (nargin<3 | isempty(m_digit))
      m_digit = 0; % default: linearly elastic material
   end
   if (nargin<4 | isempty(rparam))
      disp('ERROR in cntc_setmaterialparameters: rparam is mandatory.');
      return
   end

   calllib(libname,'cntc_setmaterialparameters', ire, icp, m_digit, length(rparam), rparam);

end % cntc_setmaterialparameters

%------------------------------------------------------------------------------------------------------------

