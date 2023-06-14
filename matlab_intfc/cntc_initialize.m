
%------------------------------------------------------------------------------------------------------------
% function [ ifcver, ierror ] = cntc_initialize(ire, imodul, [outpath], [idebug])
%
% upon first call: initialize the addon internal data and initialize output channels,
%                  print version information;
% for each ire:   initialize and return the addon version number.
%
%  in:  integer    ire          - result element ID
%       integer    imodul       - module number 1=w/r contact, 3=basic contact
%       character  outpath(*)   - full path of output directory
%       integer    idebug       - show (1) or hide (0) error messages
%  out: integer    ifcver       - version of the CONTACT add-on
%       integer    ierror       - error flag
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 0: m=*, glob   - no icp needed

function [ ifcver, ierror ] = cntc_initialize(ire, imodul, c_outpath, idebug)
   global libname;
   CNTC_err_allow = -12;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(imodul))
      disp(sprintf('cntc_initialize: please select module 1 or 3'));
      return;
   end
   if (nargin<3 | isempty(c_outpath))
      c_outpath = ' ';
   end
   if (nargin<4 | isempty(idebug))
      idebug = 1;
   end

   p_ierr = libpointer('int32Ptr',-1);
   p_ver  = libpointer('int32Ptr',-1);

   p_ierr.value = 0;
   calllib(libname,'cntc_initialize', ire, imodul, p_ver, p_ierr, c_outpath, length(c_outpath));
   % disp(sprintf('test_caddon: obtained ver=%d, ierr=%d', p_ver.value, p_ierr.value));

   ifcver = double(p_ver.value);
   ierror = double(p_ierr.value);

   if (idebug>=1 & ierror==CNTC_err_allow)
      disp(sprintf('cntc_initialize: no license found or license invalid, check output-file (%d).',ierror));
   elseif (idebug>=1 & ierror<0)
      disp(sprintf('cntc_initialize: an error occurred in the CONTACT library (%d).',ierror));
   end

end % cntc_initialize

%------------------------------------------------------------------------------------------------------------

