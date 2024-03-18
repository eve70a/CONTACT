
%------------------------------------------------------------------------------------------------------------
% function [ ierror ] = cntc_readinpfile(ire, inp_type, fname)
%
% read settings from inp-file
%
%  in:  integer    ire          - result element ID
%       integer    inp_type     - type of inp-file: CNTC_inp_spck, ...
%       character  fname(*)     - filename
%  out: integer    ierror       - error flag
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 2: m=1, wtd    - no icp needed

function [ ierror ] = cntc_readinpfile(ire, inp_type, c_fname)
   global libname;
   CNTC = cntc_getmagicnumbers();

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2)
      inp_type = [];
   end
   if (~any(inp_type==[CNTC.inp_spck]))
      disp('ERROR(readinpfile): inp_type must be CNTC.inp_spck');
      return;
   end
   if (nargin<3 | isempty(c_fname))
      disp('ERROR(readinpfile): fname is mandatory')
      return;
   end

   p_ierr = libpointer('int32Ptr',-1);

   p_ierr.value = 0;
   calllib(libname,'cntc_readinpfile', ire, inp_type, c_fname, length(c_fname), p_ierr);

   ierror = double(p_ierr.value);

   if (ierror~=0)
      disp(sprintf('cntc_readinpfile: an error occurred in the CONTACT library (%d).',ierror));
   end

end % cntc_readinpfile

%------------------------------------------------------------------------------------------------------------

