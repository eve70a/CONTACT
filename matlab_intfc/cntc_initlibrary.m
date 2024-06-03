
%------------------------------------------------------------------------------------------------------------
% function [ CNTC, ifcver, ierror ] = cntc_initlibrary(wrkdir, outdir, expnam, idebug);
%
% load the library into Matlab, initialize its internal data and output channels
%
%  in:  character  wrkdir(*)    - [optional] effective working folder
%       character  outdir(*)    - [optional] output folder
%       character  expnam(*)    - [optional] experiment name, default 'contact_addon'
%       integer    idebug       - [optional] show (1) or hide (0) error messages
%  out: integer    CNTC         - struct with 'magic numbers' for configuring CONTACT
%       integer    ifcver       - version of the CONTACT add-on
%       integer    ierror       - error flag
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 0: m=*, glob   - no icp needed

function [ CNTC, ifcver, ierror ] = cntc_initlibrary(c_wrkdir, c_outdir, c_expnam, idebug);

   global libname;

   % retrieve current directory, form library name
   [pathstr, name, ext] = fileparts(which('cntc_initlibrary'));

   if (isunix())
      libname='contact_addon_linux64';
      lib_ext='.so';
   else
      libname=['contact_addon_', computer('arch')];
      lib_ext='.dll';
   end

   % unload the library if it was loaded before
   if (1==1 & libisloaded(libname))
      calllib(libname,'cntc_finalizelast');
      unloadlibrary(libname);
   end

   % load the library into Matlab
   if (~libisloaded(libname))
      pathstr = [deblank(pathstr) filesep '..' filesep 'bin'];
      fullname = fullfile(pathstr, [libname lib_ext]);

      if (exist(fullname)~=2)
         disp(['ERROR: cant find library: ', fullname]);
         return
      else
         loadlibrary(fullname, 'contact_addon.h');
      end
   end

   % initialize the internal data of the library, open its output streams
   if (nargin<1 | isempty(c_wrkdir))
      c_wrkdir = ' ';
   end
   if (nargin<2 | isempty(c_outdir))
      c_outdir = ' ';
   end
   if (nargin<3 | isempty(c_expnam))
      c_expnam = ' ';
   end
   if (nargin<4 | isempty(idebug))
      idebug = 1;
   end
   len_wrkdir = length(c_wrkdir);
   len_outdir = length(c_outdir);
   len_expnam = length(c_expnam);
   ioutput = 0;

   p_ierr = libpointer('int32Ptr',-1);
   p_ver  = libpointer('int32Ptr',-1);

   p_ierr.value = 0;
   calllib(libname,'cntc_initializefirst_new', p_ver, p_ierr, ioutput, c_wrkdir, c_outdir, c_expnam, ...
                                                                    len_wrkdir, len_outdir, len_expnam);
   % disp(sprintf('test_caddon: obtained ver=%d, ierr=%d', p_ver.value, p_ierr.value));

   ifcver = double(p_ver.value);
   ierror = double(p_ierr.value);

   % return a struct with 'magic numbers' for setting flags later on

   CNTC = cntc_getmagicnumbers();

   if (idebug>=1 & ierror==CNTC.err_allow)
      disp(sprintf('cntc_initlibrary: no license found or license invalid, check output-file (%d).',ierror));
   elseif (idebug>=1 & ierror<0)
      disp(sprintf('cntc_initlibrary: an error occurred in the CONTACT library (%d).',ierror));
   end

end % function cntc_initlibrary

%------------------------------------------------------------------------------------------------------------

