
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_closelibrary();
%
% clean-up, close files and unload the library from Matlab
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2026 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 0: m=*, glob   - no icp needed

function [ ] = cntc_closelibrary();

   global libname;

   if (libisloaded(libname))
      calllib(libname,'cntc_finalizelast');

      unloadlibrary(libname);
   end

end % function cntc_closelibrary

%------------------------------------------------------------------------------------------------------------

