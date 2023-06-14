function contains = myisfield(STRUC,fieldname)
%
% function [contains] = myisfield(STRUC,fieldname)
%
% An extension of the function isfield, for which the fieldname may contain '.'
% Example:
%    STRUC.substruc.field = [0 1];
%    myisfield(STRUC,'substruc.field')   returns 1
%    myisfield(STRUC,'substruc')         returns 1,
%    myisfield(STRUC,'')                 returns 1,
%    myisfield(STRUC,any_other_string)   returns 0.
%
%   1.0     2004-07-28   BvtH (VORtech)   initial version

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% Initialize: assume that the field is in te struct.

contains   = 1;
fieldname  = deblank(fieldname);
structname ='STRUC';

% Check all the 'words' in the fieldname

while (contains & ~strcmp(fieldname,''))
   [word,fieldname] = strtok(fieldname,'.');

   % if the field name contains words after this one,
   %  remove the '.' at the start (it is not removed by strtok)
   if (length(fieldname)>0)
      fieldname = fieldname(2:end);
   end

   fieldname = deblank(fieldname);
   word      = deblank(word);

   % Check if the struct contains a field whose name corresponds to
   % the word
   contains = contains & eval(['isfield(',structname,',word)']);

   % Continue checking inside the field.
   structname = [structname,'.',word];
end

