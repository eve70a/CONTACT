function s = mysetfield(s,varargin)
%
% function [s] = mysetfield(s, varargin)
%
% Set structure field contents.
%
%   S = MYSETFIELD(S,'field',V) sets the contents of the specified
%   field to the value V.  This is equivalent to the syntax S.field = V.
%   S must be an n-by-1 array of structures with n-by-1 array V.
%   The changed structures S are returned.
%
%   S = MYSETFIELD(S,{i,j},'field',{k},V) is equivalent to the syntax
%       S(i,j).field(k) = V;
%   In other words, S = SETFIELD(S,sub1,sub2,...,V) sets the
%   contents of the structure S to V using the subscripts or field
%   references specified in sub1,sub2,etc.  Each set of subscripts in
%   parentheses must be enclosed in a cell array and passed to
%   SETFIELD as a separate input.  Field references are passed as
%   strings.
%
%   See also MYGETFIELD, GETFIELD, RMFIELD, MYISFIELD, ISFIELD, FIELDNAMES.

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% Check for sufficient inputs
if (isempty(varargin) | length(varargin) < 2)
    error('MYSETFIELD:InsufficientInputs', 'Not enough input arguments.');
end

% The most common case
arglen = length(varargin);
strField = varargin{1};
if (arglen==2)
    len = length(s);
    if (len<=1)
        eval(['s.',deblank(strField),' = varargin{end};']);
    else
       val = varargin{end};
       for i = 1 : len
          eval(['s(i).',deblank(strField),'= val(i);']);
       end
    end
    return
end


subs = varargin(1:end-1);
for i = 1:arglen-1
    index = varargin{i};
    if (isa(index, 'cell'))
        types{i} = '()';
    elseif isstr(index)
        types{i} = '.';
        subs{i} = deblank(index); % deblank field name
    else
        error('MYSETFIELD:InvalidType','Inputs must be either cell arrays or strings.');
    end
end

% Perform assignment
try
   s = builtin('subsasgn', s, struct('type',types,'subs',subs), varargin{end});
catch
   error('MYSETFIELD', lasterr)
end

