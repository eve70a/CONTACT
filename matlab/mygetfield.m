function result = mygetfield(struct,field);
%
% function [result] = mygetfield(struct,field);
%
% returns the field <field> from the struct <struct>,
% even if this field contains special characters.
%
% INTENDED USE:
%      a = mygetfield(struct,'substruct.field');
%
% Other possible uses (not intended):
%      a = mygetfield(struct,'substruct.field(3,4)');
%

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

result = eval(['struct.',field]);

