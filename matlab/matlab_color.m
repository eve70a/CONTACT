function [ rgb ] = matlab_color( num )

% function [ rgb ] = matlab_color( num )
%
% return the rgb-triplet for color order index num

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

matlab_colors = [
         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840
];
num_colors = size(matlab_colors,1);

if (nargin<1 | isempty(num) | ~isnumeric(num))
   num = 1;
end

% reshape to column vector
num = reshape(num, prod(size(num)), 1);

% clip
num = mod(num-1, num_colors) + 1;

% get rows
rgb = matlab_colors(num, :);

