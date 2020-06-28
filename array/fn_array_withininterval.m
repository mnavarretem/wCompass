function [vt_id] = fn_array_withininterval(vt_x,mx_int)
% [vt_id] = fn_withininterval(vt_x,mx_int) checks if vt_x is within
% intervals of mx_int. mx_int must be a n x 2 matrix with mx_int(:,1) 
% indicating the begining and mx_int(:,1) the end of intervals  
%
% Copyright (C) <2020>  <Miguel Navarrete>
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

%% Code starts here
vt_x	= vt_x(:);
nm_x    = numel(vt_x);
nm_int	= size(mx_int,1);

mx_x    = repmat(vt_x,1,nm_int);
mx_beg	= repmat(mx_int(:,1)',nm_x,1);
mx_end	= repmat(mx_int(:,2)',nm_x,1);

vt_id   = mx_beg <= mx_x & mx_x <= mx_end;
vt_id   = any(vt_id,2);
