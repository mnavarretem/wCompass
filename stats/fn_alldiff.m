function mx_diff = fn_alldiff(vt_x,vt_y)
% fn_alldiff computes all pairwise differences vt_x - vt_y. mx_diff is a
% matrix with numel(vt_x) x numel(vt_y) elements.

%% GNU licence,
% Copyright (C) <2017>  <Miguel Navarrete>
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

if isempty(vt_x) || isempty(vt_y) 
    mx_diff	= nan;
    return
end
    
vt_x	= vt_x(:);
vt_y	= vt_y(:)';

nm_numX	= numel(vt_x);
nm_numY	= numel(vt_y);

vt_x    = repmat(vt_x,1,nm_numY);
vt_y    = repmat(vt_y,nm_numX,1);

mx_diff = vt_x - vt_y;