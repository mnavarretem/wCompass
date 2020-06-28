function vt_coef = fn_linecoef(vt_x,vt_y)
% vt_coef = fn_linecoef(vt_x,vt_y)) computes linear coeficients from two
% points in the space to satisfy y = mx + b from the form 
% [(y - y1) / (x - x1)] = [(y2 - y1) / (x2 - x1)].
% 
% Inputs 
% vt_x	= [x1,x2]
% vt_y	= [y1,y2]
% 
% Output
% vt_coef	= [m,b]

%% GNU licence,
% Copyright (C) <2015>  <Miguel Navarrete>
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
vt_coef(1) = (vt_y(2)-vt_y(1))/(vt_x(2)-vt_x(1));
vt_coef(2) = vt_y(1) - vt_x(1)*vt_coef(1);
