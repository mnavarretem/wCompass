function [mx_dat16,vt_scale]	= fn_scale_int16(mx_dat)
% [mx_dat,vt_scale]	= fn_scale_int16(mx_dat) converts mx_dat to int16 by 
% mataining a linear relationship suc as mx_dat = polyval(vt_scale,mx_dat16)
% 
%% GNU licence,
% Copyright (C) <2021>  <Miguel Navarrete>
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
nm_resolution	= 2^16;

mx_scale	= [min(mx_dat(:)),max(mx_dat(:));0,nm_resolution];
vt_linCoef  = fn_linecoef(mx_scale(1,:),mx_scale(2,:));
mx_dat      = mx_dat.*vt_linCoef(1) + vt_linCoef(2);

mx_dat16	= uint16(mx_dat);
vt_scale	= fn_linecoef(mx_scale(2,:),mx_scale(1,:));