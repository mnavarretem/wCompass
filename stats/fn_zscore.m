function [mx_m,st_stat] = fn_zscore(mx_x,nm_dim,nm_robust)
% mx_m = fn_zscore(mx_x,nm_dim,nm_robust) computes the zscore along nm_dim 
% dimmension with nm_robust indicating the use of fobust statistics 
% (defalut = 0, 'no robust')

%% GNU licence,
% Copyright (C) <2018>  <Miguel Navarrete>
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
if nargin < 2
    nm_dim	= 1;
end
if nargin < 3
    nm_robust	= 0;
end

if nm_robust
    mx_median	= median(mx_x,nm_dim,'omitnan');
    mx_mad      = median(abs(mx_x - mx_median),nm_dim,'omitnan');

    mx_m	= (mx_x - mx_median) ./ (mx_mad * 1.4826);
    
    nm_zcenter      = mx_median;
    nm_zdeviation	= mx_mad;
else    
    mx_mean	= mean(mx_x,nm_dim,'omitnan');
    mx_std	= std(mx_x,0,nm_dim,'omitnan');

    mx_m	= (mx_x - mx_mean) ./ mx_std;
    
    nm_zcenter      = mx_mean;
    nm_zdeviation	= mx_std;
end

if nargout > 1
    st_stat.zcenter         = nm_zcenter;
    st_stat.zdeviation	= nm_zdeviation;
end