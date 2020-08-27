function [mx_p,mx_t,mx_df] = fn_pttest(mx_x,mx_y,nm_tails,nm_dim)
% [mx_p,mx_t,mx_df] = fn_pttest(mx_x,mx_y,nm_tails,nm_dim) computes
% a paired t-test between variables mx_x and mx_y. Tails of the
% test are determined by nm_tails (One-tailed: left = -1,rigth = 1; 
% Two-tailed = 0; Default = 0). nm_dim states the dimension to compute the
% statisics (Default: 1). Outputs: mx_p: p-value; mx_t: t-value; mx_df:
% degrees of freedom.
%

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
if nargin < 3
    nm_tails = 0;
end

if nargin < 4
    nm_dim = 1;
end

mx_diff	=  mx_x - mx_y;
mx_nans = isnan(mx_diff);
vt_size = sum(~mx_nans,nm_dim);    
mx_df   = vt_size - ones(size(vt_size));

mx_df(mx_df < 0)	= 0;

mx_dmean	= mean(mx_diff,nm_dim,'omitnan');
mx_dstd     = std(mx_diff,[],nm_dim,'omitnan');
mx_sqrtn	= sqrt(vt_size);

mx_t = mx_dmean ./ (mx_dstd ./ mx_sqrtn);

switch nm_tails
    case -1 % left one-tailed test
        mx_p	= tcdf(mx_t,mx_df);
    case 0 % two-tailed test
        mx_p	= 2*tcdf(-abs(mx_t),mx_df);
    case 1 % right one-tailed test
        mx_p	= tcdf(-mx_t,mx_df);
    otherwise
        error('[fn_pttest] - tails error')
end