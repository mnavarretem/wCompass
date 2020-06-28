function [mx_p,mx_t,mx_df] = fn_uttest(mx_x,mx_y,nm_tails,nm_dim)
% [mx_p,mx_t,mx_df] = fn_uttest(mx_x,mx_y,nm_tails,nm_dim) computes
% an unequal variances t-test between variables mx_x and mx_y. Tails of the
% test are determined by nm_tails (One-tailed: left = -1,rigth = 1; 
% Two-tailed = 0; Default = 0). nm_dim states the dimension to compute the
% statisics (Default: 1). Outputs: mx_p: p-value; mx_t: t-value; mx_df:
% degrees of freedom computed as Moser & Stevens (1992) and Ruxton (2006).
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

[mx_t,mx_df]   = fn_twelch(mx_x,mx_y,nm_dim);
switch nm_tails
    case -1 % left one-tailed test
        mx_p	= tcdf(mx_t,mx_df);
    case 0 % two-tailed test
        mx_p	= 2*tcdf(-abs(mx_t),mx_df);
    case 1 % right one-tailed test
        mx_p	= tcdf(-mx_t,mx_df);
    otherwise
        error('[fn_uttest] - tails error')
end