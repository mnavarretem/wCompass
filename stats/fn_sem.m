function mx_sem = fn_sem(mx_data,nm_dim)
% mx_sem = fn_sem(mx_data,nm_dim) 
% Computes de squared error of the means of mx_data along the dimension
% nm_dim.

%% GNU licence,
% Copyright (C) <2016>  <Miguel Navarrete>
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
    if size(mx_data,1) == 1
        nm_dim = 2;
    else
        nm_dim = 1;
    end
end

if any(isnan(mx_data(:)))
    mx_sem	= std(mx_data,0,nm_dim,'omitnan')./sqrt(sum(~isnan(mx_data),nm_dim));
else
    mx_sem	= std(mx_data,0,nm_dim)./sqrt(size(mx_data,nm_dim));
end