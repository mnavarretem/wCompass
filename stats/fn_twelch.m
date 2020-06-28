function [mx_t,mx_df] = fn_twelch(mx_x,mx_y,nm_dim)
% [mx_t,mx_df] = fn_twelch(mx_x,mx_y,nm_dim) computes an unequal variances 
% t-test between variables mx_x and mx_y. nm_dim states the dimension to 
% compute the statisics (Default: 1). NaN values are ommited and discarded
% from the evaluation of the test. Outputs: mx_t: t-value; mx_df: degrees
% of freedom computed as Moser & Stevens (1992) and Ruxton (2006).

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
if nargin < 3
    nm_dim = 1;
end

if isempty(mx_x) || isempty(mx_x)
    mx_t = nan;
    mx_df = nan;
    return
end

if any(isnan(mx_x(:))) || any(isnan(mx_y(:)))    
    mx_uX1	= mean(mx_x,nm_dim,'omitnan');
    mx_uX2	= mean(mx_y,nm_dim,'omitnan');
    mx_vX1	= var(mx_x,0,nm_dim,'omitnan');
    mx_vX2	= var(mx_y,0,nm_dim,'omitnan');   
    mx_nX1	= sum(~isnan(mx_x),nm_dim);
    mx_nX2	= sum(~isnan(mx_y),nm_dim); 
else    
    mx_uX1	= mean(mx_x,nm_dim);
    mx_uX2	= mean(mx_y,nm_dim);
    mx_vX1	= var(mx_x,0,nm_dim);
    mx_vX2	= var(mx_y,0,nm_dim);
    mx_nX1	= size(mx_x,nm_dim);
    mx_nX2	= size(mx_y,nm_dim);
end

mx_u        = mx_vX2 ./ mx_vX1;

mx_t	= (mx_uX1-mx_uX2)./sqrt((mx_vX1./mx_nX1)+(mx_vX2./mx_nX2));
mx_df      = (((1./mx_nX1) + mx_u./mx_nX2).^2) ./...
            ((1./(mx_nX1.^2.*(mx_nX1-1))) + ((mx_u.^2)./(mx_nX2.^2.*(mx_nX2-1))));