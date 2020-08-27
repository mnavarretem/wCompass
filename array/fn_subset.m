function mx_subset = fn_subset(mx_X,nm_dim,vt_id)
%fn_subset	obtains a subset of a matrix mx_X form the index vt_id
%	mx_subset = fn_subset(mx_X,nm_dim,vt_id) obtains a subset mx_subset
%	of the matrix mx_X from the index vt_id along the dimension nm_dim 
%
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


if islogical(vt_id)
    vt_id = find(vt_id);
end

s_NumDim	= numel(size(mx_X));
v_StrCol    = repmat([{':'};{','}],1,s_NumDim);

v_StrCol{1,nm_dim}	= 'vt_id';

v_StrCol	= cell2mat(v_StrCol(:)');
v_StrCol	= v_StrCol(1:end-1);

mx_subset    = eval(sprintf('mx_X(%s)',v_StrCol));

end