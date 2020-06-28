function vt_binCenter = fn_edge2bin(vt_binEdges)
% vt_binCenter = fn_edge2bin(vt_binEdges) converts bin edges to bin centers
% ignoring nans 

%% GNU licence,
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
vt_binEdges	= vt_binEdges(~isinf(vt_binEdges));
vt_binDiff	= diff(vt_binEdges)/2;
vt_binEdges = vt_binEdges(1:(numel(vt_binEdges)-1));
vt_binCenter= vt_binEdges + vt_binDiff;
