function mx_perm = fn_uperms(nm_nvar,nm_numper,nm_numS1)
% fn_uperms computes non-repeated nm_numper permutations for nm_numel
% variables and a bllock of nm_numS1 values

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

mx_perm	= zeros(nm_numper, nm_nvar, 'uint32');

vt_rIdx        = 1:nm_nvar;
mx_perm(1,:)	= vt_rIdx;

nm_count     = 1;
nm_permBlock = ceil(nm_numper/10);

while nm_count < nm_numper
    [~,mx_id]	= sort(rand(nm_permBlock,nm_nvar),2);
    mx_newIdx	= [mx_perm(1:nm_count, :); vt_rIdx(mx_id)];
    
    [~,vt_idC1]	= unique(sort(mx_newIdx(:,1:nm_numS1),2), 'rows','stable');
    [~,vt_idC2] = unique(sort(mx_newIdx(:,nm_numS1+1:end),2), 'rows','stable');
    
    vt_idStore          = false(size(mx_newIdx,1),1);
    vt_idStore(vt_idC1)	= true;
    vt_idC1             = vt_idStore;
    
    vt_idStore          = false(size(mx_newIdx,1),1);
    vt_idStore(vt_idC2)	= true;
    vt_idC2             = vt_idStore;
    
    clear v_IdxStore
    
    mx_newIdx    = mx_newIdx(vt_idC1 & vt_idC2,:);
    
    if size(mx_newIdx,1) > nm_count
        
        mx_perm	= mx_newIdx;
        nm_count     = size(mx_newIdx,1);
        
    end
end

if size(mx_perm,1) > nm_numper
    mx_perm  = mx_perm(1:nm_numper,:);
end

end