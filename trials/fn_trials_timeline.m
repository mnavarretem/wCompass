function vt_trials	= fn_trials_timeline(vt_signal,st_cnf)
%  mx_trials = fn_trials_timeline(vt_signal,st_cnf) get timeline trials of
%  vt_signal as specified in st_cnf structure. vt_trials is a cell vector
%  composed by the trials selected
%
% st_cnf corresponds to a structure with the fields:
%
%   - st_cnf.markers: 	vector indicating the time reference to each
%                      	trial
%
%   - st_cnf.fsample:	sampling frequency of the signal
%   
%   - st_cnf.window:    1 x 2 vector indicating the begining and end
%                       time (in seconds) of each trial referred to the
%                       marker position (marker = 0)
%   
%   - st_cnf.timetrial:	relative time for trials (in seconds) 

%% GNU licence,
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


%% Code starts here

%%	- Check default inputs

if ~isvector(vt_signal)
    error('fn_trials_timeline:signalIsNotVector','Input must be a vector;'); 
end
    
if nargin < 2   
    error('fn_detectsleepSO:wrongStruct','Second Input must be an structure;'); 
end
    
% Check whether st_cnf fileds need defaults
if ~isfield(st_cnf,'timetrial') 
    st_cnf.timetrial	= st_cnf.window(1):1/st_cnf.fsample:st_cnf.window(2);
end

%%	- Extract trials

vt_signal	= single(vt_signal(:)');
vt_markers	= round(st_cnf.markerTime .* st_cnf.fsample);

if isempty(vt_markers)
    vt_trials	= [];
    return
end

vt_tTrial   = st_cnf.window(1):1/st_cnf.fsample:st_cnf.window(2);

mx_window	= round(st_cnf.window .* st_cnf.fsample);
mx_window   = mx_window(1):mx_window(2);

mx_window   = repmat(mx_window,numel(vt_markers),1);
vt_markers  = repmat(vt_markers,1,size(mx_window,2));

mx_window   = mx_window + vt_markers;

vt_idx      = sum(numel(vt_signal) < mx_window,2) == 0;
mx_window   = mx_window(vt_idx,:); 

vt_trials   = vt_signal(mx_window);

if ~isequal(vt_tTrial,st_cnf.timetrial)
    vt_trials	= vt_trials';
    vt_trials	= interp1(vt_tTrial,vt_trials,st_cnf.timetrial);
    vt_trials   = vt_trials';
end

vt_trials	= mat2cell(vt_trials,ones(size(vt_trials,1),1),size(vt_trials,2));

