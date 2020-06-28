function vt_trials	= fn_trials_timeline(vt_signal,st_Cnf)
%  mx_trials = fn_extracttimetrials(vt_signal,st_Cnf) get time trials from
%  vt_signal as specified in st_Cnf structure. vt_trials is a cell vector
%  composed by the trials selected
%
% st_Cnf corresponds to a structure with the
% fields:
%
%   - st_Cnf.markers:       vector indicating the time reference to each
%                           trial
%
%   - st_Cnf.fsampling:     sampling frequency of the signal
%   
%   - st_Cnf.window:        1 x 2 vector indicating the begining and end
%                           time (in seconds) of each trial referred to the
%                           marker position
%   
%   - st_Cnf.timetrial:     time vector (in seconds) in which the output
%                           should be referenced. 

% Miguel Navarrete
% CUBRIC
% 2017

%% Code starts here

%%	- Check default inputs

if ~isvector(vt_signal)
    error('fn_extracttimetrials:signalIsNotVector','Input must be a vector;'); 
end
    
if nargin < 2   
    error('fn_detectsleepSO:wrongStruct','Second Input must be an structure;'); 
end
    
% Check whether st_Cnf fileds need defaults
if ~isfield(st_Cnf,'timetrial') 
    st_Cnf.timetrial	= st_Cnf.window(1):1/st_Cnf.fsampling:st_Cnf.window(2);
end

%%	- Extract trials

vt_signal	= single(vt_signal(:)');
vt_markers	= round(st_Cnf.markerTime .* st_Cnf.fsampling);

if isempty(vt_markers)
    vt_trials	= [];
    return
end

vt_timeTrial= st_Cnf.window(1):1/st_Cnf.fsampling:st_Cnf.window(2);

mx_window	= round(st_Cnf.window .* st_Cnf.fsampling);
mx_window   = mx_window(1):mx_window(2);

mx_window   = repmat(mx_window,numel(vt_markers),1);
vt_markers  = repmat(vt_markers,1,size(mx_window,2));

mx_window   = mx_window + vt_markers;

vt_Idx      = sum(numel(vt_signal) < mx_window,2) == 0;
mx_window   = mx_window(vt_Idx,:); 

vt_trials   = vt_signal(mx_window);

if ~isequal(vt_timeTrial,st_Cnf.timetrial)
    vt_trials	= vt_trials';
    vt_trials	= interp1(vt_timeTrial,vt_trials,st_Cnf.timetrial);
    vt_trials   = vt_trials';
end

vt_trials	= mat2cell(vt_trials,ones(size(vt_trials,1),1),size(vt_trials,2));

