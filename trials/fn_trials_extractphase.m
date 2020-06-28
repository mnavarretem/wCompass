function vt_trials	= fn_trials_extractphase(vt_signal,vt_phase,st_cnf)
%  mx_trials = fn_trials_extractphase(vt_signal,vt_phase,st_cnf) get phase 
%  trials from vt_signal depending in phase specified in vt_phase and
%  settings form st_cnf structure. vt_trials is a cell vector composed 
%  by the selected trials
%
% st_cnf corresponds to a structure with the fields:
%
%   - st_cnf.markers:       Vector indicating the sample markers defining
%                           each trial <triggers>
%
%   - st_cnf.fsample:       Sampling frequency of the signal
%   
%   - st_cnf.window:        1 x 2 vector indicating the begining and end
%                           phase (in radians) of each trial surrounding
%                           each marker position
%   
%   - st_cnf.alignPhase:	Array indicating the samples in wich the phase
%                           of the trial should be aligned
%
%   - st_cnf.phasebins:     Phase array determining the output bins  

%   - st_cnf.border:        Phase 'lag' to include around the phase window

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

%%	- Check default inputs
    
if nargin < 3   
    error('fn_trials_extractphase:wrongInputs','function requires three inputs'); 
end

% Check whether st_cnf fileds need defaults
if ~isfield(st_cnf,'alignPhase') 
    st_cnf.alignPhase	= st_cnf.markers;
end

%%	- Extract trials
vt_signal	= vt_signal(:)';
vt_uwPhase	= unwrap(vt_phase);
vt_phBins	= st_cnf.phasebins;
nm_phRes	= diff(vt_phBins);
nm_phRes	= nm_phRes(1);

vt_trials	= cell(numel(st_cnf.markers),1);
 
for kk = 1:numel(st_cnf.markers)
    
    nm_curPhSample	= unique(st_cnf.alignPhase(kk,:));
    nm_curPhValue   = vt_phase(nm_curPhSample);
    
    % NOTE: Only uses one value without linearizing multiple markers, it
    % should take in account different points for linearizing
    
    vt_curPhase     = vt_uwPhase - vt_uwPhase(nm_curPhSample) ...
                    + nm_curPhValue;
                                
    vt_trialSamples = find(...
                    st_cnf.window(1) - st_cnf.border <= vt_curPhase & ...
                    vt_curPhase <= st_cnf.window(2) + st_cnf.border);
                
    vt_trialSegment	= vt_signal(vt_trialSamples);
    vt_trialPhase   = vt_curPhase(vt_trialSamples);
    
    % Phase boucle correction
    vt_phBoucle     = findextrema(vt_trialPhase);
    
    if ~isempty(vt_phBoucle)
        
        for  tt = 1:numel(vt_phBoucle)
            nm_start	= vt_phBoucle(tt);
            nm_end      = find(vt_trialPhase > vt_trialPhase(vt_phBoucle(tt)),...
                        1,'first');
                    
            if isempty(nm_end)
                nm_end = numel(vt_trialPhase);
            end
            
            nm_Id   = nm_start:nm_end;
            vt_x    = [nm_start,nm_end];
            vt_y    = [nm_start,nm_end];
            vt_line	= polyfit(vt_x(:),vt_y(:),1);
            
            vt_trialPhase(nm_Id)	= vt_line(1)*nm_Id+vt_line(2);
        end % Boucle correction for each found
    end % if phase boucles are present in the current trial
        
    % Convert phase trial into phase bins
    vt_cSegm	= zeros(size(vt_phBins));
    
    for tt = 1:numel(vt_phBins)
        vt_idx   = vt_trialPhase >= vt_phBins(tt) & ...
                vt_trialPhase < vt_phBins(tt) + nm_phRes;
        
        if sum(vt_idx) == 0
            if tt == 1
                vt_idx	= find(vt_trialPhase < vt_phBins(tt),1,'last');
            elseif tt == numel(vt_phBins)
                vt_idx	= find(vt_trialPhase > vt_phBins(tt),1,'first');
            end
        end
        
        vt_cSegm(tt)	= mean(vt_trialSegment(vt_idx));
                
    end
    
    % Find empty bins and make linear interpolation depending on
    % neighborhoods
    vt_idx	= isnan(vt_cSegm);
    
    if any(vt_idx)
        vt_line     = find(vt_idx==0);        
        vt_cSegm	= interp1(...
                    vt_phBins(vt_line),vt_cSegm(vt_line),vt_phBins);        
    end
    
    vt_trials{kk}   = vt_cSegm(:)';
end

