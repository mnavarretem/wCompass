function [vt_aTrl,vt_phTrl] = fn_trials_phaseamplitude(vt_ampl,vt_phase,st_cnf)
% [vt_aTrl,vt_phTrl] = fn_trials_phaseamplitude(vt_ampl,vt_phase,st_cnf) 
% gets the single trial phase-amplitude lineseries. vt_ampl values are 
% referenced to phases of vt_phase timeseries. Trials settings are 
% configured in the st_cnf structure. vt_aTrl is a cell vector composed 
% by the selected trials
% 
% Use as
% vt_aTrl = fn_trials_phaseamplitude(vt_ampl,vt_phase,st_cnf) 
% [vt_aTrl,vt_phTrl] = fn_trials_phaseamplitude(vt_ampl,vt_phase,st_cnf) 
%
% st_cnf corresponds to a structure with the fields:
%
%   - st_cnf.markers:       Vector indicating the sample markers defining
%                           each trial <triggers>
%   - st_cnf.window:        1 x 2 vector indicating the begining and end
%                           phase (in radians) of each trial surrounding
%                           each marker position
%   - st_cnf.align:         Array indicating the samples in wich the phase
%                           of the trial should be aligned
%   - st_cnf.binarize:      Convert phase timeseries to timebins (default:
%                           true)
%   - st_cnf.nbins:         Number of bins for computing the phase
%                           distribution (default: 100)
%   - st_cnf.border:        Phase 'lag' to include around the phase window
%   - st_cnf.correction:	Phase loop correction ('polyfix', fix linear
%                           phase to the best fit in the LMS sense from
%                           polymials from order 1 to 5. 'flat', flattens
%                           phase during negative slope sections. 'linear',
%                           linear fitting during loop. 'none', no
%                           correction. Default: 'polyfix')
%                           

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
if nargin < 3   
    error('fn_trials_phaseamplitude:wrongInputs','function requires three inputs'); 
end

% Check whether st_cnf fileds need defaults
if ~isfield(st_cnf,'align') 
    st_cnf.align	= st_cnf.markers;
end
if ~isfield(st_cnf,'correction') 
    st_cnf.correction	= 'polyfix';
end
if ~isfield(st_cnf,'binarize') 
    st_cnf.binarize     = true;
end
if ~isfield(st_cnf,'nbins') 
    st_cnf.nbins        = 100;
end
if ~isfield(st_cnf,'border') 
    st_cnf.border       = pi/6;
end

vt_ORDER    = 1:5;

%%	- Extract trials
vt_ampl     = vt_ampl(:)';
vt_uwPhase	= unwrap(vt_phase);

vt_aTrl     = cell(numel(st_cnf.markers),1);
vt_phTrl	= cell(numel(st_cnf.markers),1);
 
for kk = 1:numel(st_cnf.markers)
    
    nm_curPhSample	= unique(st_cnf.align(kk,:));
    nm_curPhValue   = vt_phase(nm_curPhSample);
    
    % NOTE: Only uses one value without linearizing multiple markers, it
    % should take in account different points for linearizing
    
    vt_curPhase     = vt_uwPhase - vt_uwPhase(nm_curPhSample) ...
                    + nm_curPhValue;
                                
    vt_aTrlSamples  = find(...
                    st_cnf.window(1) - st_cnf.border <= vt_curPhase & ...
                    vt_curPhase <= st_cnf.window(2) + st_cnf.border);
                
    vt_trialAmp     = vt_ampl(vt_aTrlSamples);
    vt_trialPh      = vt_curPhase(vt_aTrlSamples);
    
    % Phase boucle correction
    vt_phBoucle     = findextrema(vt_trialPh);
    
    if ~isempty(vt_phBoucle)  
        switch st_cnf.correction
            case 'polyfix'                
                vt_err	= zeros(size(vt_ORDER));
                vt_fit	= cell(size(vt_ORDER));
                vt_nTrl = 1:numel(vt_trialPh);
                
                for od = vt_ORDER                    
                    vt_polyFit	= polyfix(vt_nTrl(:),vt_trialPh(:),od,...
                                [1,numel(vt_trialPh)],vt_trialPh([1,end]));
                    vt_fitPh	= polyval(vt_polyFit,vt_nTrl);
                    
                    if any(diff(vt_fitPh)<0)
                        vt_err(od)	= inf;
                    else
                        vt_err(od)	= sum((vt_trialPh(:)-vt_fitPh(:)).^2);
                    end
                    vt_fit{od}	= vt_fitPh;
                end
                
                [~,od]   	= min(vt_err);
                vt_trialPh	= vt_fit{od};

            case 'flat'
                for  tt = 1:numel(vt_phBoucle)
                    nm_beg	= vt_phBoucle(tt);
                    nm_end	= find(vt_trialPh > vt_trialPh(vt_phBoucle(tt)),...
                            1,'first');

                    if isempty(nm_end)
                        nm_end = numel(vt_trialPh);
                    end

                    nm_Id   = nm_beg:nm_end;

                    vt_trialPh(nm_Id)	= vt_trialPh(nm_beg);
                end % tt
                
            case 'linear'
                for  tt = 1:numel(vt_phBoucle)
                    [vt_hi,vt_lo]	= findextrema(vt_trialPh);
                    vt_hi   = vt_hi(vt_phBoucle(tt) == vt_hi);
                    vt_lo   = vt_lo(vt_lo > vt_phBoucle(tt));
                    vt_lo   = vt_lo(1);
                    vt_hiPh	= vt_trialPh(vt_hi);
                    vt_loPh	= vt_trialPh(vt_lo); 
                    
                    nm_beg	= find(vt_trialPh < vt_loPh,1,'last');
                    nm_end	= find(vt_trialPh > vt_hiPh,1,'first');

                    if isempty(nm_end)
                        nm_end = numel(vt_trialPh);
                    end

                    nm_Id   = nm_beg:nm_end;
                    vt_x    = [nm_beg,nm_end];
                    vt_y    = vt_trialPh([nm_beg,nm_end]);
                    vt_line	= polyfit(vt_x(:),vt_y(:),1);

                    vt_trialPh(nm_Id)	= vt_line(1)*nm_Id+vt_line(2);
                end % tt
                
            otherwise
                % Do nothing
                
        end % case
    end % if phase boucles are present in the current trial
    
    if st_cnf.binarize
        vt_phBins	= linspace(...
                    st_cnf.window(1),st_cnf.window(2),st_cnf.nbins + 1);

        % Convert phase trial into phase bins
        vt_binPh	= discretize(vt_trialPh,vt_phBins);
        vt_mAmp     = vt_trialAmp(~isnan(vt_binPh)); 
        vt_mBins	= vt_binPh(~isnan(vt_binPh)); 
        
        vt_cAmpl    = accumarray(vt_mBins(:),vt_mAmp(:),...
                    [st_cnf.nbins 1],@mean);
                
        % Find empty bins and make linear interpolation depending on
        % neighborhoods
        vt_idx	= vt_cAmpl == 0 | isnan(vt_cAmpl);

        if any(vt_idx)
            vt_line     = find(vt_idx == 0);        
            vt_cAmpl	= interp1(...
                        vt_phBins(vt_line),vt_cAmpl(vt_line),vt_phBins,...
                        'linear','extrap');         
            vt_cAmpl    = vt_cAmpl(1:end-1);
        end
        vt_phBins       = fn_edge2bin(vt_phBins);
        vt_aTrl{kk}     = vt_cAmpl(:)';
        vt_phTrl{kk}    = vt_phBins(:)';
    else          
        vt_id	= find(st_cnf.window(1) <= vt_trialPh & ...
                vt_trialPh <= st_cnf.window(2));

        vt_trialAmp     = vt_trialAmp(vt_id);
        vt_trialPh      = vt_trialPh(vt_id);
        
        vt_aTrl{kk}     = vt_trialAmp(:)';
        vt_phTrl{kk}    = vt_trialPh(:)';
    end    
end

