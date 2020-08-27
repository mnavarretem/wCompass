function [nm_pac,nm_angle,st_stat]	= fn_trials_pac(vt_ampl,vt_phase,st_cnf)
%  nm_pac = fn_trials_pac(vt_amplitude,vt_phase,st_cnf) calculates the 
% trial-based phase-amplitude coupling (PAC). vt_ampl is the amplitude 
% evolving timeseries, vt_phase is the phase evolving timeseries and trials 
% settings are configured in the st_cnf structure. 
%
% Use as:
%  nm_pac	= fn_trials_pac(vt_ampl,vt_phase,st_cnf); 
%  [nm_pac,nm_angle]	= fn_trials_pac(vt_ampl,vt_phase,st_cnf);
%  [nm_pac,nm_angle,st_stat]	= fn_trials_pac(vt_ampl,vt_phase,st_cnf);
% 
% OUTPUT: 
%   nm_pac:     Coupling strength
%   nm_angle:   Preferred coupling angle
%   st_stats:   Coupling statistics
% 
% INPUT
%  vt_ampl:     Time evolving amplitude vector
%  vt_phase:    Time evolving phase vector
%  st_cnf:      Structure with fields:
%
%   - st_cnf.markers:       Vector indicating the sample markers defining
%                           each trial <triggers>
%   - st_cnf.window:        1 x 2 vector indicating the begining and end
%                           phase (in radians) of each trial surrounding
%                           each marker position (default: [-pi,pi])
%   - st_cnf.method:        Method for computing PAC: 
%                           'direct': Averaging over trial samples (Özkurt et al 2011),     
%                           'mvl': Mean vector length (Canolty et al. 2006),      
%                           'umvl': Unbiased mean vector length (Kutil 2012),
%                           'kl': Kullback–Leibler divergence (Tort 2008),
%                           'glm': General linear model (Penny 2008),
%                           (default: 'umvl')
%   - st_cnf.normalize:     Normalize amplitude (default: true)
%   - st_cnf.binarize:      Convert phase timeseries to timebins (default:
%                           true)
%   - st_cnf.nbins:         Number of bins for computing the phase
%                           distribution (default: 100)
%   - st_cnf.surrugates:    Compute data surrogates (dafault: true)
%   - st_cnf.nSurrugates:	Number of permutations to calculate surrugates
%                           (default: 200)
%   - st_cnf.rOutliers:     Remove outliers (default: true)
%   - st_cnf.correction:	Phase loop correction ('polyfix', fix linear
%                           phase to the best fit in the LMS sense from
%                           polymials from order 1 to 5. 'flat', flattens
%                           phase during negative slope sections. 'linear',
%                           linear fitting during loop. 'none', no
%                           correction - default: 'polyfix')
%   - st_cnf.border:        Phase 'lag' to include around the phase window
%                           (default: pi/6)
% 
% Example:
% % Simulate coupled and uncupled data
% fs = 250;t = 0:1/fs:900;
% lowDat  = (5*sin(2*pi*0.8*t)).*(0.5*sin(2*pi*1e-4*t)+2);
% highDat1 = 0.25*sin(2*pi*15*t).*((sin(2*pi*0.8*t))+1).*(0.15*sin(2*pi*1e-3*t)+2);
% highDat2 = 0.25*sin(2*pi*15*t).*((sin(2*pi*0.5*t))+1).*(0.15*sin(2*pi*1e-3*t)+2);
% highDat3 = 0.25*sin(2*pi*15*t).*((sin(2*pi*3.5*t))+1).*(0.15*sin(2*pi*1e-2*t)+2);
% plot(t,lowDat + highDat1,t,lowDat + highDat2,t,lowDat + highDat3)
% lowPhase = angle(hilbert(lowDat));
% highAmp1 = abs(hilbert(highDat1));
% highAmp2 = abs(hilbert(highDat2));
% highAmp3 = abs(hilbert(highDat3));
% markers  = sort(round(rand(120,1)*numel(lowPhase)));
% markers(diff(markers) < 2 * fs) = [];
% 
% % Try default values
% st_cnf = [];
% st_cnf.markers = markers;
% [nm_c1,nm_a1,st_st1]	= fn_trials_pac(highAmp1,lowPhase,st_cnf);
% [nm_c2,nm_a2,st_st2]	= fn_trials_pac(highAmp2,lowPhase,st_cnf);
% [nm_c3,nm_a3,st_st3]	= fn_trials_pac(highAmp3,lowPhase,st_cnf);
% 
% % Try different methods
% st_cnf = [];
% st_cnf.markers = markers;
% st_cnf.method	= 'kl';
% [nm_c1,nm_a1,st_st1]	= fn_trials_pac(highAmp1,lowPhase,st_cnf);
% [nm_c2,nm_a2,st_st2]	= fn_trials_pac(highAmp2,lowPhase,st_cnf);
% [nm_c3,nm_a3,st_st3]	= fn_trials_pac(highAmp3,lowPhase,st_cnf);
% 
% % Try different methods
% st_cnf = [];
% st_cnf.markers = markers;
% st_cnf.method	= 'glm';
% [nm_c1,nm_a1,st_st1]	= fn_trials_pac(highAmp1,lowPhase,st_cnf);
% [nm_c2,nm_a2,st_st2]	= fn_trials_pac(highAmp2,lowPhase,st_cnf);
% [nm_c3,nm_a3,st_st3]	= fn_trials_pac(highAmp3,lowPhase,st_cnf);
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
    error('fn_trials_pac:wrongInputs','function requires three inputs'); 
end

% Check whether st_cnf fileds need defaults
if ~isfield(st_cnf,'window') 
    st_cnf.window       = [-pi,pi];
end
if ~isfield(st_cnf,'method') 
    st_cnf.method       = 'umvl';
end
if ~isfield(st_cnf,'normalize') 
    st_cnf.normalize	= true;
end
if ~isfield(st_cnf,'binarize') 
    st_cnf.binarize  	= true;
end
if ~isfield(st_cnf,'nbins') 
    st_cnf.nbins        = 100;
end
if ~isfield(st_cnf,'surrugates') 
    st_cnf.surrugates	= true;
end
if ~isfield(st_cnf,'nSurrugates') 
    st_cnf.nSurrugates	= 200;
end
if ~isfield(st_cnf,'rOutliers') 
    st_cnf.rOutliers	= true;
end
if ~isfield(st_cnf,'correction') 
    st_cnf.correction	= 'polyfix';
end

%% Pre-process trials
% Obtain phase-amplitude trials
[vt_aTrl,vt_phTrl] = fn_trials_phaseamplitude(vt_ampl,vt_phase,st_cnf);

% Remove outliers
if st_cnf.rOutliers
    vt_dat      = cell2mat(vt_aTrl(:)');
    nm_thres	= prctile(abs(vt_dat(:)),99.8);
    vt_remove   = cellfun(@(x) any(abs(x) > nm_thres | x < 0) ,vt_aTrl);
    vt_aTrl     = vt_aTrl(~vt_remove);
    vt_phTrl 	= vt_phTrl(~vt_remove);
end

if st_cnf.normalize
    vt_aTrl	= cellfun(@(x) x ./ mean(x), vt_aTrl,'UniformOutput',false);    
end

%% Compute PAC

vt_ampli	= cell2mat(vt_aTrl(:)');
vt_phase	= cell2mat(vt_phTrl(:)');

switch st_cnf.method
    case 'direct'
        [nm_pac,nm_angle]	= fn_local_direct(vt_ampli,vt_phase);
    case 'mvl'
        [nm_pac,nm_angle]   = fn_local_mvl(vt_ampli,vt_phase);
    case 'umvl'
        [nm_pac,nm_angle]   = fn_local_umvl(vt_ampli,vt_phase);
    case 'kl'
        [nm_pac,nm_angle]   = fn_local_kl(vt_ampli,vt_phase);
    case 'glm'
        [nm_pac,nm_angle]   = fn_local_glm(vt_ampli,vt_phase);
end

%% Compute surrogates
if st_cnf.surrugates
    % Calculate the shift surrogate for each trial
    vt_numel    = cellfun(@numel,vt_aTrl);
    mx_rSteps	= rand(numel(vt_aTrl),st_cnf.nSurrugates); 
    
    vt_cSurr    = nan(st_cnf.nSurrugates,1);
    vt_aSurr    = nan(st_cnf.nSurrugates,1);
    
    % Loop surrugates
    for ss = 1:st_cnf.nSurrugates
        vt_shift    = mx_rSteps(:,ss);
        vt_shift    = round(vt_numel .* vt_shift);
        vt_shift    = num2cell(vt_shift);
        
        vt_aShift   = cellfun(@(x,y) circshift(x,y,2),...
                    vt_aTrl,vt_shift,'UniformOutput',false);
        
        vt_ampli	= cell2mat(vt_aShift(:)');
        vt_phase	= cell2mat(vt_phTrl(:)');

        switch st_cnf.method
            case 'direct'
                [nm_cSurr,nm_aSurr]	= fn_local_direct(vt_ampli,vt_phase);
            case 'mvl'
                [nm_cSurr,nm_aSurr]	= fn_local_mvl(vt_ampli,vt_phase);
            case 'umvl'
                [nm_cSurr,nm_aSurr]	= fn_local_umvl(vt_ampli,vt_phase);
            case 'kl'
                [nm_cSurr,nm_aSurr]	= fn_local_kl(vt_ampli,vt_phase);
            case 'glm'
                [nm_cSurr,nm_aSurr]	= fn_local_glm(vt_ampli,vt_phase);
        end
        
        vt_cSurr(ss)	= nm_cSurr;
        vt_aSurr(ss)    = nm_aSurr;
    end %ss
    
    st_stat.surrugates.pac      = vt_cSurr;
    st_stat.surrugates.angle	= vt_aSurr;
    
    % Compute statistics for couplig
    st_stat.pac.zScore	= (nm_pac - mean(vt_cSurr)) / std(vt_cSurr);
    
    st_stat.pac.pValue      = 1 - sum(nm_pac >= vt_cSurr) /...
                            st_cnf.nSurrugates;
    st_stat.pac.diffMean	= mean(nm_pac - vt_cSurr);    
    st_stat.pac.diff95CI    = [prctile(nm_pac - vt_cSurr,5),...
                            prctile(nm_pac - vt_cSurr,95)];   
                        
else
    % If surrogated are not computed, then st_stat.pac is an empty structure
    st_stat.pac	= struct;
    
    st_stat.surrugates.pac      = nan;
    st_stat.surrugates.angle	= nan;
end

% Compute angular statistics if circ_stat exist in MATLB path
nm_isCirc	= exist('circ_stats.m','file');
if nm_isCirc == 2
    vt_trialAngle       = cellfun(@circ_mean,vt_phTrl,vt_aTrl);
    [nm_a,nm_aH,vt_aL]  = circ_mean(vt_trialAngle);
    
	% Rayleigh test for non-uniformity of circular data
	[nm_raylP,nm_raylZ] = circ_rtest(vt_trialAngle);
    
    if ~isnan(st_stat.surrugates.angle)
        % Watson-Williams multi-sample test for equal means
        [nm_wwTest,mx_tab]	= circ_wwtest(vt_trialAngle,st_stat.surrugates.angle); 
    else
        nm_wwTest   = nan;
    end
    
    st_stat.angle.trials	= vt_trialAngle;
    st_stat.angle.mean      = nm_a;
    st_stat.angle.CI        = [vt_aL,nm_aH];
    
    st_stat.angle.wwTest.p      = nm_wwTest;
    st_stat.angle.wwTest.table	= mx_tab;
    
    st_stat.angle.rayleigh.p	= nm_raylP;
    st_stat.angle.rayleigh.z    = nm_raylZ;
else
    % If circ_stat does not exist, then st_stat.angle is an empty structure
    st_stat.angle	= struct;
end

%% Nested functions :::::::::::::::::::::::::::::::::::::::::::::::::::::::
    function [nm_c,nm_a]= fn_local_direct(vt_am,vt_ph)
        % (Özkurt et al 2011) divides by the power of amplitude vector
        nm_num	= numel(vt_am);
        nm_c	= abs(sum(vt_am .* exp(1i*vt_ph)));
        nm_c    = (1 / sqrt(nm_num)) * (nm_c ./ sqrt(sum(vt_am.^2)));
        % Compute preferred angle
        nm_a    = angle(sum(vt_am .* exp(1i*vt_ph)));
    end
%--------------------------------------------------------------------------
    function [nm_c,nm_a]= fn_local_mvl(vt_am,vt_ph)
        % (Canolty et al. 2006)
        nm_c	= sum(vt_am .* exp(1i*vt_ph));
        nm_c  	= abs(nm_c) ./ sum(vt_am);
        % Compute preferred angle
        nm_a    = angle(sum(vt_am .* exp(1i*vt_ph)));
        
    end
%--------------------------------------------------------------------------
    function [nm_c,nm_a]= fn_local_umvl(vt_am,vt_ph)
        % (Kutil 2012)
        vt_mpac	= (vt_am .* exp(1i*vt_ph)).^2;
        nm_c	= abs(sum(vt_mpac)) ./ sum(abs(vt_mpac));
        % Compute preferred angle
        nm_a    = angle(sum(vt_am .* exp(1i*vt_ph)));
        
    end
%--------------------------------------------------------------------------
    function [nm_c,nm_a]= fn_local_kl(vt_am,vt_ph)
        % (Tort 2008)
        % Convert phase trial into phase bins
        nm_bins	= 36; % start with 10 degrees bins
        
        while true
            vt_phBins	= linspace(...
                        st_cnf.window(1),st_cnf.window(2),nm_bins + 1);
            vt_binPh	= discretize(vt_ph,vt_phBins);
            vt_mAmp     = vt_am(~isnan(vt_binPh)); 
            vt_mBins	= vt_binPh(~isnan(vt_binPh)); 

            vt_ampPdf	= accumarray(vt_mBins(:),vt_mAmp(:),...
                        [nm_bins,1],@mean);
                    
            % Check for empty bins, KL canot be computed because of log(0) 
            if any(vt_ampPdf == 0)
               %  decrease bin number and continue
               nm_bins	= nm_bins - 1;
            else
                % compute pdf
                vt_ampPdf	= vt_ampPdf ./ sum(vt_ampPdf);
                break
            end
        end
        
        % Compute KL distance
        nm_Dkl 	= log(numel(vt_ampPdf)) + sum(vt_ampPdf .* log(vt_ampPdf));
        nm_c	= nm_Dkl / log(numel(vt_ampPdf));
        
        % Compute preferred angle
        [~,nm_id]	= max(vt_ampPdf);
        nm_a        = mean(vt_phBins([nm_id,nm_id+1]));
    end
%--------------------------------------------------------------------------
    function [nm_c,nm_a]= fn_local_glm(vt_am,vt_ph)
        % (Penny 2008)
        % Build regressors
        vt_am   = vt_am(:);
        vt_ph   = vt_ph(:);
        mx_dat	= [cos(vt_ph), sin(vt_ph) ones(size(vt_ph))];
        % Fit glm model 
        [vt_b,~,st_st]	= glmfit(mx_dat,vt_am,'normal','constant','off');
        
        nm_c	= 1 - sum(st_st.resid.^2) / sum((vt_am - mean(vt_am)).^2);  
        
        % Compute preferred angle
        % Set 360 phase model values        
        vt_phBins	= linspace(st_cnf.window(1),st_cnf.window(2),360);
        vt_phBins   = vt_phBins(:);
        mx_dat      = [cos(vt_phBins), sin(vt_phBins) ones(size(vt_phBins))];
        % Predict amplitude values. 
        % Link: 'identity' (as for 'normal' distribution)
        vt_ampPred  = glmval(vt_b,mx_dat,'identity','constant','off'); 
        
        % Get preferred angle
        [~,nm_id]	= max(vt_ampPred);
        nm_a        = vt_phBins(nm_id);
    end
end