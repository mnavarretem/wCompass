% function f_Matrix2Gain.m
% 
function mx_gain = ...
    f_Matrix2Gain(...
    pm_InMatrix, ...
    ps_FirstIndRef, ...
    ps_LastIndRef,...
    ps_IsLog)

    if nargin < 1
        error('[f_Matrix2ZScore] - ERROR: bad number of parameters!')
    end

    mx_gain = pm_InMatrix;
    
    if ~exist('ps_FirstIndRef', 'var') || isempty(ps_FirstIndRef)
        ps_FirstIndRef = 1;
    end
    if ~exist('ps_LastIndRef', 'var') || isempty(ps_LastIndRef)
        ps_LastIndRef = size(mx_gain, 2);
    end
    if ~exist('ps_IsLog', 'var') || isempty(ps_IsLog)
        ps_IsLog = false;
    end
    
    if ps_FirstIndRef < 1
        ps_FirstIndRef = 1;
    end
    if ps_LastIndRef > size(mx_gain, 2)
        ps_LastIndRef = size(mx_gain, 2);
    end

    clear v_Mean v_Std
    if numel(size(pm_InMatrix)) == 2

        vt_mean = mean(mx_gain(:, ps_FirstIndRef:ps_LastIndRef), 2);
        mx_gain	= mx_gain ./ repmat(vt_mean, 1, size(mx_gain, 2));
        
    elseif numel(size(pm_InMatrix)) == 3

        vt_mean = mean(mx_gain(:,ps_FirstIndRef:ps_LastIndRef,:),2);
        mx_gain	= mx_gain ./ repmat(vt_mean,1,size(mx_gain, 2),1);
        
    end
    if ps_IsLog
        mx_gain = 10 * log10(mx_gain);
    end
end
