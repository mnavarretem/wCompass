function mx_outNorm	= fn_stimuluslock_trialnormalization(mx_trials,st_cfg)
                
%% Read variables
if ~isfield(st_cfg,'dimension')
    if numel(size(mx_trials)) == 3
        st_cfg.dimension	= 3;
    else
        st_cfg.dimension	= 1;
    end
end

if iscell(mx_trials)
    mx_trials	= cell2mat(mx_trials);
end

%% process trials

switch st_cfg.approach
    case 'single'
        vt_bLine	= [];
        mx_outNorm	= fn_baselinenormalization(mx_trials,...
                    vt_bLine,'abs');
        mx_outNorm	= mean(mx_outNorm,st_cfg.dimension,'omitnan');
        mx_outNorm	= fn_baselinenormalization(mx_outNorm,...
                    st_cfg.baseline,st_cfg.method);
    case 'classical'
        mx_outNorm	= mean(mx_trials,st_cfg.dimension,'omitnan');
        mx_outNorm	= fn_baselinenormalization(mx_outNorm,...
                    st_cfg.baseline,st_cfg.method);
end