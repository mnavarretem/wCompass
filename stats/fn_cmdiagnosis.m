function st_stat = fn_cmdiagnosis(vt_true,vt_pred)
%cfdiagnosis Binary classification diagnosis based on confusion matrix
%   st_stat	= cfdiagnosis(vt_true,vt_pred) takes logical vectors,
%   and returns the statistics from a confusion matrix. 
%   
%    Each row of the matrix represents the instances in a predicted class 
% while each column represents the instances in an actual class 
% (or vice versa)
%
% Output:
%
%	st_stat.TPR: Recall, sensitivity or true positive rate (TPR)
%   st_stat.TNR: Specificity, selectivity or true negative rate (TNR)
%   st_stat.PPV: precision or positive predictive value (PPV)
%   st_stat.NPV: negative predictive value (NPV)
%   st_stat.FNR: miss rate or false negative rate (FNR)
%   st_stat.FPR: fall-out or false positive rate (FPR)
%   st_stat.FDR: false discovery rate (FDR)
%   st_stat.FOR: false omission rate (FOR)
%   st_stat.TS:  Threat score (TS) or Critical Success Index (CSI)
%   st_stat.ACC: accuracy (ACC)
%   st_stat.BA:  balanced accuracy (BA)
%   st_stat.F1:  F1 score
%   st_stat.MCC: Matthews correlation coefficient (MCC)
%   st_stat.BM:  informedness or bookmaker informedness (BM)
%   st_stat.MK:  markedness (MK) or deltaP
%   st_stat.MAE: mean absolute error (MAE)

%% code starts here
% check size
if numel(vt_true) ~= numel(vt_pred)
    error('[fn_cmdiagnosis] both x and y vectors should have the same numel')
end

%Discard NaN;s
vt_id	= ~(isnan(vt_true) | isnan(vt_pred));
vt_true    = vt_true(vt_id);
vt_pred    = vt_pred(vt_id);

mx_cm	= fn_confusionmatrix(vt_true,vt_pred);

nm_TP	= mx_cm(1,1);
nm_FP	= mx_cm(1,2);
nm_FN	= mx_cm(2,1);
nm_TN	= mx_cm(2,2);

% Recall, sensitivity or true positive rate (TPR)
st_stat.TPR = nm_TP/(nm_TP+nm_FN);
% Specificity, selectivity or true negative rate (TNR)
st_stat.TNR = nm_TN/(nm_TN+nm_FP);
% precision or positive predictive value (PPV)
st_stat.PPV	= nm_TP/(nm_TP+nm_FP);
% negative predictive value (NPV)
st_stat.NPV	= nm_TN/(nm_TN+nm_FN);
% miss rate or false negative rate (FNR)
st_stat.FNR	= 1 - st_stat.TPR;
% fall-out or false positive rate (FPR)
st_stat.FPR	= 1 - st_stat.TNR;
% false discovery rate (FDR)
st_stat.FDR	= 1 - st_stat.PPV;
% false omission rate (FOR)
st_stat.FOR	= 1 - st_stat.PPV;
% Threat score (TS) or Critical Success Index (CSI)
st_stat.TS  = nm_TP/(nm_TP+nm_FN+nm_FP);
% accuracy (ACC)
st_stat.ACC	= (nm_TP + nm_TN)/(nm_TP + nm_TN + nm_FP + nm_FN);
% balanced accuracy (BA)
st_stat.BA  = (st_stat.TPR + st_stat.TNR) / 2;
% F1 score
st_stat.F1  = 2 * (st_stat.PPV * st_stat.TPR)/(st_stat.PPV + st_stat.TPR);
% Matthews correlation coefficient (MCC)
st_stat.MCC	= (nm_TP * nm_TN - nm_FP * nm_FN) / sqrt((nm_TP + nm_FP) * ...
            (nm_TP + nm_FN) * (nm_TN + nm_FP) * (nm_TN + nm_FN));
% informedness or bookmaker informedness (BM)
st_stat.BM  = st_stat.TPR + st_stat.TNR - 1;
% markedness (MK) or deltaP
st_stat.MK  = st_stat.PPV + st_stat.NPV - 1;
% mean absolute error (MAE)
st_stat.MAE	= mean(abs(vt_true-vt_pred));
