function mx_mat	= fn_confusionmatrix(vt_true,vt_pred)
%fn_confusionmatrix Confusion matrix of logical inputs
%   mx_mat	= fn_confusionmatrix(vt_true,vt_pred) takes logical vectors,
%   and returns the Confusion matrix. 
%   
%    Each row of the matrix represents the instances in a predicted class 
% while each column represents the instances in an actual class 
% (or vice versa)

vt_true = logical(vt_true);
vt_pred = logical(vt_pred);

vt_idTP	= vt_true == 1 & vt_pred == 1;
vt_idFP	= vt_true == 0 & vt_pred == 1;
vt_idFN	= vt_true == 1 & vt_pred == 0;
vt_idTN	= vt_true == 0 & vt_pred == 0;

mx_mat(1,1) = sum(vt_idTP);
mx_mat(1,2) = sum(vt_idFP);
mx_mat(2,1) = sum(vt_idFN);
mx_mat(2,2) = sum(vt_idTN);

