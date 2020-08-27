% function f_Matrix2ComplexNorm.m
% 
function m_NormMat = ...
    f_Matrix2ComplexNorm(...
    pm_InMatrix,ps_IsRobust)

    if nargin < 1
        error('[f_Matrix2ComplexNorm] - ERROR: bad number of parameters!')
    elseif nargin < 2
        ps_IsRobust	= false;
    end

    if ~iscomplex(pm_InMatrix)
        error('[f_Matrix2ComplexNorm] - ERROR: input matrix must be complex!')
    end
    if ps_IsRobust
        m_RealMat   = f_RobustMatrixZScore(real(pm_InMatrix));
        m_ImagMat   = f_RobustMatrixZScore(imag(pm_InMatrix));
    else        
        m_RealMat	= real(pm_InMatrix);
        m_ImagMat	= imag(pm_InMatrix);
        
        m_MAD       = repmat(mad(m_RealMat',0)',1,size(m_RealMat,2));
        m_Sigma     = 1.4826*m_MAD;
        
        m_RealMat   = (m_RealMat - m_MAD)./m_Sigma;
        
        m_MAD       = repmat(mad(m_ImagMat',0)',1,size(m_ImagMat,2));
        m_Sigma     = 1.253*m_MAD;
%         m_Sigma     = 1.4826*m_MAD;
        
        m_ImagMat   = (m_ImagMat - m_MAD)./m_Sigma;
%         
%         s_Numel     = size(m_RealMat,2);
%         v_QtlPos	= round([0.25 0.75] * s_Numel);
%         
%         m_Val	= sort(m_RealMat,2);
%         m_Qtl   = m_Val(:,v_QtlPos);
%         m_Rng   = 1.5*diff(m_Qtl,1,2);
%         
%         m_QtlLo = repmat(m_Qtl(:,1),1,s_Numel) - repmat(m_Rng(:,1),1,s_Numel);
%         m_QtlHi = repmat(m_Qtl(:,2),1,s_Numel) + repmat(m_Rng(:,1),1,s_Numel);
%         m_Idx   = m_Val >= m_QtlLo & m_Val <= m_QtlHi;
%         m_Avg	= repmat(sum(m_Val.*m_Idx,2) ./ sum(m_Idx,2),1,s_Numel);
%         m_Std   = (m_Val-m_Avg).^2;
%         m_Std	= sqrt(repmat(sum(m_Std.*m_Idx,2) ./ sum(m_Idx,2),1,s_Numel));
%         
%         m_ReAvg = m_Avg;
%         m_ReStd = m_Std;
%         
%         m_Val	= sort(m_ImagMat,2);
%         m_Qtl   = m_Val(:,v_QtlPos);
%         m_Rng   = diff(m_Qtl,1,2);
%         
%         m_QtlLo = repmat(m_Qtl(:,1),1,s_Numel) - repmat(m_Rng(:,1),1,s_Numel);
%         m_QtlHi = repmat(m_Qtl(:,2),1,s_Numel) + repmat(m_Rng(:,1),1,s_Numel);
%         m_Idx   = m_Val >= m_QtlLo & m_Val <= m_QtlHi;
%         m_Avg	= repmat(sum(m_Val.*m_Idx,2) ./ sum(m_Idx,2),1,s_Numel);
%         m_Std   = (m_Val-m_Avg).^2;
%         m_Std	= repmat(sum(m_Std.*m_Idx,2) ./ sum(m_Idx,2),1,s_Numel);
%         
%         m_ImAvg = m_Avg;
%         m_ImStd = m_Std;
%         
%         m_RealMat = (m_RealMat - m_ReAvg)./m_ReStd;
%         m_ImagMat = (m_ImagMat - m_ImAvg)./m_ImStd;
    end
    m_NormMat = sqrt(m_RealMat.^2+m_ImagMat.^2);
end
