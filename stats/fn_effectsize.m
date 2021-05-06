function [nm_d,vt_ci] = fn_effectsize(vt_x,vt_y)
% [nm_d,vt_ci] = fn_effectsize(vt_x,vt_y)
% Function for further development. Check:
% Nakagawa & Cuthil. Effect size, confidence interval and statistical 
% significance: A practical guide for biologists. Biological Reviews. 2007

% This case is for comparing two independent or dependent groups. (i.e.
% both  paired and unpaired t-test cases)

vt_x	= vt_x(:);
vt_y	= vt_y(:);
vt_x	= vt_x(~isnan(vt_x));
vt_y	= vt_y(~isnan(vt_y));

nm_mx   = mean(vt_x);
nm_my   = mean(vt_y);
nm_sx   = std(vt_x);
nm_sy   = std(vt_y);
nm_nx   = numel(vt_x);
nm_ny   = numel(vt_y);

nm_spooled	=  sqrt(...
            ((nm_ny - 1)*nm_sy^2 + (nm_nx - 1)*nm_sx^2)... 
            /(nm_nx + nm_ny - 2));
        
nm_d 	= (nm_my - nm_mx) / nm_spooled;

% For d(independent,unpaired)
nm_se	= sqrt(...
        ((nm_nx + nm_ny - 1)/(nm_nx + nm_ny - 3))...
        *((4/(nm_nx + nm_ny))*(1 + nm_d^2/8)));
    
vt_ci   = [(nm_d - 1.96*nm_se), (nm_d + 1.96*nm_se)];
