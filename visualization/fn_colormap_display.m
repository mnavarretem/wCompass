function ob_f	= fn_colormap_display(mx_map,vt_lim,ch_label)
% ob_f	= fn_colormap_display(mx_map,vt_lim) displays mx_map as a colorbar
% figure within limits from vt_lim = [min,max] (Default: vt_lim = [0,1]),
% and with colorbar label ch_label (Default: '').
% Output ob_f is the figure handle
% 
% Use as:
%	ob_f	= fn_colormap_display(mx_map)
%   ob_f	= fn_colormap_display(mx_map,vt_lim)
%   ob_f	= fn_colormap_display(mx_map,vt_lim,ch_label)
%       
% <2020> Miguel Navarrete

if nargin < 2
    vt_lim = [0,1];
end

if nargin < 3
    ch_label = '';
end

vt_x    = 1;
vt_y    = linspace(vt_lim(1),vt_lim(2),size(mx_map,1));
vt_z    = vt_y;

ob_f    = figure(...
        'Units','pixels',...
        'Position',[100,100,100,250]);

ob_ax   = axes(ob_f,...
        'Units','pixels',...
        'Position',[10,10,25,230]);
    
fn_imagearray(vt_z',vt_x,vt_y,'Parent',ob_ax,'colormap',mx_map);

set(ob_ax,...
    'Box','on',...
    'Xtick',[],...
    'YAxisLocation','right')

ylabel(ch_label)