function st_handles	= fn_dataplot(mx_data,st_pref)
% fn_dataplot.m
% st_handles = fn_dataplot(mx_data) plots a visualization plot for the
% values in mx_data.
% 
% Inputs:
% 
%	mx_data: If matrix, then each level is determined colomnwise, 
%           and categories/factors are defined along the 3rd 
%           dimension.
%            If cell, values in each cell are considered within the
%           same distribution. Levels are considered colomnwise, 
%           and categories/factors are defined along the 2nd 
%           dimension.
%
%	st_pref: Structure indicating plot settings:
% 
%       FIGURE TYPES >>
%       st_pref.dotplot:    Display dot plot [true/false] (default: true)
%       st_pref.boxplot:	Display box plot [true/false] (default: false)
%       st_pref.violinplot:	Display violin plot [true/false] (default: false)
%       st_pref.bargraph:	Display bar graph [true/false] (default: false)
%       
%       AXIS SETTINGS >>
%       st_pref.axes:       Parent axes
%       st_pref.xPos:       Tick positions for each distribution. 
%       st_pref.xLabels:	Tick labels for each distribution. Cell array
%                           (factor x level)
%                          
%       DOT PLOT SETTINGS >>
%       st_pref.dot.width:          Dot width (default: 5)
%       st_pref.dot.marker:         Matlab key for dots marks (default: .),
%                                   if empty, then no marker
%       st_pref.dot.color:          Matlab color vector or string
%                                   indicating the dots color. For
%                                   different colors, then input is a cell
%                                   matrix as the values in mx_data.
%       st_pref.dot.fillDots:       Fill dots (default: false)
%       st_pref.dot.jitter:         Apply dot jitter (default: true)
%       st_pref.dot.jitterWidth:	Relative jitter width (default: 0.9)
%       st_pref.dot.connect:        Connect dots between levels of the same
%                                   factor/category 
%       st_pref.dot.connectColor:	Matlab color vector or string
%                                   indicating the connection line color
%       st_pref.dot.connectLine:    Matlab LineStyle for connections (default: -),
%       st_pref.dot.meanColor:      Matlab color vector or string
%                                   indicating the mean color
%       st_pref.dot.meanMarker:     Matlab key for mean marks (default: -),
%                                   if empty, then no marker
%       st_pref.dot.meanSize:       Size for mean marks (default: 12),
%       st_pref.dot.errorBar:       Display eror bar (default: false)
%       st_pref.dot.errorColor:     Matlab color vector or string
%                                   indicating the error bar color
%       st_pref.dot.errorConnect:	Connect error bars beween levels
%       st_pref.dot.errorType:      Error to plot 'std','sem','95CI'
%                                   (default: '95CI')
% 
%       BOX PLOT SETTINGS >>
%       st_pref.box.width:          Relative box width (default: 0.9)
%       st_pref.box.color:          Matlab color vector or string
%                                   indicating the box color. For
%                                   different colors, then input is a cell
%                                   matrix as the values in mx_data.
%       st_pref.box.fill:           Fill box with color (default: false)
%       st_pref.box.wiskerColor:    Matlab color vector or string
%                                   indicating the wisker color
%       st_pref.box.medianColor:    Matlab color vector or string
%                                   indicating the median color
%       st_pref.box.meanColor:      Matlab color vector or string
%                                   indicating the mean color
%       st_pref.box.outlierColor:	Matlab color vector or string
%                                   indicating the outliers color
%       st_pref.box.outlierMarker:	Matlab key for outlier marks (default: o)
%       st_pref.box.meanMarker:     Matlab key for mean marks (default: x),
%                                   if empty, then no marker
%       st_pref.box.wiskerLine:     Matlab LineStyle for wiskers (default: --),
%       st_pref.box.isNotched:      Notch the box if st_pref.ConfidenceI is
%                                   set (default: false)
%       st_pref.box.notchWidth:     Relative notch width (default: 0.6)
%       st_pref.box.isFenced:       Plot wisker fences (default: true)
%       st_pref.box.fenceWidth:     Relative box width (default: 0.2)
% 
%       VIOLIN PLOT SETTINGS >>
%       st_pref.violin.width:       Relative box width (default: 0.9)
%       st_pref.violin.color:       Matlab color vector or string
%                                   indicating the box color. For
%                                   different colors, then input is a cell
%                                   matrix as the values in mx_data.
%       st_pref.violin.fill:           Fill violin with color (default: false)
%       st_pref.violin.split:       Violin sides 'left','rigth','full'
%                                   (default: 'full')
%       st_pref.violin.statsBar:    Display stats bar with violin 
%                                   (default: true)
%       st_pref.violin.statsColor:	Matlab color vector or string
%                                   indicating the stats bar color
%       st_pref.violin.medianColor:	Matlab color vector or string
%                                   indicating the median color
% 
%       BAR GRAPH SETTINGS >>
%       st_pref.bar.stat:           Stat to display (default: 'mean')
%       st_pref.bar.width:          Relative box width (default: 0.9)
%       st_pref.bar.color:          Matlab color vector or string
%                                   indicating the box color. For
%                                   different colors, then input is a cell
%                                   matrix as the values in mx_data.
%       st_pref.bar.meanMarker:     Matlab key for mean marks (default: []),
%                                   if empty, then no marker
%       st_pref.bar.meanColor:      Matlab color vector or string
%                                   indicating the mean color
%       st_pref.bar.meanMarker:     Matlab key for mean marks (default: +),
%                                   if empty, then no marker
%       st_pref.bar.meanSize:       Size for mean marks (default: 12),
%       st_pref.bar.fill:           Fill bars with color (default: false)
%       st_pref.bar.hide:           Hide bars  (default: true)
%       st_pref.bar.errorBar:       Display eror bar (default: true)
%       st_pref.bar.errorColor:     Matlab color vector or string
%                                   indicating the error bar color
%       st_pref.bar.errorConnect:	Connect error bars beween levels
%       st_pref.bar.errorType:      Error to plot 'std','sem','95CI'
%                                   (default: '95CI')
% 
% Miguel Navarrete
% CUBRIC
% 2017

%% GNU licence,
% Copyright (C) <2019>  <Miguel Navarrete>
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

%% Code starts here:


% Check input variables
if nargin < 2
    st_pref      = [];
end

if iscell(mx_data)
    mx_data	= cellfun(@(x) x(:),mx_data,'UniformOutput',false);
    
    [nm_factor,nm_level,nm_vals]	= size(mx_data); 
    
    if nm_vals > 1
        error(...
            '[fn_dataplot] - m x n x 1 just allowed for mx_data as a cell')
    end
else
    [~,nm_level,nm_factor]	= size(mx_data); 
    
    mx_datTmp	= mx_data;    
    mx_datTmp(mx_datTmp == Inf) = NaN;
    mx_data     = cell(nm_factor,nm_level);
    
    for ff = 1:nm_factor
        for ll = 1:nm_level
            mx_data{ff,ll}	= mx_datTmp(:,ll,ff);
        end
    end    
end

if nargin < 2 || isempty(st_pref)
    st_pref      = struct;
end

if ~isfield(st_pref,'dotplot')
    st_pref.dotplot	= true;
end
if ~isfield(st_pref,'boxplot')
    st_pref.boxplot	= false;
end
if ~isfield(st_pref,'violinplot')
    st_pref.violinplot	= false;
end
if ~isfield(st_pref,'bargraph')
    st_pref.bargraph	= false;
end

if ~isfield(st_pref,'axes')
    st_pref.axes	= gca;
end

if ~isfield(st_pref,'xPos')
    st_pref.xPos	= nan(nm_level+1,nm_factor);

    for kk = 1:numel(st_pref.xPos)
        st_pref.xPos(kk)	= kk;
    end
    
    st_pref.xPos	= st_pref.xPos(1:nm_level,:);
end

if ~isfield(st_pref,'xLabels')
    
    mx_rows	= repmat((1:size(mx_data,2))',1,size(mx_data,1));
    mx_cols	= repmat((1:size(mx_data,1)),size(mx_data,2),1);
    
    st_pref.xLabels = cell(size(mx_data));
    
    for kk = 1:numel(mx_data)
        st_pref.xLabels{kk}	= sprintf('(%i,%i)',mx_cols(kk),mx_rows(kk)) ;
    end
end
    
%% Figure settings
% cla(st_pref.axes)
mx_data         = mx_data';

st_pref.Section = min(diff(st_pref.xPos(:))) / 2; 
st_pref.xLabels	= st_pref.xLabels(:);
st_pref.xPos	= st_pref.xPos(:);
st_pref.xPos    = reshape(st_pref.xPos,nm_level,nm_factor);

if isempty(st_pref.Section)
    st_pref.Section = 1;
end

%% Defaults for dot plot
if st_pref.dotplot
    
    if ~isfield(st_pref,'dot')
        st_dot	= struct;
    else
        st_dot	= st_pref.dot;
    end
        
    if ~isfield(st_dot,'jitter')
        st_dot.jitter	= true;
    end
    if ~isfield(st_dot,'jitterWidth')
        st_dot.jitterWidth	= 0.9;
    end
    if ~isfield(st_dot,'marker')
        st_dot.marker     = '.';
    end
    if ~isfield(st_dot,'width')
        st_dot.width	= 5;
    end
    if ~isfield(st_dot,'color')
        st_dot.color	= [0    0.4510    0.7412];
    end
    if ~isfield(st_dot,'fillDots')
        st_dot.fillDots	= false;
    end
    if ~isfield(st_dot,'connect')
        st_dot.connect 	= false;
    end
    if ~isfield(st_dot,'connectLine')
        st_dot.connectLine   = '-';
    end
    if ~isfield(st_dot,'connectColor')
        st_dot.connectColor	= [0.7 0.7 0.7];
    end
    if ~isfield(st_dot,'meanColor')
        st_dot.meanColor	= [0.5 0.0 0.0];
    end
    if ~isfield(st_dot,'meanMarker')
        st_dot.meanMarker     = '-';
    end
    if ~isfield(st_dot,'meanSize')
        st_dot.meanSize	= 3;
    end
    
    if ~isfield(st_dot,'errorBar')
        st_dot.errorBar	= false;
    end
    if ~isfield(st_dot,'errorColor')
        st_dot.errorColor	= [0 0 0];
    end
    if ~isfield(st_dot,'errorConnect')
        st_dot.errorConnect	= false;
    end
    if ~isfield(st_dot,'errorType')
        st_dot.errorType	= '95CI';
    end
    
    if ~iscell(st_dot.color)        
        st_dot.color	= repmat({st_dot.color},...
                        nm_factor,nm_level);
    end
    
    st_dot.color	= st_dot.color';
end

%% Defaults for box plot
if st_pref.boxplot
    
    if ~isfield(st_pref,'box')
        st_box	= struct;
    else
        st_box	= st_pref.box;
    end
        
    if ~isfield(st_box,'width')
        st_box.width	= 0.9;
    end
    if ~isfield(st_box,'color')
        st_box.color	= [0.5 0.5 0.5];
    end
    if ~isfield(st_box,'fill')
        st_box.fill     = true;
    end
    if ~isfield(st_box,'wiskerColor')
        st_box.wiskerColor	= [0.5 0.5 0.5];
    end
    if ~isfield(st_box,'medianColor')
        st_box.medianColor	= [0.5 0.0 0.0];
    end
    if ~isfield(st_box,'meanColor')
        st_box.meanColor	= [0.0 0.0 0.5];
    end
    if ~isfield(st_box,'outlierColor')
        st_box.outlierColor	= [0.0 0.0 0.5];
    end
    if ~isfield(st_box,'outlierMarker')
        st_box.outlierMarker	= 'o';
    end
    if ~isfield(st_box,'meanMarker')
        st_box.meanMarker     = 'x';
    end
    if ~isfield(st_box,'wiskerLine')
        st_box.wiskerLine   = '--';
    end
    if ~isfield(st_box,'isNotched')
        st_box.isNotched 	= false;
    end
    if ~isfield(st_box,'fenceWidth')
        st_box.notchWidth	= 0.5;
    end
    if ~isfield(st_box,'isFenced')
        st_box.isFenced     = true;
    end
    if ~isfield(st_box,'fenceWidth')
        st_box.fenceWidth	= 0.2;
    end
    
    if ~iscell(st_box.color)        
        st_box.color	= repmat({st_box.color},...
                        nm_factor,nm_level);
    end
    
    st_box.color	= st_box.color';
end

%% Defaults for violin plot
if st_pref.violinplot
    
    if ~isfield(st_pref,'violin')
        st_violin	= struct;
    else
        st_violin	= st_pref.violin;
    end
        
    if ~isfield(st_violin,'width')
        st_violin.width	= 0.9;
    end
    if ~isfield(st_violin,'color')
        st_violin.color	= [0.1 0.1 0.5];
    end
    if ~isfield(st_violin,'fill')
        st_violin.fill      = true;
    end
    if ~isfield(st_violin,'split')
        st_violin.split     = 'full';
    end
    if ~isfield(st_violin,'statsBar')
        st_violin.statsBar	= true;
    end
    if ~isfield(st_violin,'statsColor')
        st_violin.statsColor	= 'k';
    end
    if ~isfield(st_violin,'medianColor')
        st_violin.medianColor	= [0.5 0.0 0.0];
    end
    
    if ~iscell(st_violin.color)        
        st_violin.color	= repmat({st_violin.color},...
                        nm_factor,nm_level);
    end
    
    st_violin.color	= st_violin.color';
end

%% Defaults for bar graph
if st_pref.bargraph
    
    if ~isfield(st_pref,'bar')
        st_bar	= struct;
    else
        st_bar	= st_pref.bar;
    end
        
    if ~isfield(st_bar,'width')
        st_bar.width	= 0.9;
    end
    if ~isfield(st_bar,'color')
        st_bar.color	= [0.5 0.5 0.5];
    end
    if ~isfield(st_bar,'meanColor')
        st_bar.meanColor	= [0.5 0.0 0.0];
    end
    if ~isfield(st_bar,'meanMarker')
        st_bar.meanMarker     = [];
    end
    if ~isfield(st_bar,'meanSize')
        st_bar.meanSize	= 6;
    end
    if ~isfield(st_bar,'fill')
        st_bar.fill     = true;
    end
    if ~isfield(st_bar,'hide')
        st_bar.hide     = true;
    end
    if ~isfield(st_bar,'errorBar')
        st_bar.errorBar	= true;
    end
    if ~isfield(st_bar,'errorColor')
        st_bar.errorColor	= [0 0 0];
    end
    if ~isfield(st_bar,'errorConnect')
        st_bar.errorConnect	= true;
    end
    if ~isfield(st_bar,'errorType')
        st_bar.errorType	= '95CI';
    end
    if ~isfield(st_bar,'stat')
        st_bar.stat	= 'mean';
    end
    
    if ~iscell(st_bar.color)        
        st_bar.color	= repmat({st_bar.color},...
                        nm_factor,nm_level);
    end
    
    st_bar.color	= st_bar.color';
end

%% Layer 4 - Bar graph

if st_pref.bargraph    
    mx_measures	= nan(size(mx_data));
    mx_error 	= nan(size(mx_data));

    for kk = 1:numel(mx_data)
        vt_curData      = mx_data{kk}(:);
        vt_curData      = vt_curData(~isnan(vt_curData));
        
        switch lower(st_bar.stat)
            case 'mean' 
                mx_measures(kk)	= mean(vt_curData);
            case 'median'
                mx_measures(kk)	= median(vt_curData);
        end

        switch lower(st_bar.errorType)
            case 'std' 
                mx_error(kk)	= std(vt_curData);
            case 'sem'
                mx_error(kk)	= fn_sem(vt_curData);
            case '95ci'
                mx_error(kk)	= 1.96*fn_sem(vt_curData);
        end
        
    end
    
	fn_figure_bargraph(mx_measures,mx_error);
end

%% Layer 3 - Violin plot

if st_pref.violinplot    
    mx_measures	= nan(numel(mx_data),5);
    vt_violin	= cell(numel(mx_data),1);

    for kk = 1:numel(mx_data)
        vt_curData              = mx_data{kk}(:);
        vt_curData              = vt_curData(~isnan(vt_curData));
        [~,vt_freqs,vt_values]	= kde(vt_curData,2^6);
        
        vt_freqs	= vt_freqs/max(vt_freqs);
        
        switch st_violin.split
            case 'left'
                vt_freqs	= -vt_freqs(:);
                vt_values	= vt_values(:);
            case 'right'
                vt_freqs	= vt_freqs(:);
                vt_values	= vt_values(:);
            otherwise
                vt_freqs	= [vt_freqs(:);flipud(-vt_freqs(:))];
                vt_values	= [vt_values(:);flipud(vt_values(:))];
        end
        
        mx_measures(kk,:)	= fn_compute_statistics(vt_curData);
        vt_violin{kk}       = [vt_freqs(:),vt_values(:)];
    end
    
	fn_figure_violinplot(mx_measures,vt_violin);
end

%% Layer 2 - Box plot

if st_pref.boxplot    
    mx_measures	= nan(numel(mx_data),5);
    mx_outliers	= cell(numel(mx_data),2);
    vt_means    = nan(numel(mx_data),1);
    vt_95CI     = nan(numel(mx_data),1);

    for kk = 1:numel(mx_data)
        vt_curData      = mx_data{kk}(:);
        vt_curData      = vt_curData(~isnan(vt_curData));

        [vt_curMeas,vt_curOutl]	= fn_compute_statistics(vt_curData);

        mx_measures(kk,:)       = vt_curMeas;
        mx_outliers(kk,:)       = vt_curOutl;
        vt_means(kk,1)          = mean(vt_curData);
        vt_95CI(kk,1)           = 3.14 * diff(vt_curMeas([2,4])) / ...
                                sqrt(numel(vt_curData)) ;
    end

    st_box.means        = vt_means;
    st_box.ConfidenceI  = vt_95CI;
    
	fn_figure_boxplot(mx_measures,mx_outliers);
end

%% Layer 1 - Dot plot

if st_pref.dotplot    
    vt_dataFactor	= cell(nm_factor,1);
    vt_dataPos      = cell(nm_factor,1);
        
    nm_dotsBin      = st_pref.Section * st_dot.jitterWidth;
    
    for ff = 1:numel(vt_dataFactor)
        
        vt_curFactor	= mx_data(:,ff);
        vt_numFactor    = cellfun(@numel,vt_curFactor);
        nm_maxNumel     = max(vt_numFactor);
        
        mx_curLevels    = nan(nm_maxNumel,nm_level);
        
        for lv = 1:nm_level
            mx_curLevels(1:vt_numFactor(lv),lv)	= vt_curFactor{lv}(:);
        end
        
        mx_xLevels	= nan(size(mx_curLevels));
        
        for lv = 1:nm_level
            nm_xTick    =  st_pref.xPos(lv,ff);
            if st_dot.jitter
                
                nm_iqr      = prctile(vt_curFactor{lv},75) - ...
                            prctile(vt_curFactor{lv},25);
                        
                nm_bin      = 2*nm_iqr/(numel(vt_curFactor{lv}))^(1/3);
                
                if nm_bin == 0                    
                    vt_binedges = [min(vt_curFactor{lv}),...
                                max(vt_curFactor{lv})+0.5];                            
                else
                    vt_binedges = min(vt_curFactor{lv}):nm_bin:...
                                max(vt_curFactor{lv});
                                            
                end
                
                if numel(vt_binedges) < 2
                    vt_binedges = [min(vt_curFactor{lv}),...
                                max(vt_curFactor{lv})];
                end
                
                if vt_binedges(end) ~= max(vt_curFactor{lv})
                    vt_binedges(end+1)	= vt_binedges(end) + nm_bin;    %#ok<AGROW>
                end
                
                vt_bin      = discretize(mx_curLevels(:,lv),vt_binedges);
                vt_binLevel	= unique(vt_bin);
                vt_binLevel = vt_binLevel(~isnan(vt_binLevel));
                vt_numCount = histcounts(mx_curLevels(:,lv),vt_binedges);
                nm_maxCount = max(vt_numCount);
                
                nm_trVal    = roots([1,1,-nm_maxCount]);
                nm_trVal    = ceil(nm_trVal(nm_trVal > 0));
                mx_pos      = rot90(triu(ones(nm_trVal)),2);
                mx_binPos   = linspace(0,1,nm_trVal+1);
                mx_binPos   = mx_binPos(2:end);
                nm_binShift = mx_binPos(1)/2;
                mx_binPos   = mx_binPos - nm_binShift;
                mx_binPos   = repmat(mx_binPos,nm_trVal,1);
                mx_binPos   = mx_binPos .* mx_pos;
                
                nm_shiftStep	= nm_binShift/1.5;
                vt_shiftValues  = [nm_shiftStep,0,-nm_shiftStep,0];
                nm_shCount      = 1;
                nm_shift        = vt_shiftValues(nm_shCount);
                
                for bd = 1:numel(vt_binLevel)
                    
                    vt_id   = vt_bin == vt_binLevel(bd);
                    [~,vt_idSort]	= sort(mx_curLevels(vt_id,lv));
                    [~,vt_idSort]   = sort(vt_idSort);
                    
                    vt_pos  = nan(size(vt_idSort));
                    nm_row  = 1;
                    nm_col  = 1;
                    for pd = 1:numel(vt_pos)
                        if mod(pd,2)
                            vt_pos(pd)  = mx_binPos(nm_row,nm_col)...
                                        - nm_shift;
                        else
                            vt_pos(pd) = -mx_binPos(nm_row,nm_col)...
                                        - nm_shift;
                            nm_col      = nm_col + 1;
                        end
                                                
                        if nm_row > size(mx_binPos,1)
                            nm_row	= 1;
                        end
                        if nm_col > size(mx_binPos,2)
                            nm_col	= 1;
                        end
                        
                        if mx_binPos(nm_row,nm_col) == 0
                            nm_row  = nm_row + 1;
                            nm_col  = 1;
                            
                            nm_shCount	= nm_shCount + 1;
                            if nm_shCount > numel(vt_shiftValues)
                                nm_shCount	= 1;
                            end
                                
                            nm_shift= vt_shiftValues(nm_shCount);
                        end
                        
                    end
                                      
                    vt_pos	= nm_dotsBin .* vt_pos(vt_idSort) + nm_xTick;
                    
                    mx_xLevels(vt_id,lv)	= vt_pos;
                end
                
            else
                mx_xLevels(:,lv)	= nm_xTick;
                
            end
        end
        vt_dataFactor{ff}	= mx_curLevels;
        vt_dataPos{ff}      = mx_xLevels;
        
    end

	fn_figure_dotplot(vt_dataPos,vt_dataFactor);
end

%% Set axes

set(st_pref.axes,...
    'xTick',st_pref.xPos(:),...
	'xtickLabel',st_pref.xLabels)

%% Nested functions :::::::::::::::::::::::::::::::::::::::::::::::::::::::
    function fn_figure_bargraph(mx_measures,mx_error)
    % Display fn_figure_bar graph 
    
    nm_boxLim	= st_pref.Section * st_bar.width ;
    vt_barPos   = nan(size(mx_measures));

    for bb = 1:numel(mx_measures)

        % plot the current box
        nm_curBin   = st_pref.xPos(bb);
        nm_curValue = mx_measures(bb);
        
        vt_barPos(bb)	= nm_curBin;
        
        if any(isnan(mx_measures(bb)))
            continue
        end

        if st_bar.hide
            continue
        end

        if st_bar.fill
            
            st_hBox.Patch	= patch(st_pref.axes,...
                            'XData',[nm_curBin - nm_boxLim,...
                                    nm_curBin - nm_boxLim,...
                                    nm_curBin + nm_boxLim,...
                                    nm_curBin + nm_boxLim,...
                                    nm_curBin - nm_boxLim],...
                            'YData',[0,...
                                    nm_curValue,...
                                    nm_curValue,...
                                    0,0],...
                            'EdgeColor','none',...
                            'FaceAlpha',0.8,...
                            'FaceColor',st_bar.color{bb});
        end
        
        st_hBox.Bar	= line(st_pref.axes,...
                            'XData',[nm_curBin - nm_boxLim,...
                                    nm_curBin - nm_boxLim,...
                                    nm_curBin + nm_boxLim,...
                                    nm_curBin + nm_boxLim,...
                                    nm_curBin - nm_boxLim],...
                            'YData',[0,...
                                    nm_curValue,...
                                    nm_curValue,...
                                    0,0],...
                            'Color',st_bar.color{bb});

    end
    
    if st_bar.errorBar
        if st_bar.errorConnect
            ch_line	= '-';
        else
            ch_line	= 'none';            
        end
        
        for fb = 1:size(mx_error,2)
            
            vt_xVal	= vt_barPos(:,fb);
            vt_yVal	= mx_measures(:,fb);
            vt_eVal	= mx_error(:,fb);
            hold on
            st_hBox.error{fb}	= errorbar(st_pref.axes,...                
                                vt_xVal,...
                                vt_yVal,...
                                vt_eVal,...
                                'Color',st_bar.errorColor,...
                                'LineStyle',ch_line,...
                                'LineWidth',1.5);
            
            if ~isempty(st_bar.meanMarker)
                st_hBox.mean{fb}	= plot(st_pref.axes,vt_xVal,vt_yVal,...
                                    'Marker',st_bar.meanMarker,...
                                    'MarkerSize',st_bar.meanSize,...
                                    'MarkerEdgeColor',st_bar.meanColor,...
                                    'MarkerFaceColor',st_bar.meanColor,...
                                    'LineStyle','none');

            end
            hold off
        end
    end
    
    end
%--------------------------------------------------------------------------
    function fn_figure_dotplot(vt_dataPos,vt_dataFactor)
    % Display dot plot diagram
    
    nm_dotsLim      = st_pref.Section * st_dot.jitterWidth;
    
    if ~st_dot.connect
        st_dot.connectLine = 'none';
    end
    
    for fd = 1:numel(vt_dataFactor)
                
        mx_dval = vt_dataFactor{fd};
        mx_dpos = vt_dataPos{fd};
        
        for pp = 1:size(mx_dval,2)
            
            if st_dot.fillDots
                vt_fill = st_dot.color{pp,fd};
            else
                vt_fill = 'none';
            end

            hold on
            st_hDot.Dots{fd}	= plot(st_pref.axes,...
                                mx_dpos(:,pp),mx_dval(:,pp),...%
                                'Marker',st_dot.marker,...
                                'MarkerSize',st_dot.width,...
                                'MarkerEdgeColor',st_dot.color{pp,fd},...
                                'MarkerFaceColor',vt_fill,...
                                'LineStyle','none',...
                                'Color',st_dot.connectColor);

            if ~isempty(st_dot.meanMarker) && strcmpi(st_dot.meanMarker,'-')

                nm_mean	= mean(mx_dval(:,pp),'omitnan');
                nm_xPos = st_pref.xPos(pp,fd);

                st_hDot.mean{fd}{pp}= line(st_pref.axes,...
                                    'XData',[nm_xPos - nm_dotsLim,...
                                            nm_xPos + nm_dotsLim],...
                                    'YData',[nm_mean,nm_mean],...
                                    'LineWidth',st_dot.meanSize,...
                                    'Color',st_dot.meanColor);

            elseif  ~isempty(st_dot.meanMarker) && ~strcmpi(st_dot.connectLine,'-')
                vt_mean	= mean(vt_dataFactor{fd},'omitnan');
                vt_xPos = st_pref.xPos(:,fd);

                st_hDot.mean{fd}	= plot(st_pref.axes,vt_xPos,vt_mean,...
                                    'Marker',st_dot.meanMarker,...
                                    'MarkerSize',st_dot.meanSize,...
                                    'MarkerEdgeColor',st_dot.meanColor,...
                                    'LineStyle','none');            

            end         
        end
        
        if st_dot.connect
            
            st_hDot.Connect{fd}	= plot(st_pref.axes,...
                                vt_dataPos{fd}',mx_dval',...
                                'Marker','none',...
                                'Color',st_dot.connectColor,...
                                'LineStyle',st_dot.connectLine);
        end
        
        
        if st_dot.errorBar  
            if st_dot.errorConnect
                ch_line	= '-';
            else
                ch_line	= 'none';            
            end          
            vt_dotmean	= mean(mx_dval,'omitnan');
            switch lower(st_dot.errorType)
                case 'std' 
                    vt_doterror	= std(mx_dval,'omitnan');
                case 'sem'
                    vt_doterror	= fn_sem(mx_dval);
                case '95ci'
                    vt_doterror	= 1.96*fn_sem(mx_dval);
            end
            
            vt_xVal	= st_pref.xPos(:,fd);
            vt_yVal	= vt_dotmean;
            vt_eVal	= vt_doterror;
            
            st_hDot.error{fd}	= errorbar(st_pref.axes,...                
                                vt_xVal,...
                                vt_yVal,...
                                vt_eVal,...
                                'Color',st_dot.errorColor,...
                                'LineStyle',ch_line,...
                                'LineWidth',3);
        end
        hold off
        
    end
       
    end
%--------------------------------------------------------------------------
    function fn_figure_boxplot(mx_measures,mx_outliers)
    % Display box plot diagram
    
    nm_isMean	= ~isempty(st_box.meanMarker);

    mx_colors   = st_box.color;
    nm_boxLim	= st_pref.Section * st_box.width ;
    nm_fenceLim = st_pref.Section * st_box.fenceWidth;
    nm_notchLim	= st_pref.Section * st_box.notchWidth * st_box.width;

    if isempty(mx_outliers)
        mx_outliers     = cell(size(mx_measures,1),2);
    end

    for bb = 1:size(mx_measures,1)

        if any(isnan(mx_measures(bb,:)))
            continue
        end

        % plot the current box
        nm_curBin       = st_pref.xPos(bb);
        nm_curLoFence   = mx_measures(bb,1);
        nm_curQ1        = mx_measures(bb,2);
        nm_curQ2        = mx_measures(bb,3);
        nm_curQ3        = mx_measures(bb,4);
        nm_curUpFence	= mx_measures(bb,5);

        nm_curLoNotch   = nm_curQ2 - st_box.ConfidenceI(bb)/2;
        nm_curUpNotch   = nm_curQ2 + st_box.ConfidenceI(bb)/2;

        if st_box.fill
            
            st_hBox.Patch	= patch(st_pref.axes,...
                            'XData',[nm_curBin - nm_boxLim,...
                                    nm_curBin + nm_boxLim,...
                                    nm_curBin + nm_boxLim,...
                                    nm_curBin + nm_notchLim,...
                                    nm_curBin + nm_boxLim,...
                                    nm_curBin + nm_boxLim,...
                                    nm_curBin - nm_boxLim,...
                                    nm_curBin - nm_boxLim,...
                                    nm_curBin - nm_notchLim,...
                                    nm_curBin - nm_boxLim,...
                                    nm_curBin - nm_boxLim],...
                            'YData',[nm_curQ1,...
                                    nm_curQ1,...
                                    nm_curLoNotch,...
                                    nm_curQ2,...
                                    nm_curUpNotch,...
                                    nm_curQ3,...
                                    nm_curQ3,...
                                    nm_curUpNotch,...
                                    nm_curQ2,...
                                    nm_curLoNotch,...
                                    nm_curQ1],...
                            'FaceAlpha',.5,...
                            'FaceColor',mx_colors{bb});
        end
        
        st_hBox.Box     = line(st_pref.axes,...
                        'XData',[nm_curBin - nm_boxLim,...
                                nm_curBin + nm_boxLim,...
                                nm_curBin + nm_boxLim,...
                                nm_curBin + nm_notchLim,...
                                nm_curBin + nm_boxLim,...
                                nm_curBin + nm_boxLim,...
                                nm_curBin - nm_boxLim,...
                                nm_curBin - nm_boxLim,...
                                nm_curBin - nm_notchLim,...
                                nm_curBin - nm_boxLim,...
                                nm_curBin - nm_boxLim],...                            
                        'YData',[nm_curQ1,...
                                nm_curQ1,...
                                nm_curLoNotch,...
                                nm_curQ2,...
                                nm_curUpNotch,...
                                nm_curQ3,...
                                nm_curQ3,...
                                nm_curUpNotch,...
                                nm_curQ2,...
                                nm_curLoNotch,...
                                nm_curQ1],...
                        'Color',mx_colors{bb});

        st_hBox.Median	= line(st_pref.axes,...
                        'XData',[nm_curBin - nm_notchLim,...
                                nm_curBin + nm_notchLim],...                            
                        'YData',[nm_curQ2,...
                                nm_curQ2],...
                        'LineWidth',2,...
                        'Color',st_box.medianColor);

        st_hBox.LoWisker	= line(st_pref.axes,...
                            'XData',[nm_curBin,nm_curBin],...                            
                            'YData',[nm_curLoFence,nm_curQ1],...
                            'Color',st_box.wiskerColor,...
                            'LineStyle',st_box.wiskerLine);

        st_hBox.UpWisker	= line(st_pref.axes,...
                            'XData',[nm_curBin,nm_curBin],...                            
                            'YData',[nm_curQ3,nm_curUpFence],...
                            'Color',st_box.wiskerColor,...
                            'LineStyle',st_box.wiskerLine);

        if st_box.isFenced

            st_hBox.LoFence	= line(st_pref.axes,...
                            'XData',[nm_curBin - nm_fenceLim,...
                                    nm_curBin + nm_fenceLim],...                            
                            'YData',[nm_curLoFence,nm_curLoFence],...
                            'Color',st_box.wiskerColor,...
                            'LineStyle','-');

            st_hBox.UpFence	= line(st_pref.axes,...
                            'XData',[nm_curBin - nm_fenceLim,...
                                    nm_curBin + nm_fenceLim],...                            
                            'YData',[nm_curUpFence,nm_curUpFence],...
                            'Color',st_box.wiskerColor,...
                            'LineStyle','-');
        end     

        if nm_isMean
            st_hBox.Mean	= line(st_pref.axes,...
                            'XData',nm_curBin,...                            
                            'YData',st_box.means(bb),...
                            'Color',st_box.meanColor,...
                            'LineStyle','none',...
                            'Marker',st_box.meanMarker);
        end

        if ~isempty(st_box.outlierMarker) && ~isempty(mx_outliers{bb,1})
            vt_curOutliers	= mx_outliers{bb,1};
            vt_curBinValues	= (2 * nm_fenceLim).*rand(size(vt_curOutliers)) + ...
                            (nm_curBin - nm_fenceLim);

            st_hBox.LoOut	= line(st_pref.axes,...
                            'XData',vt_curBinValues,...                            
                            'YData',vt_curOutliers,...
                            'Color',st_box.outlierColor,...
                            'LineStyle','none',...
                            'Marker',st_box.outlierMarker);

        end

        if ~isempty(st_box.outlierMarker) && ~isempty(mx_outliers{bb,2})
            vt_curOutliers	= mx_outliers{bb,2};
            vt_curBinValues	= (2 * nm_fenceLim).*rand(size(vt_curOutliers)) + ...
                            (nm_curBin - nm_fenceLim);

            st_hBox.LoOut	= line(st_pref.axes,...
                            'XData',vt_curBinValues,...                            
                            'YData',vt_curOutliers,...
                            'Color',st_box.outlierColor,...
                            'LineStyle','none',...
                            'Marker',st_box.outlierMarker);

        end
    end
    end
%--------------------------------------------------------------------------
    function fn_figure_violinplot(mx_measures,vt_violin)
    % Display box plot diagram
    
    if isempty(st_pref.Section)
        st_pref.Section = 0.5;
    end

    nm_boxLim	= st_pref.Section * st_violin.width ;

    for bb = 1:size(mx_measures,1)

        if any(isnan(mx_measures(bb,:)))
            continue
        end

        % plot the current box
        nm_curBin       = st_pref.xPos(bb);
        nm_curLoFence   = mx_measures(bb,1);
        nm_curQ1        = mx_measures(bb,2);
        nm_curQ2        = mx_measures(bb,3);
        nm_curQ3        = mx_measures(bb,4);
        nm_curUpFence	= mx_measures(bb,5);

        vt_xData = vt_violin{bb}(:,1) .* nm_boxLim + nm_curBin;
        vt_yData = vt_violin{bb}(:,2);

        if st_violin.fill
            
            st_hViolin.Patch     = patch(st_pref.axes,...
                                'XData',vt_xData,...                            
                                'YData',vt_yData,...
                                'FaceAlpha',.5,...
                                'FaceColor',st_violin.color{bb});
        end
        
        st_hViolin.Violin     = line(st_pref.axes,...
                            'XData',vt_xData,...                            
                            'YData',vt_yData,...
                            'Color',st_violin.color{bb});

        if st_violin.statsBar
            
            st_hViolin.Wiskers	= line(st_pref.axes,...
                                'XData',[nm_curBin,nm_curBin],...                            
                                'YData',[nm_curLoFence,nm_curUpFence],...
                                'Color',st_violin.statsColor,...
                                'LineStyle','-');

            st_hViolin.lineIQR	= line(st_pref.axes,...
                                'XData',[nm_curBin,nm_curBin],...                            
                                'YData',[nm_curQ1,nm_curQ3],...
                                'Color',st_violin.statsColor,...
                                'LineWidth',6,...
                                'LineStyle','-');

            st_hViolin.Median	= line(st_pref.axes,...
                            'XData',nm_curBin,...                            
                            'YData',nm_curQ2,...
                            'Marker','o',...
                            'MarkerFaceColor','w',...
                            'MarkerSize',6,...
                            'Color',st_violin.medianColor);
        end

    end
    end
%--------------------------------------------------------------------------
    function [vt_stats,st_outl]	= fn_compute_statistics(vt_dat)
        
        
        vt_dat	= vt_dat(:);
        vt_dat	= vt_dat(~isnan(vt_dat));

        % Computing measures for first click
        nm_lowQtle  = prctile(vt_dat,25);
        nm_midQtle  = prctile(vt_dat,50);
        nm_uppQtle  = prctile(vt_dat,75);

        nm_minVal   = min(vt_dat);
        nm_maxVal   = max(vt_dat);

        nm_minExt   = nm_lowQtle - 1.5 * (nm_uppQtle - nm_lowQtle);
        nm_maxExt   = nm_uppQtle + 1.5 * (nm_uppQtle - nm_lowQtle);

        vt_lowOutl  = vt_dat(vt_dat < nm_minExt);
        vt_uppOutl  = vt_dat(vt_dat > nm_maxExt);

        if ~isempty(vt_lowOutl)
            nm_minVal = nm_minExt;
        end

        if ~isempty(vt_uppOutl)
            nm_maxVal = nm_maxExt;
        end

        if isempty(nm_minVal)
            nm_minVal = nan;
        end

        if isempty(nm_maxVal)
            nm_maxVal = nan;
        end

        vt_stats(1,1)	= nm_minVal;
        vt_stats(1,2)	= nm_lowQtle;
        vt_stats(1,3)	= nm_midQtle;
        vt_stats(1,4)	= nm_uppQtle;
        vt_stats(1,5)	= nm_maxVal;

        st_outl{1,1}	= vt_lowOutl;
        st_outl{1,2}    = vt_uppOutl;

    end
end