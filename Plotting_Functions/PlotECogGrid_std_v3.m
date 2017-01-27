function PlotECogGrid_std_v3(ERPs, timerange, varargin)
%% This function takes a input matricies of (arbitrary n_channels) x timepts containing
% various statisitics that applly over the ECoG channels over time,
% and plots the time evolution of these parameters in each channel.
%
% If there are not 256 channels, electrodes will be plotted in an arbitrary
% ordering.
% It is also possible to alter this order by specifying a grid order with
% the variable flag (...,'grid_order',grd,...)
%
% Inputs:
% ERPs - standard ERPs strucutured array that contains the information on
% the time axis, bad channels and bad trials
%
% timerange - 2x1 array that contains the start and stop times in SECONDS for ERP
% plots (will default to the maximum range of ERPs.time_axis)
% (ex: [-1, 2], sets the time window from 1 second before event to 2
% seconds after
%
% Variable inputs:
% Flag INPUTS: 
% (inputs that must be preceded by a flag)
%
% 'SigPts', sig_timepts_mat: FLAG INPUT - (boolean (nchans x ntimepts) matrix)) to add shading on the basis of a (nchans x ntimepts)
% boolean matrix indicating significant time points in each channel
% ex:
% PlotECogGrid_std(ERPs, [-1 1], 'SigPts', sig_timepts_mat, erps_1(:,:,is_1), erps_1(:,:,is_2), erps_2(:,:,is_1))
%
% 
% The tag normalize, will normalize all plots to a maximum absolute value
% of 1. If empty, this feature is turned off.
%
% 'Normalize', true: FLAG INPUT - (boolean) normalize all plots maximum
% absolute value of 1 (not including error).
%
% 'grid_order', matrix (with num_grid integers): FLAG INPUT specfies the order
% or plots in the tiling (...,'grid_order',grd,...)
%
% 'label_anatomy' - n_chans cell array: FLAG INPUT specifies the
% corresponding anatomical area for a given electrode
%
% 'Title', 'title_description': FLAG INPUT - (string) will add a title to
% be centered above the grid (Formatting this text is difficult, so keep it
% simple)
% Variable Inputs (Non-Flag):
% data matricies: 256xtimeptsxnum_trials data matricies of ECoG data
% (can be any number of them)
% if they contain many trials (size(dat,3) > 1), shaded error plots
% otherwise will make a simply line plot (useful for statistics)
%
%

%% Normalize plots 'Normalize'
% Sets maximum displacement to 1
Normalize = false;
norm_flag_ind = find(strcmpi(varargin, 'Normalize'));
if ~isempty(norm_flag_ind)
    if length(norm_flag_ind) ~= 1
        warning('Multiple Normalization Flags Given')
    end
    Normalize = varargin{norm_flag_ind(1) + 1};
    varargin(norm_flag_ind(1):(norm_flag_ind(1)+1)) = [];
end

%% significant time points for shading ('SigPts')
sig_flag_ind = find(strcmpi(varargin, 'SigPts'));
include_sig_timpts = false;
if ~isempty(sig_flag_ind)
    if length(sig_flag_ind) ~= 1
        warning('Multiple Significant TimePoint Matricies Given')
    end
    include_sig_timpts = true;
    sig_timepts_mat = varargin{sig_flag_ind(1) + 1};
    varargin(sig_flag_ind(1):(sig_flag_ind(1)+1)) = [];
end

%% Variable Grid Tiling Order
grid_flag_ind = find(strcmpi(varargin,'grid_order'));
if ~isempty(grid_flag_ind)
    if length(grid_flag_ind) ~= 1
        warning('Multiple Significant Grid Flags given')
    end
    grd = varargin{grid_flag_ind+1};
    ERPs.grd = grd;
    varargin(grid_flag_ind:(grid_flag_ind+1)) = [];
end


%% Label Anatomy (Color Code too)
anatomy_flag_ind = find(strcmpi(varargin,'label_anatomy'));
label_anatomy = false;
color_data = false;
if ~isempty(anatomy_flag_ind)
    color_data = true; % Potentially a variable with additional flags later...
    
    if length(anatomy_flag_ind) ~= 1
        warning('Multiple Anatomyy Flags given')
    end
    label_anatomy = true;
    anatomy_data = varargin{anatomy_flag_ind+1};
    varargin(anatomy_flag_ind:(anatomy_flag_ind+1)) = [];
    
    anatomical_areas = unique(anatomy_data);
    if color_data
        anatomy_colors = 0.8*distinguishable_colors(length(anatomical_areas));
    end
    if length(anatomy_data) ~= size(varargin{1},1)
        warning('Anatomy Labels Dimensions Inconsistent with ECoG Channels')
        label_anatomy = false;
    end
end


%% add Title
add_title_flag_ind = find(strcmpi(varargin, 'Title'));
add_title = false;
if ~isempty(add_title_flag_ind)
    if length(sig_flag_ind) ~= 1
        warning('Multiple Title Flags Given')
    end
    add_title = true;
    plot_title = varargin{add_title_flag_ind(1) + 1};
    varargin(add_title_flag_ind(1):(add_title_flag_ind(1)+1)) = [];
end

%% PLOTTTING:




num_plots = length(varargin);

%% Get Time Axis:
[~,min_time_ind] = min(abs(ERPs.time_axis - timerange(1)));
[~,max_time_ind] = min(abs(ERPs.time_axis - timerange(2)));
timerange = min_time_ind:max_time_ind;
time_axis = ERPs.time_axis(timerange);



trial_size = zeros(num_plots,1);
is_shaded_error_bar = true(num_plots,1); % plot shaded error lines or just line
for i = 1:num_plots
    trial_size(i) = size(varargin{i},3);
    if trial_size(i) == 1
        varargin{i} = squeeze(varargin{i});
        is_shaded_error_bar(i) = false;
    end
end
trial_size

%% Clear Bad Channels:
clear_bad_chan = true;
if clear_bad_chan
    Bad_Channels = ERPs.BadChans(ERPs.BadChans <= size(varargin{1},1)); % only include bad channels in range
    for i = 1:num_plots
        tmp = varargin{i};
        if is_shaded_error_bar(i)
            tmp(Bad_Channels,:,:) = 0;
        else
            tmp(Bad_Channels,:) = 0;
        end
        varargin{i} = tmp;
    end
    if include_sig_timpts
        sig_timepts_mat(Bad_Channels,:) = false;
    end
end
    
%% Get Global axis scaling for plots- (Normalize if Flagged)
maxbounds = zeros(num_plots,1);
minbounds = zeros(num_plots,1);

if Normalize
    for i = 1:num_plots
        dat = varargin{i};
        if is_shaded_error_bar(i)
            dat = dat/max(max(abs(mean(dat,3))));
        else
            dat = dat/max(max(abs(dat)));
        end
        varargin{i} = dat;
        
        if is_shaded_error_bar(i)
            maxbounds(i) = max(max(mean(dat,3)+nansem(dat,3)));
            minbounds(i) = min(min(mean(dat,3) - nansem(dat,3)));
        else
            maxbounds(i) = max(max(dat));
            minbounds(i) = min(min(dat));
        end
    end
else
    for i = 1:num_plots
        dat = varargin{i};
        if is_shaded_error_bar(i)
            maxbounds(i) = max(max(mean(dat,3)+nansem(dat,3)));
            minbounds(i) = min(min(mean(dat,3) - nansem(dat,3)));
        else
            maxbounds(i) = max(max(dat));
            minbounds(i) = min(min(dat));
        end
    end
end
maxbound = max(maxbounds);
minbound = min(minbounds);
bounds = [minbound, maxbound]
    
    
%% Plot Grid
n_chans = size(varargin{i},1);
% plotting parameters:
colorlist = [0.8 0 0; 0 0 0.8;0.8 0.8 0; 0 0.8 0.8; 0 0.8 0; 0.8 0 0.8];    % Plots of colors listed as Red, Blue, Yellow, Teal, Green, Purple
grid_plot = true;
grd = rot90(reshape(1:256,16,16)',3);
if isfield(ERPs, 'grd')
    grd = ERPs.grd;
else
    warning 'grid field not found in ERPs structure'
end
if n_chans ~= 256
    grid_plot = false;
end

figure('units','normalized','outerposition',[0 0 1 1]);
for k = 1:num_plots
    dat = varargin{k};
    if is_shaded_error_bar(k)
        dat_mean = mean(varargin{k},3);
        dat_sem = nansem(varargin{k},3);
    else
        dat = varargin{k};
    end
    for i = 1:n_chans
        if grid_plot
            p = plotGridPosition(grd(i)); 
        else
            p = plotGridPosition(grd(i),n_chans, ceil(sqrt(n_chans)));
        end
        subplot('position',p);
        if is_shaded_error_bar(k)
            shadedErrorBar(time_axis,dat_mean(i,timerange), dat_sem(i,timerange),{'color', colorlist(k,:)},1);
        else
            plot(time_axis, dat(i,timerange), 'Color', colorlist(k,:));
        end
        hold on;
        %% Plot annotations after all data is plotted 
        if k == num_plots
            axis tight;
            set(gca,'YLim',[minbound maxbound]...
                ,'XTickLabel',[],'YTickLabel',[]); %'YTick',[minbound maxbound]);
            line([0 0],get(gca,'YLim'),'Color','k');
            line(get(gca,'XLim'),[0 0],'Color','k');
            text(0,(max(get(gca,'YLim')) - (max(get(gca,'YLim')) * 0.1)),num2str(i));

            %% Add shading
            if include_sig_timpts
                sig_timepts_ch = sig_timepts_mat(i,timerange);
                shaded_patch_significant_timepoints(time_axis, sig_timepts_ch)
                axis tight; set(gca,'YLim',[minbound maxbound]);
            end
            %% Add Anatomy Labels:
            if color_data
                anatomy_area_ind = find(strcmpi(anatomical_areas, anatomy_data(i)),1);
                anatomy_colors(anatomy_area_ind,:);
                set(gca, 'XColor', anatomy_colors(anatomy_area_ind,:));
                set(gca, 'YColor', anatomy_colors(anatomy_area_ind,:));
                set(gca, 'box', 'off')
                set(gca, 'LineWidth', 4);
                set(gca, 'XTick', []); set(gca, 'YTick',[]);
            end
            
        end
    end
    %% Title Anatomy:
    if label_anatomy
        for i = 1:n_chans
            if grid_plot
                p = plotGridPosition(grd(i)); 
            else
                p = plotGridPosition(grd(i),n_chans,ceil(n_chans/floor(sqrt(n_chans))))
            end
            subplot('position',p);
            anatomy_area_ind = find(strcmpi(anatomical_areas, anatomy_data(i)),1);
            title(anatomical_areas{anatomy_area_ind}, 'FontSize',8);
        end
    end

end
%% Add title:
if add_title
    annotation('textbox', 'String', plot_title, ...
        'Position', [0, 0.98, 1, 0.02], ...
        'HorizontalAlignment','center','LineStyle', 'none');
end
a = 1;

end

