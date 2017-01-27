function fig = PlotECogGrid_std_v2(ERPs, timerange, varargin)
%% This function takes a input matricies of 256 x timepts containing
% various statisitics that applly over the ECoG grid channels over time,
% and plots the time evolution of these parameters in each channel.
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
% 'Title', 'title_description': FLAG INPUT - (string) will add a title to
% be centered above the grid (Formatting this text is difficult, so keep it
% simple)
%
%
%
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
inclue_sig_timpts = false;
if ~isempty(sig_flag_ind)
    if length(sig_flag_ind) ~= 1
        warning('Multiple Significant TimePoint Matricies Given')
    end
    inclue_sig_timpts = true;
    sig_timepts_mat = varargin{sig_flag_ind(1) + 1};
    varargin(sig_flag_ind(1):(sig_flag_ind(1)+1)) = [];
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
    Bad_Channels = ERPs.BadChans(ERPs.BadChans <= 256); % only include bad channels in range
    for i = 1:num_plots
        tmp = varargin{i};
        if is_shaded_error_bar(i)
            tmp(Bad_Channels,:,:) = 0;
        else
            tmp(Bad_Channels,:) = 0;
        end
        varargin{i} = tmp;
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
% plotting parameters:
colorlist = [0.8 0 0; 0 0 0.8;0.8 0.8 0; 0 0.8 0.8; 0 0.8 0; 0.8 0 0.8];    % Plots of colors listed as Red, Blue, Yellow, Teal, Green, Purple
grd = rot90(reshape(1:256,16,16)',3);
if isfield(ERPs, 'grd')
    grd = ERPs.grd;
else
    warning 'grid field not found in ERPs structure'
end

fig = figure('units','normalized','outerposition',[0 0 1 1]);
for k = 1:num_plots
    dat = varargin{k};
    if is_shaded_error_bar(k)
        dat_mean = mean(varargin{k},3);
        dat_sem = nansem(varargin{k},3);
    else
        dat = varargin{k};
    end
    for i = 1:256
        p = plotGridPosition(grd(i)); 
        subplot('position',p);
        if is_shaded_error_bar(k)
            shadedErrorBar(time_axis,dat_mean(i,timerange), dat_sem(i,timerange),{'color',colorlist(k,:)},1);
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
            if inclue_sig_timpts
                sig_timepts_ch = sig_timepts_mat(i,timerange);
                shaded_patch_significant_timepoints(time_axis, sig_timepts_ch)
                axis tight; set(gca,'YLim',[minbound maxbound]);
            end
        end
    end

end

%% Add Title
if add_title
    annotation('textbox', 'String', plot_title, ...
        'Position', [0, 0.98, 1, 0.02], ...
        'HorizontalAlignment','center','LineStyle', 'none');
end
set(gcf, 'PaperPositionMode', 'auto')
a = 1;

end

