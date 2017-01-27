function plot_spectral_tstat(ERPs, spectral_dat_dir, is_class_1, is_class_2)
% Calculate the T-Statistic between class1 and class2 for each channel
% and plot this information.


% spectral_dat_dir =
% '/Users/changlab/Documents/data/EC125/ecog_data/spectrogram_data';
% is_class_1 = ERPs.is_good_trial & ERPs.is_high_fa;
% is_class_2 = ERPs.is_good_trial & ERPs.is_low_fa;
n_chans = size(ERPs.ecog_targets,1);

is_bad_chan = false(n_chans,1);
is_bad_chan(ERPs.BadChans) = true;
grd = ERPs.grd;

if n_chans ~= 256 % implies a non-grid configuration (
    grid_offset = 256; % assumes chan 1 is chan 257
else
    grid_offset = 0; % assumes chan 1 is chan 1
end

t_thresh = 1.5;
cmaxes = zeros(n_chans,1);
cmins = zeros(n_chans,1);
%load_basic_parameters:
load([spectral_dat_dir filesep 'freq_axis']);
load([spectral_dat_dir filesep 'time_axis']);
figure('units','normalized','outerposition',[0 0 1 1]);
for n = 1:n_chans
    if n_chans == 256
        p = plotGridPosition(grd(n)); 
    else
        p = plotGridPosition(grd(n), n_chans,ceil(sqrt(n_chans)));
    end
    subplot('position',p);
    
    if is_bad_chan(n)
        tstat_ch = zeros(length(time_axis), length(freq_axis));
    else
        load([spectral_dat_dir filesep 'spectral_erps_ch_' num2str(n) '.mat'])
        dat1 = spectral_erps_ch(:,:,is_class_1);
        dat2 = spectral_erps_ch(:,:,is_class_2);
        tstat_ch = (mean(dat1,3) - mean(dat2,3))./sqrt((var(dat1,[],3)./size(dat1,3)) + (var(dat2,[],3)./size(dat2,3)));
        % zero intermediate values:
        tstat_ch(abs(tstat_ch) < t_thresh) = 0;
    end
    
    % Get cmins and cmaxes:
    cmaxes(n) = max(1.5, max(tstat_ch(:)));
    cmins(n) = min(-1.5, min(tstat_ch(:)));
    %color_mapping = generate_red_blue_cmap(cmin, cmax, t_thresh)

    %imagesc(tstat_ch'); 
    imagesc(time_axis, freq_axis, tstat_ch'); 
    set(gca, 'YDir', 'normal');
    set(gca, 'YTick', []); set(gca, 'XTick', []);
    %colormap(color_mapping)
    hold on; plot([0 0], get(gca, 'YLim'), 'k');
    text(0,(max(get(gca,'YLim')) - (max(get(gca,'YLim')) * 0.1)),num2str(n));
end
%% Add Colormap:
[color_mapping, vals] = generate_red_blue_cmap(min(cmins), max(cmaxes), t_thresh);
for n = 1:n_chans
    if n_chans == 256
        p = plotGridPosition(grd(n)); 
    else
        p = plotGridPosition(grd(n), n_chans,ceil(sqrt(n_chans)));
    end
    ax = subplot('position',p);
    
    local_cmap = color_mapping((vals>=cmins(n)) & (vals<=cmaxes(n)),:);
    colormap(ax, local_cmap);
end

max_c = max(cmaxes)
min_c = min(cmins)
a = 1;
end

%% Helper Functions

%% ColorMap Function:
function [color_map, val_range] = generate_red_blue_cmap(cmin, cmax, white_thresh)
map_size = 101;
val_range = linspace(cmin, cmax, map_size);
max_red = abs(cmax)/max(abs([cmin cmax])); % percentage of saturation allowed in red
max_blue = abs(cmin)/max(abs([cmin cmax])); % percentage of saturation allowed in blue

base_intensity = 0.9; % Base intensity (white is 1)
saturation_val = 0.8; % highest color value permitted (relative to base)


color_map = base_intensity*ones(map_size,3);
% set red:
is_positive = (val_range > white_thresh);
color_map(is_positive,2) = linspace(base_intensity, (base_intensity-saturation_val*max_red), sum(is_positive));
color_map(is_positive,3) = linspace(base_intensity, (base_intensity-saturation_val*max_red), sum(is_positive));

% set blue:
is_negative = (val_range < -white_thresh);
color_map(is_negative,1) = linspace((base_intensity-saturation_val*max_blue), base_intensity, sum(is_negative));
color_map(is_negative,2) = linspace((base_intensity-saturation_val*max_blue), base_intensity, sum(is_negative));

end
