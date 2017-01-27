function [] = plot_erp_comparison_select_electrodes(time_axis, ecog1, ecog2, select_electrodes, is_sig_mat, anatomy_elecs)
%% This is a specialized function designed to plot the ERPs for two ecog grids
% in a subset of electrodes;
% by default the function will attempt to plot a square grid as best a possible

% Inputs:
% time_axis - the time range of the second dimension of ecog1 &2 (must
% match the size of second dimension of ecog1&2 (these must be scaled prior
% to input)
%
% ecog1, ecog2 - 256xtimeptsxn ecog matricies (must be clipped on time
% axis)
%
% select_electrodes - list of indecies (of 1:256), that are significant
%
% is_sig_mat - 256xtimepts boolean matrix for whether timeoints show
% significant difference between ecog1 and ecog2
%
% (OPTIONAL):
% anatomy_elecs = 256x1 list of the anatomical region of each electrode

if ~exist('anatomy_elecs', 'var') || isempty(anatomy_elecs)
    anatomy_elecs = repmat({'brain'},256,1);
end

numplots = length(select_electrodes);
grid_dim = ceil(sqrt(numplots));

mean_mat = cat(3, mean(ecog1(select_electrodes,:,:),3), mean(ecog2(select_electrodes,:,:),3));
sem_mat = cat(3, nansem(ecog1(select_electrodes,:,:),3), nansem(ecog2(select_electrodes,:,:),3));

maxbound = max(mean_mat(:) + sem_mat(:));
minbound = min(mean_mat(:) - sem_mat(:));







figure;
electrode_colorlist = [0.8 0 0; 0 0 0.8;0.8 0.8 0];
anatomy_colorlist = [0.8 0 0; 0 0 0.8; 0.8 0.8 0; 0 0.8 0.8; 0 0.8 0; 0.8 0 0.8; 0.2 0.6 0; 0.6 0.2 0; 0 0.6 0.2];    % Plots of colors listed as Red, Blue, Yellow, Teal, Green, Purple
grd = reshape(1:(grid_dim^2),grid_dim,grid_dim); % grid plot...
anatomy_regions = unique(anatomy_elecs(select_electrodes));
if length(anatomy_regions) == 1
    anatomy_colorlist = [0.2 0.2 0.2];
end
for i = 1:numplots
    anatomy_region = find(strcmpi(anatomy_regions, anatomy_elecs{select_electrodes(i)}));

    p = plotGridPosition(grd(i), grid_dim^2, grid_dim); 
    subplot('position',p);
    %subplot(grid_dim, grid_dim, i)
    hold on;

    %% Plot ERPs for Channel
    for j = 1:2
        shadedErrorBar(time_axis,mean_mat(i,:,j), sem_mat(i,:,j), {'color', electrode_colorlist(j,:)},1);
    end
    
    axis tight;
    set(gca,'YLim',[minbound maxbound]...
          ,'XTickLabel',[],'YTickLabel',[],'YTick',[minbound maxbound]);
        line([0 0],get(gca,'YLim'),'Color','k');
        line(get(gca,'XLim'),[0 0],'Color','k');
    text(0,(max(get(gca,'YLim')) - (max(get(gca,'YLim')) * 0.1)),num2str(select_electrodes(i)), 'Color', anatomy_colorlist(anatomy_region,:));
        

    
    %% Add shaded patches:
    sig_ch = is_sig_mat(select_electrodes(i),:);
    sig_ch(1) = false;
    sig_ch(end) = false; % boundary conditions - simplify the situation
    onsets = find(diff(sig_ch) == 1);
    offsets = find(diff(sig_ch) == -1);
    if ~isempty(onsets)
        for j = 1:length(onsets)
            alph = 0.1;
            patch([time_axis(onsets(j)), time_axis(offsets(j))+0.005, time_axis(offsets(j))+0.005, time_axis(onsets(j))], [minbound minbound maxbound maxbound],'black', 'FaceAlpha',alph,'EdgeAlpha',alph)
        end
    end

    %% color plot axis for anatomy:
    box on;
    ax = gca;
    ax.XColor = anatomy_colorlist(anatomy_region,:);
    ax.YColor = anatomy_colorlist(anatomy_region,:);
 
end
    
    
    
    
    
    
end
