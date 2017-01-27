spectral_erps = spectral_erps_rel - spectral_erps_unrel;
%time_axis
%freq_axis
sub_plot_coords = position_subplot_grid(size(spectral_erps,1));
figure;
for k = 1:size(spectral_erps,1);
    subplot('position', sub_plot_coords(k,:)); hold on;
    ecog_ch = squeeze(spectral_erps(k,:,:));
    imagesc(time_axis, freq_axis, ecog_ch'); axis tight;
    set(gca,'YDir', 'normal')
    hold on;
    plot([0 0], get(gca,'YLim'),'r')
    set(gca, 'XTick', []);
    set(gca, 'YTick',[]);
    axis tight;
end
