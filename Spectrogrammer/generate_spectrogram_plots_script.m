function [] = generate_spectrogram_plots_script(ERPs,spectral_dat_dir,save_fig_dir)
% Plots a standard set of Homophone Comparisons
%
% EPRs - any erp structure for a subject (is necessary to get the proper
% timing data on stimuli)
%
% spectral_dat_dir = '/Users/changlab/Documents/data/EC133/ecog_data/spectrogram_data';
% 
mkdir(save_fig_dir)
is_rel = ERPs.is_good_trial & (ERPs.is_related_dominant | ERPs.is_related_subordinant);
homophones = unique(ERPs.target_names(is_rel));

%% High vs Low FA:
HighLowFA_dir = [save_fig_dir filesep 'HighFA_vs_LowFA_Spectrograms'];
mkdir(HighLowFA_dir);
is_highFA = ERPs.is_good_trial & ERPs.is_high_fa;
is_lowFA = ERPs.is_good_trial & ERPs.is_low_fa;
plot_spectral_tstat(ERPs, spectral_dat_dir, is_highFA, is_lowFA);
    drawnow
    pause(15)
    set(gcf, 'PaperPositionMode', 'auto')
    annotation('textbox', 'String', ['HighFA Primes (n = ' num2str(sum(is_highFA)) ') vs LowFA Primes (n = ' num2str(sum(is_lowFA)) ') for ALL Homophones'], ...
        'Position', [0, 0.98, 1, 0.02], 'HorizontalAlignment','center','LineStyle', 'none');
    saveas(gcf,[HighLowFA_dir filesep 'HighFA_vs_LowFA_ALL.png']);
    
    close gcf
    

for k = 1:length(homophones)
    [homophones{k} ':']
    is_homo = strcmpi(ERPs.target_names, homophones{k});
    is_highFA_homo = is_highFA & is_homo;
    is_lowFA_homo = is_lowFA & is_homo;
    plot_spectral_tstat(ERPs, spectral_dat_dir, is_highFA_homo, is_lowFA_homo)
    
    drawnow
    pause(10)
    set(gcf, 'PaperPositionMode', 'auto')
    annotation('textbox', 'String', ['HighFA Prime (' ERPs.prime_names{find(is_highFA_homo,1)} '; n = ' num2str(sum(is_highFA_homo)) ') vs LowFA Prime (' ERPs.prime_names{find(is_lowFA_homo,1)} '; n = ' num2str(sum(is_lowFA_homo)) ') for ' upper(homophones{k})], ...
        'Position', [0, 0.98, 1, 0.02], 'HorizontalAlignment','center','LineStyle', 'none');

    saveas(gcf,[HighLowFA_dir filesep 'HighFA_vs_LowFA_' upper(homophones{k}) '.png']);
    close gcf
    a = 1;
end

%% Related vs Unrelated:
relunrel_dir = [save_fig_dir filesep 'Related_vs_Unrelated_Spectrograms'];
mkdir(relunrel_dir)
is_rel = ERPs.is_good_trial & (ERPs.is_high_fa | ERPs.is_low_fa);
is_unrel = ERPs.is_good_trial & ~(ERPs.is_high_fa | ERPs.is_low_fa);
plot_spectral_tstat(ERPs, spectral_dat_dir, is_rel, is_unrel);
    drawnow
    pause(10)
    set(gcf, 'PaperPositionMode', 'auto')
    annotation('textbox', 'String', ['Related Primes (n = ' num2str(sum(is_rel)) ') vs Unrelated Primes (n = ' num2str(sum(is_unrel)) ') for ALL Homophones'], ...
        'Position', [0, 0.98, 1, 0.02], 'HorizontalAlignment','center','LineStyle', 'none');

    saveas(gcf,[relunrel_dir filesep 'Rel_vs_Unrel_ALL.png']);
    close gcf
for k = 1:length(homophones)
    is_homo = strcmpi(ERPs.target_names, homophones{k});
    is_rel_homo = is_rel & is_homo;
    is_unrel_homo = is_unrel & is_homo;
    plot_spectral_tstat(ERPs, spectral_dat_dir, is_rel_homo, is_unrel_homo)
    drawnow
    pause(10)
    set(gcf, 'PaperPositionMode', 'auto')
        annotation('textbox', 'String', ['Related Primes (n = ' num2str(sum(is_rel_homo)) ') vs Unrelated Primes (n = ' num2str(sum(is_rel_homo)) ') for ' upper(homophones{k})], ...
        'Position', [0, 0.98, 1, 0.02], 'HorizontalAlignment','center','LineStyle', 'none');
    saveas(gcf,[relunrel_dir filesep 'Rel_vs_Unrel_' upper(homophones{k}) '.png']);
    close gcf
    homophones{k}
    a=1;
end

%% dominant vs subordinate
domsub_dir = [save_fig_dir filesep 'Dominant_vs_Subordinate_Spectrograms'];
mkdir(domsub_dir)

is_dom = ERPs.is_good_trial & ERPs.is_related_dominant;
is_sub = ERPs.is_good_trial & ERPs.is_related_subordinant;
plot_spectral_tstat(ERPs, spectral_dat_dir, is_dom, is_sub);
    drawnow
    pause(10)
    set(gcf, 'PaperPositionMode', 'auto')
        annotation('textbox', 'String', ['Dominant Primes (n = ' num2str(sum(is_dom)) ') vs SubordinatePrimes (n = ' num2str(sum(is_sub)) ') for ALL Homophones'], ...
        'Position', [0, 0.98, 1, 0.02], 'HorizontalAlignment','center','LineStyle', 'none');
    saveas(gcf,[domsub_dir filesep 'Dom_vs_Sub_ALL.png']);
    close gcf
for k = 1:length(homophones)
    ['DomSub ' homophones{k} ':']
    is_homo = strcmpi(ERPs.target_names, homophones{k});
    is_dom_homo = is_dom & is_homo;
    is_sub_homo = is_sub & is_homo;
    plot_spectral_tstat(ERPs, spectral_dat_dir, is_dom_homo, is_sub_homo)
        drawnow
    pause(10)
    set(gcf, 'PaperPositionMode', 'auto')
      annotation('textbox', 'String', ['Dominant Prime (' ERPs.prime_names{find(is_dom_homo,1)} '; n = ' num2str(sum(is_dom_homo)) ') vs Subordinate Prime (' ERPs.prime_names{find(is_sub_homo,1)} '; n = ' num2str(sum(is_sub_homo)) ') for ' upper(homophones{k})], ...
        'Position', [0, 0.98, 1, 0.02], 'HorizontalAlignment','center','LineStyle', 'none');

    saveas(gcf,[domsub_dir filesep 'Dom_vs_Sub_' upper(homophones{k}) '.png']);
    close gcf

    a=1;
end


