function [] = generate_erp_plots_script(ERPs, save_fig_dir)
% Plots a standard set of Homophone Comparisons for an given ERP strucutre
% Uses target locked ERPs and a [-1 1] time window
%
% EPRs - any erp structure for a subject (is necessary to get the proper
% timing data on stimuli)
%
% directory to save figures into
% /data/EC131/basic_erp_comparisons
% 
mkdir(save_fig_dir)
is_rel = ERPs.is_good_trial & (ERPs.is_related_dominant | ERPs.is_related_subordinant);
homophones = unique(ERPs.target_names(is_rel));

erp_twin = [-1 1];
ecog = ERPs.ecog_targets;

%% High vs Low FA:
HighLowFA_dir = [save_fig_dir filesep 'HighFA_vs_LowFA_ERPs'];
mkdir(HighLowFA_dir);
is_highFA = ERPs.is_good_trial & ERPs.is_high_fa;
is_lowFA = ERPs.is_good_trial & ERPs.is_low_fa;
is_null = ERPs.is_good_trial & ~is_highFA & ~is_lowFA;

is_sigpts = ttest2(ecog(:,:,is_highFA), ecog(:,:,is_lowFA), 'Dim',3, 'Alpha', 0.05);
fig_title = ['HighFA Primes (Red; n = ' num2str(sum(is_highFA)) ') vs LowFA Primes (Blue; n = ' num2str(sum(is_lowFA)) ') for ALL Homophones'];

PlotECogGrid_std_v3(ERPs, erp_twin, ecog(:,:,is_highFA), ecog(:,:,is_lowFA), ecog(:,:,is_null), ...
    'SigPts', is_sigpts, 'Title', fig_title);
    drawnow
    pause(5)
    set(gcf, 'PaperPositionMode', 'auto')

    saveas(gcf,[HighLowFA_dir filesep 'HighFA_vs_LowFA_ALL.png']);
    
    close gcf
    

for k = 1:length(homophones)
    [homophones{k} ':']
    is_homo = strcmpi(ERPs.target_names, homophones{k});
    is_highFA_homo = is_highFA & is_homo;
    is_lowFA_homo = is_lowFA & is_homo;
    is_null_homo = is_null & is_homo;
    
    is_sigpts = ttest2(ecog(:,:,is_highFA_homo), ecog(:,:,is_lowFA_homo), 'Dim',3, 'Alpha', 0.05);
    fig_title = ['HighFA Prime (' ERPs.prime_names{find(is_highFA_homo,1)} '; Red; n = ' num2str(sum(is_highFA_homo)) ') vs LowFA Prime (' ERPs.prime_names{find(is_lowFA_homo,1)} '; blue; n = ' num2str(sum(is_lowFA_homo)) ') for ' upper(homophones{k})];

    PlotECogGrid_std_v3(ERPs, erp_twin, ecog(:,:,is_highFA_homo), ecog(:,:,is_lowFA_homo), ecog(:,:,is_null_homo), ...
        'SigPts', is_sigpts, 'Title', fig_title);

    drawnow
    pause(5)
    set(gcf, 'PaperPositionMode', 'auto')

    saveas(gcf,[HighLowFA_dir filesep 'HighFA_vs_LowFA_' upper(homophones{k}) '.png']);
    close gcf
    a = 1;
end





%% Related vs Unrelated:
relunrel_dir = [save_fig_dir filesep 'Related_vs_Unrelated_ERPs'];
mkdir(relunrel_dir)
is_rel = ERPs.is_good_trial & (ERPs.is_high_fa | ERPs.is_low_fa);
is_unrel = ERPs.is_good_trial & ~(ERPs.is_high_fa | ERPs.is_low_fa);

is_sigpts = ttest2(ecog(:,:,is_rel), ecog(:,:,is_unrel), 'Dim',3, 'Alpha', 0.05);
fig_title = ['Related Primes (Red; n = ' num2str(sum(is_rel)) ') vs Unrelated Primes (Blue; n = ' num2str(sum(is_unrel)) ') for ALL Homophones'];

PlotECogGrid_std_v3(ERPs, erp_twin, ecog(:,:,is_rel), ecog(:,:,is_unrel), ...
    'SigPts', is_sigpts, 'Title', fig_title);
    drawnow
    pause(5)
    set(gcf, 'PaperPositionMode', 'auto')

    saveas(gcf,[relunrel_dir filesep 'Rel_vs_Unrel_ALL.png']);
    close gcf
for k = 1:length(homophones)
    ['Rel-UnRel ' homophones{k} ':']
    is_homo = strcmpi(ERPs.target_names, homophones{k});
    is_rel_homo = is_rel & is_homo;
    is_unrel_homo = is_unrel & is_homo;
    
    is_sigpts = ttest2(ecog(:,:,is_rel_homo), ecog(:,:,is_unrel_homo), 'Dim',3, 'Alpha', 0.05);
    fig_title = ['Related Primes (Red; n = ' num2str(sum(is_rel_homo)) ') vs Unrelated Primes (Blue; n = ' num2str(sum(is_unrel_homo)) ') for ' upper(homophones{k})];

    PlotECogGrid_std_v3(ERPs, erp_twin, ecog(:,:,is_rel), ecog(:,:,is_unrel), ...
         'SigPts', is_sigpts, 'Title', fig_title);

    drawnow
    pause(5)
    set(gcf, 'PaperPositionMode', 'auto')
    saveas(gcf,[relunrel_dir filesep 'Rel_vs_Unrel_' upper(homophones{k}) '.png']);
    close gcf
    homophones{k}
    a=1;
end

%% dominant vs subordinate
domsub_dir = [save_fig_dir filesep 'Dominant_vs_Subordinate_ERPs'];
mkdir(domsub_dir)

is_dom = ERPs.is_good_trial & ERPs.is_related_dominant;
is_sub = ERPs.is_good_trial & ERPs.is_related_subordinant;
is_null = ERPs.is_good_trial & ~is_dom & ~is_sub;

is_sigpts = ttest2(ecog(:,:,is_dom), ecog(:,:,is_sub), 'Dim',3, 'Alpha', 0.05);
fig_title = ['Dominant Primes (Red; n = ' num2str(sum(is_dom)) ') vs Subordinate Primes (Blue; n = ' num2str(sum(is_sub)) ') for ALL Homophones'];

PlotECogGrid_std_v3(ERPs, erp_twin, ecog(:,:,is_dom), ecog(:,:,is_sub), ecog(:,:,is_null),...
             'SigPts', is_sigpts, 'Title', fig_title);

drawnow
pause(5)
set(gcf, 'PaperPositionMode', 'auto')
saveas(gcf,[domsub_dir filesep 'Dom_vs_Sub_ALL.png']);
close gcf

for k = 1:length(homophones)
    ['DomSub ' homophones{k} ':']
    is_homo = strcmpi(ERPs.target_names, homophones{k});
    is_dom_homo = is_dom & is_homo;
    is_sub_homo = is_sub & is_homo;
    is_null_homo = is_null & is_homo;
    
    is_sigpts = ttest2(ecog(:,:,is_dom_homo), ecog(:,:,is_sub_homo), 'Dim',3, 'Alpha', 0.05);
    is_sigpts(ERPs.BadChans
    fig_title = ['Dominant Primes (Red; n = ' num2str(sum(is_dom_homo)) ') vs Subordinate Primes (Blue; n = ' num2str(sum(is_sub_homo)) ') for ' upper(homophones{k})];

    PlotECogGrid_std_v3(ERPs, erp_twin, ecog(:,:,is_dom_homo), ecog(:,:,is_sub_homo), ecog(:,:,is_null_homo),...
       'SigPts', is_sigpts, 'Title', fig_title);
    
    drawnow
    pause(10)
    set(gcf, 'PaperPositionMode', 'auto')
    saveas(gcf,[domsub_dir filesep 'Dom_vs_Sub_' upper(homophones{k}) '.png']);
    close gcf

    a=1;
end


