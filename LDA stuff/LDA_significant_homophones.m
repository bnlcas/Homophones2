function [lda_accuracies, lda_deviances] = LDA_significant_homophones(ERPs)
%% The propose of this function to find significant electrodes encoding differences
% in the primes of each homophone and running an LDA classifier on this
% data.


%% Operator Flag to stack adjacent time points for classifier
stack_data = false; % aggregate data over a time window for the purpose of training a classifier
window_data = true;
if stack_data
    stack_size = 10;
end
if window_data
    window_size = 10; % number of time points to smooth data by.
end

%% First plot primes as a matter of semantic related need
is_good = ERPs.is_good_trial;
is_dom = (ERPs.is_related_dominant == 1) & is_good;
is_sub = (ERPs.is_related_subordinant == 1) & is_good;
is_rel = (is_dom | is_sub) & is_good;
is_unr = ~is_dom & ~is_sub & is_good;

%% Set Time Range of Analysis:
time_range = [-1 1]; % in seconds
[~,time_win(1)] = min(abs(ERPs.time_axis - time_range(1)));
[~,time_win(2)] = min(abs(ERPs.time_axis - time_range(2)));


%% Second, plot the relatedness for each homophone
homophones = unique(ERPs.target_names);
lda_accuracies = zeros(length(time_win(1):time_win(2)), length(homophones));
lda_deviances = zeros(length(time_win(1):time_win(2)), length(homophones));

for i = 1:length(homophones)
%i = 4;
    homophone = homophones{i};
    is_homo = strcmpi(ERPs.target_names, homophone);

    
    %% Compare P1 P2
    %ecog_p1 = ERPs.ecog_targets(:,time_win(1):time_win(2),is_homo&is_dom);
    %ecog_p2 = ERPs.ecog_targets(:,time_win(1):time_win(2), is_homo&is_sub);
    
    ecog_p1 = ERPs.ecog_targets(:,time_win(1):time_win(2),is_homo&is_rel);
    ecog_p2 = ERPs.ecog_targets(:,time_win(1):time_win(2), is_homo&is_unr);
    time_axis = ERPs.time_axis(time_win(1):time_win(2));
    
    % zero Bad Trials
    ecog_p1(ERPs.BadChans,:,:) = 0;
    ecog_p2(ERPs.BadChans,:,:) = 0;
    
    if window_data
        ecog_p1 = smooth_3d(ecog_p1, window_size, 2);
        ecog_p2 = smooth_3d(ecog_p2, window_size, 2);
    end
    
    
    %% Run through each time point and channel to make statistical comparison for significance
    is_sig_ttest = false(size(ecog_p1,1), size(ecog_p1,2));
    stat_tstat = zeros(size(ecog_p1,1), size(ecog_p1,2));
    
    pval = 0.01;
    %pval = 0.01/100;
    bf_correct = false; % adjust p-thresh by a boniferroni correction of # of channels
    if bf_correct
        pval = pval/(size(ecog_p1,1) - length(ERPs.BadChans));
    end
    
    for j = 1:size(ecog_p1,1)
        if (sum(j == ERPs.BadChans) == 0) % ignore badchans (t-test gives nans)
            for k = 1:size(ecog_p1,2)
                dat_p1 = squeeze(ecog_p1(j,k,:));
                dat_p2 = squeeze(ecog_p2(j,k,:));

                %% T-Test
                [is_sig_ttest(j,k),~,~,tmp] = ttest2(dat_p1, dat_p2, 'Alpha', pval);
                stat_tstat(j,k) = abs(tmp.tstat);
             end
        end
    end
    
    %% Plot the electrodes that are significant
    sig_timepts_thresh = 1; % minimum number of significant time points for a channel to be significant
    sig_window = 50:size(is_sig_ttest,2); % -0.5 to 1 second window
    sig_ch_t = find(sum(is_sig_ttest(:,sig_window ),2) >= sig_timepts_thresh);    

    
    %% Add Neal Speech Responsive
    is_speech_responsive = get_speech_responsive_elecs(ERPs);
    sig_ch_t = intersect(sig_ch_t, find(is_speech_responsive));
    
    plot_sig_elects = false;
    if plot_sig_elects
        plot_erp_comparison_select_electrodes(time_axis, ecog_p1, ecog_p2, sig_ch_t, is_sig_ttest)
    end
    
    %% Hand Select Electrodes?:
    %sig_ch_t = [93 220];  [147 220];
    
    p_labels = [repmat({'p1'},size(ecog_p1,3),1); repmat({'p2'},size(ecog_p2,3),1)];
    if stack_data
        p_labels = repmat(p_labels,stack_size,1);
    end
    
    
    %% Generate PCA matrix for (Num_Chans x (time_pts x trials))
    dat_pca = cat(3,ecog_p1(sig_ch_t,:,:), ecog_p2(sig_ch_t,:,:));
    dat_pca = reshape(dat_pca, size(dat_pca,1), size(dat_pca,2)*size(dat_pca,3)); % flatten time x trials
    dat_pca = dat_pca - repmat(mean(dat_pca,2),1,size(dat_pca,2)); % center data
    [pca_mat,~,~,~,exp] = pca(dat_pca');
    
% % % %     pc1 = pca_mat(:,1);
% % % %     elects = zeros(256,1);
% % % %     elects(sig_ch_t) = pc1;
% % % %     sig_grid = rot90(reshape(abs(elects),16,16),2);
% % % %     figure; imagesc(sig_grid); title('Electrode Weigths of PC 1')
    num_components = 10;    % number of components to use in PCA projection
%    explained = sum(exp(1:num_components));
    
    if window_data
        ecog_p1 = smooth_3d(ecog_p1, window_size, 2);
        ecog_p2 = smooth_3d(ecog_p2, window_size, 2);
    end

    %% Run Classifier 
    %% For each time point construct a classifier & Calculate accuracy
    for j = 1:length(time_axis)
        if ~stack_data
            dat1 = squeeze(ecog_p1(sig_ch_t,j,:)); % transpose necessary later
            dat2 = squeeze(ecog_p2(sig_ch_t,j,:)); % num_ch x num_trials matrix of data for the time point
            dat_comb = [dat1'; dat2'];
        else
            dat_comb = [];
            stack_inds = (j+1-ceil(stack_size/2)):(j+floor(stack_size/2)); % adjacent time points
            stack_inds(stack_inds > length(time_axis)) = length(time_axis); % set points at edge be in bounds
            stack_inds(stack_inds < 1) = 1;                                 % Kludgy, depends on significanct of end points
            for shift = 1:length(stack_inds)
                dat1 = squeeze(ecog_p1(sig_ch_t,stack_inds(shift),:));
                dat2 = squeeze(ecog_p2(sig_ch_t,stack_inds(shift),:));
                dat_comb = [dat_comb; dat1'; dat2']; % growing array bad, fix - necessary to interleave for CV partition
            end
        end
        

        dat_comb_pca = gsubtract(dat_comb, mean(dat_comb,1));%*pca_mat(:,1:num_components);  % PCA projection of combined data    
%        lda_timept = fitcdiscr(dat_comb_pca, p_labels);
%        class_predict = predict(lda_timept, dat_comb_pca);

        folds = 12;          % number of folds in K-Fold CV
        if stack_data
            partition = cvpartition(p_labels(1:(length(p_labels)/stack_size)),'KFold',folds); % partition the first repetition of labels, and repeat test/train later
        else
            partition = cvpartition(p_labels,'KFold',folds);
            %partition = cvpartition(length(p_labels),'LeaveOut');
        end
        %partition = cvpartition(p_labels, 'Holdout', 0.3);
        %% Tune model parameters:
% %         inds = randperm(length(p_labels),floor(length(p_labels)/2));
% %         lda_timept = fitcdiscr(dat_comb_pca(inds,:), p_labels(inds));
% %         [err,gamma,delta,numpred] = cvshrink(lda_timept, 'NumGamma',num_components,'NumDelta',24);
% %         [~,min_ind] = min(err(:).*numpred(:)./(err(:)+numpred(:)));
% %         [mi, mj] = ind2sub([25 25], min_ind);
% %         lda_timept.Gamma = gamma(mj);
% %         lda_timept.Delta = delta(mi,mj);
        
        folds = partition.NumTestSets;
        accuracies = zeros(folds,1);
        for k = 1:folds
            %% Parition Data:
            train_inds = training(partition,k);
            test_inds = test(partition,k);
            if stack_data
                train_inds = repmat(train_inds,stack_size,1); % complicated: was necessary to prevent correlated test/train data
                test_inds = repmat(test_inds, stack_size,1); 
            end
            %% LDA: 
            %lda_timept = fitcdiscr(dat_comb_pca(train_inds,:), p_labels(train_inds));
            svm_timept = fitcsvm(dat_comb_pca(train_inds,:), p_labels(train_inds),'KernelFunction','linear', 'BoxConstraint', 0.0001, 'KernelScale', 1, 'prior','uniform');
% % % %             [err,gamma,delta,numpred] = cvshrink(lda_timept, 'NumGamma',num_components,'NumDelta',25);
% % % %             is_num_pred = (numpred<num_components);
% % % %             min_err = min(min(err(is_num_pred)));
% % % %             [a,b] = find(err == min_err);
% % % %             lda_timept.Gamma = gamma(a(1));
% % % %             lda_timept.Delta = delta(a(1),b(1));
% %         [~,min_ind] = min(err(:).*numpred(:)./(err(:)+numpred(:)));
% 
            
            %class_predict = predict(lda_timept, dat_comb_pca(test_inds,:));
            class_predict = predict(svm_timept, dat_comb_pca(test_inds,:));
                        
            accuracies(k) = sum(strcmpi(class_predict, p_labels(test_inds)))/sum(test_inds);
        end
        %lda_timept = fitcdiscr(dat_comb_pca, p_labels,'KFold',10);
        %class_predict = kfoldPredict(lda_timept);
        %lda_accuracies(j,i) = sum(strcmpi(class_predict, p_labels))/length(p_labels);
        lda_accuracies(j,i) = mean(accuracies); % mean accuray of the k-folds for homophone i, in time point j
        lda_deviances(j,i) = std(accuracies); % std accuracy of k-folds for ...
    end
    
    
     
 %% Plot time courses (use 'LETTER' as test case)
    
%end

%figure; plot(repmat(time_axis',1,length(homophones))', lda_accuracies');
%legend(homophones)
% figure; plot(time_axis, lda_accuracies(:,4))
% legend('Classifier Accuracy'); title(['LDA Classification for ' homophones{4}])

figure;
%gcf = h;
%subplot(1,2,2)
shadedErrorBar(time_axis, lda_accuracies(:,i), lda_deviances(:,i)); % mean class of k-
hold on
plot(get(gca,'XLim'), [0.5 0.5],'k--');
plot([0 0], get(gca,'YLim'), 'k')
axis tight;
xlabel('time (s)');
ylabel('classification accuracy');

title([num2str(folds) '-Fold Classification Accuracy for ' upper(homophones{i})])
%title(['K-Fold LDA Classification for ' upper(homophones{4})])
%title({['K-Fold LDA Classification Accuracy of 12-Component PCA for Dom/Sub']; ['Senses for Prime Locked ERPs of ' upper(homophones{4}) ' with 10 pt Stacking']})

%title({['K-Fold LDA Classification Accuracy of 6-Component PCA']; ['for Dom/Sub Senses for Target Locked ERPs of ' upper(homophones{4})]});% ' with 10 pt Stacking']})
end


figure;
subplot(2,2,[1,2]);
shadedErrorBar(time_axis, lda_accuracies(:,i), lda_deviances(:,i)); % mean class of k-
hold on
plot(get(gca,'XLim'), [0.5 0.5],'k--');
plot([0 0], get(gca,'YLim'), 'k')
ylim([0 1.2])
xlabel('time (s)')
ylabel('classification accuracy')

title(['K-Fold SVM Classification for ' upper(homophones{i}) ' for Target Locked Senses in CH ' num2str(sig_ch_t(1)) ' & CH ' num2str(sig_ch_t(2))])

subplot(2,2,3)
shadedErrorBar(time_axis, mean(squeeze(ecog_p1(sig_ch_t(1),:,:)),2),nansem(squeeze(ecog_p1(sig_ch_t(1),:,:)),2),'r',1);
hold on;
shadedErrorBar(time_axis, mean(squeeze(ecog_p2(sig_ch_t(1),:,:)),2),nansem(squeeze(ecog_p1(sig_ch_t(1),:,:)),2),'b',1);
plot([0 0], get(gca,'YLim'), 'k')
plot(get(gca,'XLim'), [0 0 ],'k');
xlabel('time (s)')
ylabel('High Gamma Intensity')
title(['High Gamma Activity in CH ', num2str(sig_ch_t(1))])



subplot(2,2,4)
shadedErrorBar(time_axis, mean(squeeze(ecog_p1(sig_ch_t(2),:,:)),2),nansem(squeeze(ecog_p1(sig_ch_t(2),:,:)),2),'r',1);
hold on;
shadedErrorBar(time_axis, mean(squeeze(ecog_p2(sig_ch_t(2),:,:)),2),nansem(squeeze(ecog_p1(sig_ch_t(2),:,:)),2),'b',1);
plot([0 0], get(gca,'YLim'), 'k')
plot(get(gca,'XLim'), [0 0 ],'k');
xlabel('time (s)')
ylabel('High Gamma Intensity')
title(['High Gamma Activity in CH ', num2str(sig_ch_t(2))])

a = 1;



end