function repeat_type = get_prompt_repeat_type(evnt)
%% This function takes an evnt structure and for each pair of events
% determines whether the prompt was to repeat event 1 or event 2

%exp_data_dir = '/Users/changlab/Documents/data/EC123/behavior';
tmp = strsplit(evnt(1).dpath,'/');
exp_data_dir = [strjoin(tmp(1:(end-2)),'/'), '/behavior'];
block = evnt(1).block;
subj = evnt(1).subject;

if exist([exp_data_dir '/' subj '_' block '.mat']) == 2
    load([exp_data_dir '/' subj '_' block '.mat']);
elseif exist([exp_data_dir '/Homophone_EXP_' subj '_' block '.mat']) == 2
    load([exp_data_dir '/Homophone_EXP_' subj '_' block '.mat']);
end
prompts = [exp_dat.primes exp_dat.targets exp_dat.cues];


num_trials = length(evnt)/2;
repeat_type = zeros(num_trials,1);
for i = 1:num_trials
    prime = evnt(2*i-1).name;
    target = evnt(2*i).name;
    pair_inds = find(strcmpi(prompts(:,1),prime) & strcmpi(prompts(:,2),target));
    [~,ind] = min(abs(pair_inds-i));
    prompt = prompts(pair_inds(ind),:);

    if strcmpi(prompt(3),prompt(2))
        repeat_type(i) = 2;
    else
        repeat_type(i) = 1;
    end
end


end
    