function ERPs = append_linguistic_stats(ERPs, lex_stats_dir)
%% This function loads and appends a structured array on lexical statistics
% onto the standard ERPs structured array used for analyzing the homophone experiment.
%
% This sub-strucuted array or linguistic information can be accessed as ERPs.lex_stats
%
% Inputs:
% EPRs - the standard data structure for the homophone experiment
%
% lex_stats_dir - the local directory that contains the file 'HomophoneLexicalStats.mat'
% - this file contains the lexical statistics for the stimuli of the homophone experiment, and 
% it must to saved locally for this funciton to operate.
%
% Output:
% ERPs - the standard homophone data structure, with the lex_stats
% structured array added on as an additional field
% 
% Ben Lucas - 1/11/16

%% This should have been the default:
ERPs.is_related_dominant = logical(ERPs.is_related_dominant);
ERPs.is_related_subordinant = logical(ERPs.is_related_subordinant);

%% Load Table of Lexical Statistics for homophones:
load([lex_stats_dir filesep 'HomophoneLexicalStats.mat']);
lex_stats = struct;

%% Relevant Flags:
is_dom = logical(ERPs.is_related_dominant);     % dominant trails
is_sub = logical(ERPs.is_related_subordinant);  % subordinate trials
is_rel = is_dom | is_sub;               % related prime-homophone trials
is_unrel = ~is_dom & ~is_sub;           % unrelated prime-homophone trials

words_table = HomophoneLexicalStats.Word;
homophones_exp = ERPs.target_names;
primes_exp = ERPs.prime_names;

%% Get Stats Relevant to Homophones:
dominance = zeros(size(homophones_exp));
freeassoc_strength = zeros(size(homophones_exp));
for k = 1:length(words_table)
    word = words_table{k};
    
    dominance_word = HomophoneLexicalStats.dominanceP1P2(k);
    freeassoc_word = HomophoneLexicalStats.FAscoreP1P2(k);
    
    dominance(strcmpi(homophones_exp, word)) = dominance_word;
    freeassoc_strength(strcmpi(homophones_exp, word)) = freeassoc_word;
end
lex_stats.dominance = dominance;
lex_stats.freeassoc_strength = freeassoc_strength;

is_dom = logical(ERPs.is_related_dominant);     % dominant trails
is_sub = logical(ERPs.is_related_subordinant);  % subordinate trials

is_high_fa_homophone = (freeassoc_strength > 0.5);
is_high_fa = (is_dom & is_high_fa_homophone) | (is_sub & ~is_high_fa_homophone);
is_low_fa = (is_sub & is_high_fa_homophone) | (is_dom & ~is_high_fa_homophone);


%% Get Stats Unique to Primes:
lsa_strength = zeros(size(primes_exp));
reaction_time = zeros(size(primes_exp));
prime_freeassoc_prob = zeros(size(primes_exp));
for k = 1:length(words_table)
    word = words_table{k};
    
    lsa_word = HomophoneLexicalStats.LSA_H_min_P(k);
    rt_word = HomophoneLexicalStats.meanRT_SPP(k);
    prime_assoc_word = HomophoneLexicalStats.FreeAssocationProbability(k);
    
    lsa_strength(strcmpi(primes_exp, word) & is_rel) = lsa_word;
    reaction_time(strcmpi(primes_exp, word) & is_rel) = rt_word;
    prime_freeassoc_prob(strcmpi(primes_exp, word) & is_rel) = prime_assoc_word;
end
lex_stats.lsa_strength = lsa_strength;
lex_stats.reaction_time = reaction_time;
lex_stats.prime_freeassoc_prob = prime_freeassoc_prob;

%%

ERPs.lex_stats = lex_stats;
%% ADD is_highfa and is_lowfa tags to main structure:
ERPs.is_high_fa = is_high_fa;
ERPs.is_low_fa = is_low_fa;
end
