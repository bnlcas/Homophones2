function [is_diff_ci, bs_median] = bootstrap_difference(set1, set2, varargin)
%% This function takes two sets of data and will run a boot strap to
% calculate the median difference bewteen these two data sets,
% in addition it will retun the bootstrapped confidence interval for these
% two data sets
%
% Inputs: 
% set1, set2 - two 1-d arrays of data
%
% Variable Inputs:
%
% input1: # of iterations (default 1000)
%
% input2: confidence interval (number 0-1) - default 0.95.

iterations = 1000;
ci = 0.95;

if length(varargin) > 0
    if varargin{1}>0
        iterations = varargin{1};
    end
end
if length(varargin) > 1
    if varargin{2}<=1 & varargin{2} >=0
        ci = varargin{2};
    end
end


%% bootstrap ci:
differences = zeros(1, iterations);
leng_set1 = length(set1);
leng_set2 = length(set2);
for i = 1:iterations
    differences(i) = set1(randperm(leng_set1,1)) - set2(randperm(leng_set2,1));
end
differences = sort(differences);

%% return output
bs_median = median(differences);
ci_ind(1) = ceil(iterations*(1-ci)/2);
ci_ind(2) = ceil(iterations*(1-(1-ci)/2));
bs_ci = [differences(ci_ind(1)), differences(ci_ind(2))];


if bs_ci(1) > 0 | bs_ci(2) < 0
    is_diff_ci = true;
else
    is_diff_ci = false;
end

end
    