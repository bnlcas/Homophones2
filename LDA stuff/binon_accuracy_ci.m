function ci_out = binon_accuracy_ci(n,p,alpha)
% function returns the Confidence interval 
% for a classifer on n trials with a probabily of p being correct
% to a confience level alpha
% EX: n = 27; p = 0.5; alpha = 0.05
% ci_out = binon_accuracy_ci(n,p,alpha);
x_cdf = binocdf(1:n, n, p);
ci(1) = find(x_cdf >= alpha/2,1);
ci(2) = find(x_cdf <= (1-alpha/2),1, 'last');
ci_out = ci/n;

end