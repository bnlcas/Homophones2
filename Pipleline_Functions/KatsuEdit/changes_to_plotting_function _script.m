%% Segments of code added to PlotingFunctions:


%% Determine is significant time points are included:
%(Necessary because there are now two kinds of variable inputs
sig_flag_ind = find(strcmpi(test, 'SigFlag'));
inclue_sig_timpts = false;
if ~isempty(sig_flag_ind)
    if length(sig_flag_ind) ~= 1
        warning('Multiple Significant TimePoint Matricies Given')
    end
    inclue_sig_timpts = true;
    sig_timepts_mat = varargin{sig_flag_ind(1) + 1};
    varargin(sig_flag_ind(1):(sig_flag_ind(1)+1)) = [];
end



%% Check for a field 'grd' for the ERP structure - use a grid permutation if it exists
if isfield(ERPs, 'grd')
    grd = ERPs.grd;
end


%% Adding Shading on each channel
if inclue_sig_timpts
    sig_timepts_ch = sig_timepts_mat(i,timerange);
    shaded_patch_significant_timepoints(time_axis, sig_timepts_ch)
    axis tight; set(gca,'YLim',[minbound maxbound]);
end
