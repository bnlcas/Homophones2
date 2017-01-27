function [] = combine_pdfs(root_dir)
% combines the pdfs for a given subject
% Assumes several folders with distinct png comparisons
%
% .png files must be converted into .pdf files in bash terminal: 
%% for i in $(ls *.png); do fne='.pdf'; sips -s format pdf $i --out $(basename $i .png)$fne; done


base_dir = pwd;
root_dir = '/Users/changlab/Documents/Homophone_Analysis/Spectrogrammar/EC125_spectra';
tmp = dir(root_dir);
tmp(1:2) = [];
for i = 1:length(tmp)
    cd([root_dir filesep tmp(i).name])
    pdf_list = {};
    sub_dir = [root_dir filesep tmp(i).name];
    figs = dir(sub_dir);
    for j = 1:length(figs)
        filename = strsplit(figs(j).name,'.');
        if ~strcmp(filename{1},'') & strcmp(filename{2},'pdf')
            pdf_list = [pdf_list figs(j).name];
        end
    end
    out_name = [sub_dir filesep tmp(i).name '_COMBINED_Figures.pdf'];
    append_pdfs(out_name, pdf_list{:})
    
end
   
cd(base_dir)

end
