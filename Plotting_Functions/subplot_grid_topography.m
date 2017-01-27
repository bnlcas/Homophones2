function [ all_subplot_nums , each_subplot_num ] = subplot_grid_topography(subplot_rows,subplot_cols,location_and_direction,grid_ordering)
% to plot subplots as oriented on the brain ...
% INPUTS:
% subplot_rows ............. positive integer corresponding to the total number of ROWS in the figure of subplots (e.g., 16)
% subplot_cols ............. positive integer corresponding to the total number of COLUMNS in the figure of subplots (e.g., 16)
% location_and_direction ... three-character string:
% % % % ... (1) t/b = location of first elec is Top/Bottom of grid
% % % % ... (2) l/r = location of first elec is Left/Right of grid
% % % % ... (3) t/b = location of last elec in the same strip as the first elec (e.g., #16 in a 16x16 grid) is Top/Bottom of grid
% % % % % % % % % % % % ... [indicates direction of the flow of numbering from first elec]
% % location_and_direction can be one of the following: { 'tlt' 'tlb' 'trt' 'trb' 'blb' 'blt' 'brb' 'brt' } (default is 'tlt')
% grid_ordering ............ cell array of size num_minigrids x  4 with columns corresponding to:
% % % % ... (1) num_subplot_rows in a given minigrid
% % % % ... (2) num_subplot_cols in a given minigrid
% % % % ... (3) location_and_direction for a given minigrid
% % % % ... (4) index corresponding to the order of that minigrid from top to bottom of the full plot
% % % % % % % % % NOTE: for row 1 of the minigrid, the electrode numbers will correspond to ...
% % % % % % % % % % % % % % electrodes 1 through num_subplot_rows*num_subplot_rows for that minigrid.
% % % % % % % % % % % % % % Then, for row 2, the electrode numbers will correspond to ...
% % % % % % % % % % % % % % electrodes [ num_subplot_rows*num_subplot_rows for minigrid #1] + 1 through ...
% % % % % % % % % % % % % % electrodes [ num_subplot_rows*num_subplot_rows for minigrid #1] + [ num_subplot_rows*num_subplot_rows for minigrid #2]
% % % % % % % % % NOTE: the 4th column allows the minigrids to be reordered in the overall plot ROW-WISE
% % % % % % % % % % % % % % which means that this assumes minigrids are appended as rows (and therefore have equal numbers of columns)
% defaults (if nargin < 4 or if an input is empty) are: 16 / 16 / 'tlt' / { subplot_rows subplot_cols location_and_direction 1 }
% 
% OUTPUTS:
% all_subplot_nums ... 2D matrix showing electrode numbering visually
% each_subplot_num ... row vector with indices corresponding to 'subplot number' (3rd argument to Matlab's subplot function)
% % % basically, the electrode number at a given index of each_subplot_num will mean that
% % % % % that electrode number will be plotted in the indexed position on the grid plot
% % % e.g., if all_subplot_nums = [ 1 3 ; 2 4 ] , then each_subplot_num = [ 1 3 2 4 ]


if nargin < 1
    subplot_rows = 16 ;
    subplot_cols = 16 ;
    location_and_direction = 'tlt' ;
    grid_ordering = { subplot_rows subplot_cols location_and_direction 1 } ;
else % then subplot_rows is defined or empty
    if nargin < 2
        subplot_cols = 16 ;
        location_and_direction = 'tlt' ;
        grid_ordering = { subplot_rows subplot_cols location_and_direction 1 } ;
    else % then subplot_cols is defined or empty
        if nargin < 3
            location_and_direction = 'tlt' ;
            grid_ordering = { subplot_rows subplot_cols location_and_direction 1 } ;
        else % then location_and_direction is defined or empty
            if nargin < 4
                grid_ordering = { subplot_rows subplot_cols location_and_direction 1 } ;
%             else % then grid_ordering is defined or empty
            end
            
            if isempty(location_and_direction)
                location_and_direction = 'tlt' ;
            end
            
        end
        
        if isempty(subplot_cols)
            subplot_cols = 16 ;
        end
        
    end
    
    if isempty(subplot_rows)
        subplot_rows = 16 ;
    end
    
end

all_subplot_nums = NaN(subplot_rows,subplot_cols) ;

num_minigrids = size(grid_ordering,1) ;

these_subplot_rows = cell2mat(grid_ordering(:,1)) ; % column vector
% these_subplot_cols = cell2mat(grid_ordering(:,2)) ; % column vector
these_grid_orderings = cell2mat(grid_ordering(:,4)) ; % column vector

if length(these_subplot_rows) < num_minigrids
    these_subplot_rows = repelem(subplot_rows/num_minigrids,num_minigrids)' ;
end
if length(these_grid_orderings) < num_minigrids
    these_grid_orderings = 1:num_minigrids ;
    these_grid_orderings = these_grid_orderings' ;
end

these_subplot_rows = [ 0 these_subplot_rows(these_grid_orderings)' ] ;
% these_subplot_cols = [ 0 these_subplot_cols(these_grid_orderings)' ] ;

% for now, assume that minigrids will always be combined row-wise (i.e., have equal number of columns)
these_subplot_rows = cumsum(these_subplot_rows) ;
% these_subplot_cols = cumsum(these_subplot_cols) ;

these_subplot_init_nums = these_subplot_rows*subplot_cols ;

% keyboard

for this_minigrid = 1:num_minigrids
    
    this_grid_ordering = grid_ordering{this_minigrid,4} ;
    
    this_subplot_rows = grid_ordering{this_minigrid,1} ;
    if isempty(this_subplot_rows)
        this_subplot_rows = subplot_rows/num_minigrids ; % if number of rows for a minigrid is not listed, assume it is half of total
    end
    
    this_subplot_cols = grid_ordering{this_minigrid,2} ;
    if isempty(this_subplot_cols)
        this_subplot_cols = subplot_cols ; % if number of cols for a minigrid is not listed, assume it is equal to total
    end
    
    this_location_and_direction = grid_ordering{this_minigrid,3} ;
    if isempty(this_location_and_direction)
        this_location_and_direction = location_and_direction ; % if location_and_direction string for a minigrid is not listed, assume it is overall one
    end
    
    this_tlt_minigrid_nums = reshape(1:(this_subplot_rows*this_subplot_cols),this_subplot_cols,this_subplot_rows)' ; % 'tlt'
    this_tlt_minigrid_nums = orient_grid(this_tlt_minigrid_nums,this_location_and_direction) ;
    
%     keyboard
    
    this_tlt_minigrid_nums = this_tlt_minigrid_nums + these_subplot_init_nums(this_minigrid) ;
    
    % for now, assume that minigrids will always be combined row-wise (i.e., have equal number of columns)
    all_subplot_nums((these_subplot_rows(this_grid_ordering)+1):these_subplot_rows(this_grid_ordering+1),:) = this_tlt_minigrid_nums ;
    
end

each_subplot_num = reshape(all_subplot_nums',1,[]) ;

end


%% Orient Grid Function:
function all_subplot_nums = orient_grid(all_subplot_nums,location_and_direction)
% all_subplot_nums ......... 2D matrix with electrode numbers (e.g., 1-256) oriented as if the grid had location_and_direction of 'tlt'
% location_and_direction ... three-character string:
% % % % ... (1) t/b = location of first elec is Top/Bottom of grid
% % % % ... (2) l/r = location of first elec is Left/Right of grid
% % % % ... (3) t/b = location of last elec in the same strip as the first elec (e.g., #16 in a 16x16 grid) is Top/Bottom of grid
% % % % % % % % % % % % ... [indicates direction of the flow of numbering from first elec]
% 
% space of possible location_and_direction strings:
% location_and_direction can be one of the following: { 'tlt' 'tlb' 'trt' 'trb' 'blb' 'blt' 'brb' 'brt' } (default is 'tlt')
% 
% NOTE: IF all_subplot_nums IS NOT A SQUARE MATRIX
% % % % AND ~strcmp(location_and_direction(1),location_and_direction(3)) 
% % % % THEN taking transpose will mess this up!
% % % % thus, for non-square grids, strips of elecs (e.g., 1-16) should run row-wise
% 
% NOTE: order of flipud/fliplr does not matter BUT transpose must precede the flipping if both happen

if ~strcmp(location_and_direction(1),location_and_direction(3))
    if size(all_subplot_nums,1) ~= size(all_subplot_nums,2)
        error( 'for non-square grids, strips of elecs (e.g., 1-16) should run row-wise' )
    end
    all_subplot_nums = all_subplot_nums' ;
end

if strcmp(location_and_direction(2),'r')
    all_subplot_nums = fliplr(all_subplot_nums) ;
end

if strcmp(location_and_direction(1),'b')
    all_subplot_nums = flipud(all_subplot_nums) ;
end

end
