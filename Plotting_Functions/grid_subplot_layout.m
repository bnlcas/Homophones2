function grd = grid_subplot_layout(corner_l, direction, varargin)
%% This function creats a 16 x 16 matrix that will permute the order of
% the subplot tiling given by plotGridPosition
%
% (Note that 16 x 16 grid is hard coded in this iteration)
%
% Inputs:
% corner_l: string indicated which corner of the grid the initial (CH 1) is
% located. Possible inputs:
% 'tr' (top right), 'tl' (top left), 
% 'br' (bottom right), 'bl' (bottom left)
% 
% direction: string indicated the tiling direction of the grid
% (direction from the first electrode (CH1)
% Possible inputs: 'x' (tile toward the right or left of the plot)
%                   'y' (tiel toward the top or bottom of the plot)
%
% Variable Inputs:
% must be preceded by a string flag:
% (Format: grid_subplot_layout(... , 'flag_name', input, ...)
% 
% 'plot_test_grid' - boolean, will plot a demo of the subplot tiling provided by grd
% (ex: , (..., 'plot_test_grid', true, ...)
%
% 'split_grid_flip' - string, will split the grid down a central axis 
% and flip the top/bottom (if 'vertical') or the left and right if 'horizontal.
% (ex: , (..., 'split_grid_flip', 'vertical',...) - flips the bottom half of electrodes with 
% and the top half of electrodes.
%
% 'rotate_split_grid' - string  ('left', 'right', 'top', 'bottom')
% will rotate half the grid by 180 degrees. This operation will be
% performed on the final configuration of the grid.
% left will rotate the left half of the grid and keep the bottom the same
% top will rotate the top half and keep the bottom half the same
% (ex: , (...,'rotate_split_grid', 'left',...)
%
%
% Output:
% grd - n_cols x n_rows matrix that delineates the tiling order for
% plotting, using the subplot tiling function (plotGridPosition.m)


%% Get Variable Inputs:
plot_test_grid = false;
plot_flag_ind = find(strcmpi(varargin, 'plot_test_grid'));
if ~isempty(plot_flag_ind)
    plot_test_grid = varargin{plot_flag_ind + 1};
end
split_grid_flip = false;
split_flag_ind = find(strcmpi(varargin, 'split_grid_flip'));
if ~isempty(split_flag_ind)
    split_grid_flip = true;
    grid_flip_dir = varargin{split_flag_ind+1};
    flip_vert = true;
    if strcmpi(grid_flip_dir, 'horizontal')
        flip_vert = false;
    elseif strcmpi(grid_flip_dir, 'vertical')
        flip_vert = true;
    else
        warning 'Invalid Axis to Flip Grid'
    end
end

rotate_split_grid = false;
split_rotate_flag_ind = find(strcmpi(varargin, 'rotate_split_grid'));
if ~isempty(split_rotate_flag_ind)
    rotate_split_grid = true;
    rotate_grid_half = varargin{split_rotate_flag_ind(1) + 1};
end

%% Set up base tiling:
grd_default = reshape(1:256,16,16);

if strcmpi(corner_l, 'tl')
    if strcmpi(direction, 'x')
        grd = grd_default;
    elseif strcmpi(direction,'y')
        grd = transpose(grd_default);
    else
        warning 'Invalid tiling direction'
    end
elseif strcmpi(corner_l, 'tr')
    if strcmpi(direction, 'x')
        grd = transpose(rot90(grd_default,3));
    elseif strcmpi(direction, 'y')
        grd = rot90(grd_default,3);
    else
        warning 'Invalid tiling direction'
    end
elseif strcmpi(corner_l, 'bl')
    if strcmpi(direction, 'x')
        grd = transpose(rot90(grd_default,1));
    elseif strcmpi(direction, 'y')
        grd = rot90(grd_default,1);
    else
        warning 'Invalid tiling direction'
    end
elseif strcmpi(corner_l, 'br')
    if strcmpi(direction, 'x')
        grd = rot90(grd_default,2);
    elseif strcmpi(direction, 'y')
        grd = transpose(rot90(grd_default,2));
    else
        warning 'Invalid tiling direction'
    end
else
    warning 'Invalid Base Corner Flag'
end

%% Flip grid on axis:
if split_grid_flip
    if flip_vert
        tmp = grd;
        grd(:,1:(size(grd,2)/2)) = tmp(:,(1+(size(grd,2)/2)):end);
        grd(:,(1+size(grd,2)/2):end) = tmp(:,1:(size(grd,2)/2));
    else
        tmp = grd;
        grd(1:(size(grd,2)/2),:) = tmp((1+size(grd,2)/2):end,:);
        grd((1+size(grd,2)/2):end,:) = tmp(1:(size(grd,2)/2),:);
    end
end

%% Flip Half of the grid:
if rotate_split_grid
    if strcmpi(rotate_grid_half,'right') | strcmpi(rotate_grid_half,'left')
        if ~strcmpi(direction, 'y')
            warning 'tiling and split directions are incompatible'
        end
        if strcmpi(rotate_grid_half,'right')
            grd(:,(1+size(grd,2)/2):end) = rot90(grd(:,(1+size(grd,2)/2):end),2);
        elseif strcmpi(rotate_grid_half, 'left')
            grd(:,1:(size(grd,2)/2)) = rot90(grd(:,1:(size(grd,2)/2)),2);
        end
    elseif strcmpi(rotate_grid_half, 'top') | strcmpi(rotate_grid_half, 'bottom')
        if ~strcmpi(direction, 'x')
            warning 'tiling and split directions are incompatible'
        end
        if strcmpi(rotate_grid_half, 'top')
            grd(:,(1+size(grd,2)/2):end) = rot90(grd(:,(1+size(grd,2)/2):end),2);
        elseif strcmpi(rotate_grid_half, 'bottom')
            grd(:,1:(size(grd,2)/2)) = rot90(grd(:,1:(size(grd,2)/2)),2);
        end
    else
        warning 'invalid half-section input ("right", "left", "top", "bottom")'
    end
end

if plot_test_grid
    figure;
    for i = 1:3:256
            p = plotGridPosition(grd(i)); 
            subplot('position',p);   
            text(0,(max(get(gca,'YLim')) - (max(get(gca,'YLim')) * 0.1)),num2str(i));
            ax = gca;
            ax.XTick = []; ax.YTick = [];
    end
end