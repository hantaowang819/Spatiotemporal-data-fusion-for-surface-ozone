function [model_obs, model_grid, model_soft] = getModelValues(obs, grid, soft_data_table) 
% ==============================================================================
% TITLE: Extract Model Background Values (Global Offset)
% DESCRIPTION: Interpolates the multi-model composite (background) to the 
%              spatiotemporal coordinates of hard observations, soft data, 
%              and the final estimation grid. Applies distance masking to 
%              prevent excessive extrapolation.
%
% INPUTS:
%   obs             : Struct containing hard data spatial/temporal coordinates
%   grid            : N x 2 matrix of estimation grid coordinates [lon, lat]
%   soft_data_table : Table containing soft data coordinates (lon, lat)
%
% OUTPUTS:
%   model_obs       : Struct with background values at hard data locations
%   model_grid      : Struct with background values at estimation grid nodes
%   model_soft      : Struct with background values at soft data locations
% ==============================================================================

DIR_INPUT = 'input data';    
CUTOFF_DISTANCE = 0.375; % Maximum interpolation radius in degrees

file_obs  = 'ModelAtObs.mat';
file_grid = 'ModelAtGrid.mat';
file_soft = 'ModelAtSoft.mat'; 

% ------------------------------------------------------------------------------
% 1. CHECK FOR PRE-COMPUTED DATA
% ------------------------------------------------------------------------------
if isfile(file_obs) && isfile(file_grid) && isfile(file_soft)
    fprintf('>>> Loading pre-computed Global Offset (GO) matrices...\n');
    load(file_obs, 'model_obs');
    load(file_grid, 'model_grid');
    load(file_soft, 'model_soft');
    
else
    fprintf('>>> Computing Global Offsets via spatial interpolation...\n');
    
    % Initialize coordinate and time arrays to match respective datasets
    model_obs.cMS   = obs.cMS;
    model_obs.tME   = obs.tME;
    model_obs.value = [];

    model_grid.cMS   = grid;
    model_grid.tME   = obs.tME;
    model_grid.value = [];
   
    soft_coords      = [soft_data_table.lon, soft_data_table.lat];
    model_soft.cMS   = soft_coords;
    model_soft.tME   = obs.tME;
    model_soft.value = [];

    % --------------------------------------------------------------------------
    % 2. ITERATE OVER TIME SLICES (YEARS)
    % --------------------------------------------------------------------------
    for i = 1:length(obs.tME)
        
        year = obs.tME(i);
        fprintf('  -> Processing Year: %d\n', year);
        
        % Load raw model composite for the current year
        composite_path = fullfile(DIR_INPUT, sprintf('ModelComposite%d.mat', year));
        if ~isfile(composite_path)
            error('Model composite file missing: %s', composite_path);
        end
        load(composite_path, 'model'); 
        
        % Optimize: Construct Delaunay Triangulation once per year for distance masking
        DT = delaunayTriangulation(model.lon, model.lat);
        
        % --- A. Interpolate at Hard Observations ---
        gridval_obs = griddata(model.lon, model.lat, model.val, obs.cMS(:,1), obs.cMS(:,2));
        [~, d_obs]  = nearestNeighbor(DT, obs.cMS(:,1), obs.cMS(:,2));
        gridval_obs(d_obs > CUTOFF_DISTANCE) = NaN;
        
        model_obs.value = [model_obs.value, gridval_obs];
        
        % --- B. Interpolate at Estimation Grid ---
        gridval_mod = griddata(model.lon, model.lat, model.val, grid(:,1), grid(:,2));
        [~, d_grid] = nearestNeighbor(DT, grid(:,1), grid(:,2));
        gridval_mod(d_grid > CUTOFF_DISTANCE) = NaN;        
        
        model_grid.value = [model_grid.value, gridval_mod];
        
        % --- C. Interpolate at Soft Data Locations ---
        gridval_soft = griddata(model.lon, model.lat, model.val, soft_coords(:,1), soft_coords(:,2));
        [~, d_soft]  = nearestNeighbor(DT, soft_coords(:,1), soft_coords(:,2));
        gridval_soft(d_soft > CUTOFF_DISTANCE) = NaN;
        
        model_soft.value = [model_soft.value, gridval_soft];
       
    end
    
    % --------------------------------------------------------------------------
    % 3. EXPORT MATRICES TO DISK
    % --------------------------------------------------------------------------
    fprintf('  -> Saving computed GO matrices to disk...\n');
    save(file_obs, 'model_obs');   
    save(file_grid, 'model_grid'); 
    save(file_soft, 'model_soft'); 
    
end

end % End of getModelValues function