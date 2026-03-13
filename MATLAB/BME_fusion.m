% ==============================================================================
% TITLE: Bayesian Maximum Entropy (BME) spatiotemporal data fusion
% DESCRIPTION: Integrates hard data, soft data and global offset from models 
%              to generate global ozone fields.
% ==============================================================================

clear;
close all;

% ==============================================================================
% 1. INITIALIZATION & WORKSPACE CLEANUP
% ==============================================================================
fprintf('>>> Initializing BME Pipeline and cleaning temporary workspace...\n');

if isfile('Covariance/Covariance.mat'), delete('Covariance/Covariance.mat'); end
if ~isempty(dir('BMEEstimation/*.mat')), delete('BMEEstimation/*.mat'); end
if ~isempty(dir('input data/*.mat')), delete('input data/*.mat'); end
if ~isempty(dir('*.mat')), delete('*.mat'); end

fprintf('    Temporary .mat workspace files cleared.\n');

% ==============================================================================
% 2. LOAD OBSERVATIONS AND SOFT DATA
% ==============================================================================
fprintf('\n>>> Loading Hard and Soft Datasets...\n');

% Define temporal scope
target_years = 2005:2022;

% Load hard data (ground/ozonesonde observations)
obs = getObs(min(target_years), max(target_years));

% Load soft data (Satellite retrievals with uncertainties)
csv_path = 'input data/ramp_satellite.csv'; 
if ~isfile(csv_path)
    error('Soft data CSV file not found: %s', csv_path);
end

fprintf('Reading soft data table: %s\n', csv_path);
soft_data_table = readtable(csv_path); % readtable automatically handles headers

% Define spatial estimation grid
grid = getEstimationGrid();

% ==============================================================================
% 3. EXTRACT MODEL COMPOSITES & CALCULATE RESIDUALS
% ==============================================================================
fprintf('\n>>> Processing Background Models and Residuals...\n');

% Load M3Fusion composite models into observation and grid formats
loadModelComposite;
[model_obs, model_grid, model_soft] = getModelValues(obs, grid, soft_data_table);

% Calculate residuals (Observation minus Model Background)
residuals.value = obs.Z - model_obs.value;
residuals.cMS   = obs.cMS;
residuals.tME   = obs.tME;

% ==============================================================================
% 4. SPATIOTEMPORAL COVARIANCE OPTIMIZATION
% ==============================================================================
fprintf('\n>>> Optimizing Spatiotemporal Covariance Parameters (Bayesian)...\n');

spaceUnit = 'degree';
varNameAndUnit = 'ppb'; 
doPlot = 1; % Set to 1 to enable interactive variogram fitting plots

% Perform Bayesian optimization for covariance modeling
% cov = getCovariance(residuals, doPlot, varNameAndUnit, spaceUnit);
cov = optimizeCovariance_Bayes(residuals);

fprintf('    Covariance optimization complete.\n');

% ==============================================================================
% 5. CONSTRUCT KNOWLEDGE BASE (HARD + SOFT DATA)
% ==============================================================================
fprintf('\n>>> Constructing BME Knowledge Base (Batch Processing)...\n');

output_folder = 'Output'; 
if ~exist(output_folder, 'dir'), mkdir(output_folder); end

nsmax = 20;        % Maximum number of soft data points to process locally
BMEprobaType = 1; % 1 = Gaussian distribution (L1 as mean, L2 as variance)

for i = 1:length(target_years)
    yr = target_years(i);
    fprintf('  -> Processing Knowledge Base for Year: %d ... ', yr);
    
    try
        % Generate Knowledge Base matrices (KG, KS)
        [KG, KS, BMEparam] = get_knowledgebase(obs, model_obs, model_soft, cov, ...
                                               soft_data_table, yr, ...
                                               nsmax, BMEprobaType);
        
        output_filename = sprintf('BME_Input_%d.mat', yr);
        save_path = fullfile(output_folder, output_filename);
        
        save(save_path, 'KG', 'KS', 'BMEparam');
        fprintf('SUCCESS | Saved to %s (Soft points: %d)\n', save_path, size(KS.ps, 1));
        
    catch ME
        fprintf('FAILED | Error: %s\n', ME.message);
    end
end

% ==============================================================================
% 6. BME SPATIOTEMPORAL ESTIMATION (SLIDING WINDOW)
% ==============================================================================
fprintf('\n>>> Executing BME Spatiotemporal Estimation (Sliding Window)...\n');

range_start_limit = min(target_years); 
range_end_limit   = max(target_years); 
window_size       = 2;

for i = 1:length(target_years)
    yr = target_years(i);
    
    try
        % Execute sliding window BME estimation
        T = runBME_SlidingWindow(yr, window_size, output_folder, ...
                                 model_grid, grid, ...
                                 range_start_limit, range_end_limit);
        
        % Export results to CSV
        csv_filename = fullfile(output_folder, sprintf('Ozone_Estimate_%d.csv', yr));
        writetable(T, csv_filename);
        fprintf('  -> Saved estimate: %s\n', csv_filename);
        
    catch ME
        warning('BME estimation failed for year %d: %s', yr, ME.message);
    end
end

% ==============================================================================
% 7. POST-PROCESSING (NETCDF EXPORT)
% ==============================================================================

% Create NetCDF out of fine resolution output
createNetCDF;

fprintf('\n>>> All BME data fusion tasks completed successfully.\n');