%% Script for running a simple lfmm3d setup to check validity of results
% This may be used to check if you version of lfmm3d runs correctly

%% Change to parent folder
cd ..

%% Add paths
if ~isunix
    addpath(strcat(pwd, '\Engine'));
    slash = "\";
else
    addpath(strcat(pwd, '/Engine'));
    slash = "/";
end

%% Load data
filename_mesh = "tests" + slash + "four_layer_surf_from_tets_3.mat";
load(filename_mesh, 'P');

%% Run FMM
    
pg              = 2;                        % potential and field are evaluated
srcinfo.sources = P';                       % source/target points
srcinfo.charges = srcinfo.sources(1, :);    % charges
prec            = 1e-2;                     % FMM precision
U               = lfmm3d(prec, srcinfo, pg);% FMM

filename_output = "tests" + slash + "results_lfmm3d.mat";
save(filename_output, 'pg', 'srcinfo', 'prec', 'U');

%% Remove added paths
warning off; rmpath(genpath(pwd)); warning on;

%% Change back to tests
cd tests