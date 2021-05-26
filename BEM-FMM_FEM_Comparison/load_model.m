function [P, t, normals, Area, Center, Indicator, tissue, cond, enclosingTissueIdx, condin, condout, contrast, eps0, mu0, tneighbor, EC, PC] = ...
            load_model(filename_model, filename_modelP)
%   Imitates commands executed in "bem0_load_model.m".
%
%   This function might not even be needed. One can simply load the data in
%   the executing script.
%   
%   Modifications by Paul Lunkenheimer
%
%%   Original Documentation:
%
%   This script loads mesh data into the MATLAB workspace. The data include
%   surface meshes and the potential integrals. It also loads the previous
%   solution (if exists)
%
%   Copyright SNM/WAW 2018-2020
%
%%   Modifications:
%   Does not add path to subfolders 'Engine' and 'Model'
%   Variables 'eps0' and 'mu0' are loaded from saved data
%   Does not load saved solutions
%   No GUI output
%   Modified verbosity
%   Different timers

    %% Load data
    load(filename_model, 'P', 't', 'normals', 'Area', 'Center', 'Indicator', 'tissue', 'cond', 'enclosingTissueIdx', 'condin', 'condout', 'contrast', 'eps0', 'mu0');
    load(filename_modelP, 'tneighbor', 'EC', 'PC');
        
end