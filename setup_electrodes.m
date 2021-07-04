function [P, t, normals, Indicator, IndicatorElectrodes, strge] = ...
             setup_electrodes(filename_mesh, filename_electrodes)
%   Imitates commands executed in "Electrodes/electrodes01_imprint.m"
%
%   Triangular mesh must be given as ".mat" file containing:
%       P - Mx3 matrix of mesh vertices in 3D (in mm!)
%       t - Nx3 matrix. Each row: 3 row-indices of P determining a triangle
%       normals - Nx3 matrix of unit outer normals of triangles in t
%       Indicator - Nx1 matrix of tissue indicators per triangle
%   The Indicators must be integers 0, 1, ..., N
%   s.t. tissue i encloses tissue i+1
%   Mesh must be in mm!
%
%   Electrodes must be given as ascii file.
%   Each electrode is seperated only by single line breaks
%   Each line contains 3 entries seperated by single spaces
%   Entries are x, y and z coordinate of electrode location
%
%   "IndicatorElectrodes" is Nx1 vector containing per triangle
%   the electrode number it belongs to or zero
%
%   Struct "strge" contains positions of electrodes read from ascii file,
%   the number of electrodes and their radius
%
%   This function does not change any units. Thus, if electrode position
%   and radius are given e.g. in mm, they will still be in mm when returned
%   in "strge"
%
%   Modifications by Paul Lunkenheimer
%
%%   Original Documentation:
%
%   This is an electrode processor script: it imprints an arbitrary number
%   of electrodes
%
%   Copyright SNM 2012-2020
%
%%   Modifications:
%
%   Execute as function not as script
%   Returns all computed data
%   Does not save any data
%   No GUI output
%   Modified verbosity
%   Different timers
%   Some Variables from Script are now Parameters
%   Works with whole mesh not with one mesh per layer

    %% Add paths
    if ~isunix
        addpath(strcat(pwd, '\io'));
        addpath(strcat(pwd, '\Electrodes'));
        addpath(strcat(pwd, '\Engine'));
    else
        addpath(strcat(pwd, '/io'));
        addpath(strcat(pwd, '/Electrodes'));
        addpath(strcat(pwd, '/Engine'));
    end
    
    %% Load mesh to imprint electrodes
    load(filename_mesh, 'P', 't', 'normals', 'Indicator');
    t_outer       = t(Indicator==0, :);
    normals_outer = normals(Indicator==0, :);
    P_outer = P;
    [P_outer, t_outer] = fixmesh(P_outer, t_outer);
    
    edges_outer               = meshconnee(t_outer);
    temp                      = P_outer(edges_outer(:, 1), :) - P_outer(edges_outer(:, 2), :);
    Mesh.AvgEdgeLength        = mean(sqrt(dot(temp, temp, 2)));
    
    %% Load electrodes    
    strge.PositionOfElectrodes = read_electrodes(filename_electrodes);
    strge.NumberOfElectrodes   = size(strge.PositionOfElectrodes, 1);
    RadE                       = 1.5*Mesh.AvgEdgeLength;    % electrode radius in mm (at least 3 triangles along the diameter) 
    strge.RadiusOfElectrodes   = RadE*ones(1, strge.NumberOfElectrodes);    % in mm here
    
    %% Prepare imprinting of electrodes - refine mesh at electrodes
    for m = 1:strge.NumberOfElectrodes    
        [P_outer, t_outer, normals_outer] = mesh_refinement(P_outer, t_outer, normals_outer, strge.PositionOfElectrodes(m, :), RadE);
    end
    
    %% Imprint electrodes
    [P_outer, t_outer, normals_outer, IndicatorElectrodes] = meshimprint(P_outer, t_outer, normals_outer, strge);
    
    %% Put electrodes up front sequentially (1, 2, 3)
    tt = [];
    nn = [];
    ie = [];
    for m = 1:strge.NumberOfElectrodes 
        index  = IndicatorElectrodes==m;
        tt     = [tt; t_outer(index, :)];
        nn     = [nn; normals_outer(index, :)];
        ie     = [ie; m*ones(sum(index), 1)];
    end
    index               = IndicatorElectrodes==0;
    t_outer             = [tt; t_outer(index, :)];
    normals_outer       = [nn; normals_outer(index, :)];
    IndicatorElectrodes = [ie; IndicatorElectrodes(index)];
    
    %% Insert into whole mesh
    P                       = [P_outer; P];   
    t(Indicator==0, :)      = [];
    t                       = t + size(P_outer, 1);
    t                       = [t_outer; t];
    normals(Indicator==0, :)= [];
    normals                 = [normals_outer; normals];
    Indicator(Indicator==0) = [];
    Indicator               = [zeros(size(t_outer, 1), 1); Indicator];
    [P, t]                  = fixmesh(P, t);
    
    %% Remove added paths
    warning off; rmpath(genpath(pwd)); warning on;

end