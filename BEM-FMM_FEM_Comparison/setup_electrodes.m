function [] = ...
                setup_electrodes(filename_mesh, filename_electrodes, filename_output_mesh, filename_output_electrodes)
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
%   No GUI output
%   Modified verbosity
%   Different timers
%   Some Variables from Script are now Parameters
%   Works with whole mesh not with one mesh per layer

    %% Add paths
    if ~isunix
        addpath(strcat(pwd, '\io'));
        addpath(strcat(pwd, '\Electrodes'));
    else
        addpath(strcat(pwd, '/io'));
        addpath(strcat(pwd, '/Electrodes'));
    end
    
    %% Load mesh to imprint electrodes
    load(filename_mesh, 'P', 't', 'normals', 'Indicator');
    t_outer = t(Indicator==0);
    normals_outer = normals(Indicator==0);
    
    %% Load electrodes
    
    
    %% Imprint electrodes
    [P_outer, t_outer, normals_outer, IndicatorElectrodes] = meshimprint(P, t_outer, normals_outer, strge);
    
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

    index               = IndicatorElectrodes>0;
    t_outer             = [t_outer(index, :); t_outer(~index, :)];
    normals_outer       = [normals_outer(index, :); normals_outer(~index, :)];
    IndicatorElectrodes = IndicatorElectrodes(index);
    IndicatorElectrodes(end+1:size(t_outer, 1)) = 0;
    
    %% Insert into whole mesh
    t_outer                 = size(P, 1) + t_outer;
    P                       = [P; P_outer];    
    t(Indicator==0)         = [];
    t                       = [t_outer; t];
    normals(Indicator==0)   = [];
    normals                 = [normals_outer; normals];
    Indicator(Indicator==0) = [];
    Indicator               = [zeros(size(t_outer, 1), 1); Indicator];
    [P, t]                  = fixmesh(P, t);
    
    %% Save
    save(filename_output_electrodes, 'IndicatorElectrodes', 'strge');
    save(filename_output_mesh, 'P', 't', 'normals');
    
    %% Remove added paths
    warning off; rmpath(genpath(pwd)); warning on;

end