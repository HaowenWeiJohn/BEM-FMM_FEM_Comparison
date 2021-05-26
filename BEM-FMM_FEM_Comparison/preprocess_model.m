function [P, t, normals, Area, Center, Indicator, tissue, cond, enclosingTissueIdx, condin, condout, contrast, eps0, mu0, tneighbor, EC, PC] = ...
            preprocess_model(filename_mesh, filename_cond, filename_tissue, filename_output, filename_outputP, numThreads)
%   Imitates commands executed in "Model/model01_main_script.m".
%
%   Please see "read_cond" and "read_tissue" for specifications of ".cond"
%   and ".tiss" files
%   Conductivity must be in Amm!
%
%   Triangular mesh must be given as ".mat" file containing:
%       P - Mx3 matrix of mesh vertices in 3D (in mm!)
%       t - Nx3 matrix. Each row: 3 row-indices of P determining a triangle
%       normals - Nx3 matrix of unit outer normals of triangles in t
%       Indicator - Nx1 matrix of tissue indicators per triangle
%   The Indices must be integers 0, 1, ..., N
%   s.t. tissue i encloses tissue i+1 and start
%
%   The output filenames must be ".mat" files
%   One is the combined mesh, the other additional precomputed results
%
%   "numThreads" is the number of cores to be used for MATLAB Parallel
%   Pools
%
%   Modifications by Paul Lunkenheimer
%
%%   Original Documentation:
%
%   This is a mesh processor script: it computes basis triangle parameters
%   and necessary potential integrals, and constructs a combined mesh of a
%   multi-object structure (for example, a head or a whole body)
%
%   Copyright SNM/WAW 2017-2020
%
%%   Modifications:
%
%   Execute as function not as script
%   Returns all computed data
%   No GUI output
%   Modified verbosity
%   Different timers
%   Meshes are already complete
%   'tissue_index.txt' files not used
%   Does not create/save variable 'name'
%   Creates and Saves 'eps0' and 'mu0'

    %% Add paths
    if ~isunix
        s = pwd;
        addpath(strcat(s, '\Engine'));
        addpath(strcat(s, '\io'));
    else
        s = pwd;
        addpath(strcat(s, '/Engine'));
        addpath(strcat(s, '/io'));
    end
    
    %%  Define EM constants
    eps0        = 8.85418782e-012;  %   Dielectric permittivity of vacuum(~air)
    mu0         = 1.25663706e-006;  %   Magnetic permeability of vacuum(~air)

    %% Load tissue names, conducitivies and mesh
    tissue              = read_tissue(filename_tissue);
    cond                = read_cond(filename_cond);
    load(filename_mesh);
    Indicator           = Indicator + 1;
    enclosingTissueIdx  = (0:max(Indicator)-1)';
    P                   = P*1e-3;     % only if the original data were in mm!
    cond                = cond*1e3;   % only if the original data were in Amm!

    %%  Fix triangle orientation (just in case, optional)
    t = meshreorient(P, t, normals);

    %%   Process other mesh data
    Center      = 1/3*(P(t(:, 1), :) + P(t(:, 2), :) + P(t(:, 3), :));
    Area        = meshareas(P, t);

    %%  Assign facet conductivity information
    condambient = 0.0; %   air
    [contrast, condin, condout] = assign_initial_conductivities(cond, condambient, Indicator, enclosingTissueIdx);

    %%  Check for and process triangles that have coincident centroids
    [P, t, normals, Center, Area, Indicator, condin, condout, contrast] = ...
        clean_coincident_facets(P, t, normals, Center, Area, Indicator, condin, condout, contrast);
    N           = size(t, 1);

    %%   Find topological neighbors
    DT = triangulation(t, P); 
    tneighbor = neighbors(DT);
    % Fix cases where not all triangles have three neighbors
    tneighbor = pad_neighbor_triangles(tneighbor);

    %%   Save base data
    save(filename_output, 'P', 't', 'normals', 'Area', 'Center', 'Indicator', 'tissue', 'cond', 'enclosingTissueIdx', 'condin', 'condout', 'contrast', 'eps0', 'mu0');

    %%   Add accurate integration for electric field/electric potential on neighbor facets
    %   Indices into neighbor triangles
    RnumberE        = 128;      %   number of neighbor triangles for analytical integration of electric field
    RnumberP        = 128;      %   number of neighbor triangles for analytical integration of electric potential
    ineighborE      = knnsearch(Center, Center, 'k', RnumberE);   % [1:N, 1:RnumberE]
    ineighborP      = knnsearch(Center, Center, 'k', RnumberP);   % [1:N, 1:RnumberP]
    ineighborE      = ineighborE';          %   do transpose  
    ineighborP      = ineighborP';          %   do transpose  

    parpool(numThreads);
    %[EC, PC] = meshneighborints(P, t, normals, Area, Center, RnumberE, RnumberP, ineighborE, ineighborP, numThreads);
    EC = meshneighborints_En(P, t, normals, Area, Center, RnumberE, ineighborE);
    PC = meshneighborints_P(P, t, normals, Area, Center, RnumberP, ineighborP);
    delete(gcp('nocreate'));
    
    %%   Normalize sparse matrix EC by variable contrast (for speed up)
    N   = size(Center, 1);
    ii  = ineighborE;
    jj  = repmat(1:N, RnumberE, 1); 
    CO  = sparse(ii, jj, contrast(ineighborE));
    EC  = CO.*EC;

    save(filename_outputP, 'tneighbor', 'EC', 'PC', '-v7.3');
    
    %% Remove added paths
    warning off; rmpath(genpath(s)); warning on;

end