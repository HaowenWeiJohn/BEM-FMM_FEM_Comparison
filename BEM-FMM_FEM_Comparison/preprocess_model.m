function [P, t, normals, Area, Center, Indicator, tissue, cond, enclosingTissueIdx, condin, condout, contrast, eps0, mu0, tneighbor, EC, PC] = ...
            preprocess_model(filename_mesh, filename_cond, filename_tissue, filename_output, filename_outputP, numThreads, RnumberE, RnumberP)
%   Imitates commands executed in "Model/model01_main_script.m"
%
%   First run "setup_electrodes" as mesh must be refined around electrodes!

%   Please see "read_cond" and "read_tissue" for specifications of ".cond"
%   and ".tiss" files
%   Conductivity must be in S/mm!
%
%   Triangular mesh must be given as ".mat" file containing:
%       P - Mx3 matrix of mesh vertices in 3D (in mm!)
%       t - Nx3 matrix. Each row: 3 row-indices of P determining a triangle
%       normals - Nx3 matrix of unit outer normals of triangles in t
%       Indicator - Nx1 matrix of tissue indicators per triangle
%   The Indicators must be integers 0, 1, ..., N
%   s.t. tissue i encloses tissue i+1
%   Mesh must be in mm!
%   First run "setup_electrodes" as mesh must be refined around electrodes!
%
%   The output filenames must be ".mat" files
%   One is the combined mesh, the other additional precomputed results
%
%   "numThreads" is the number of cores to be used for MATLAB Parallel
%   Pools
%
%   "RnumberE" is number of neighbor triangles for analytical integration
%   of electric field (suggestion for now: RunmberE=4)
%
%   "RnumberP" is number of neighbor triangles for analytical integration
%   of electric potential (suggestion for now: RnumberP=4)
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
%   Some Variables from Script are now Parameters

    %% Add paths
    if ~isunix
        addpath(strcat(pwd, '\Engine'));
        addpath(strcat(pwd, '\io'));
        addpath(strcat(pwd, '\Electrodes'));
    else
        addpath(strcat(pwd, '/Engine'));
        addpath(strcat(pwd, '/io'));
        addpath(strcat(pwd, '\Electrodes'));
    end
    
    %%  Define EM constants
    eps0        = 8.85418782e-012;  %   Dielectric permittivity of vacuum(~air)
    mu0         = 1.25663706e-006;  %   Magnetic permeability of vacuum(~air)

    %% Load tissue names, conducitivies and mesh
    tissue              = read_tissue(filename_tissue);
    cond                = read_cond(filename_cond);
    load(filename_mesh, 'P', 't', 'normals', 'Indicator');
    assert(length(tissue)==length(unique(Indicator)));
    Indicator           = Indicator + 1;
    enclosingTissueIdx  = (0:max(Indicator)-1)';
    P                   = P*1e-3;     % only if the original data were in mm!
    cond                = cond*1e3;   % only if the original data were in Amm!
    
    %% Refine two inner shells
    [P3, t3, normals3]          = meshrefiner(P, t(Indicator==3, :), normals(Indicator==3, :));
    [P4, t4, normals4]          = meshrefiner(P, t(Indicator==4, :), normals(Indicator==4, :));
    t3                          = t3 + size(P, 1);
    t4                          = t4 + size(P, 1) + + size(P3, 1);
    P                           = [P; P3; P4];
    t(Indicator==3, :)          = [];
    t(Indicator==4, :)          = [];
    t                           = [t; t3; t4];
    normals(Indicator==3, :)    = [];
    normals(Indicator==4, :)    = [];
    normals                     = [normals; normals3; normals4];
    Indicator(Indicator==3)     = [];
    Indicator(Indicator==4)     = [];
    Indicator                   = [Indicator; repmat(3, size(t3, 1), 1); repmat(4, size(t4, 1), 1)];
    [P, t]                      = fixmesh(P, t);    
    
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

    %%   Find topological neighbors
    DT = triangulation(t, P); 
    tneighbor = neighbors(DT);
    % Fix cases where not all triangles have three neighbors
    tneighbor = pad_neighbor_triangles(tneighbor);

    %%   Save base data
    save(filename_output, 'P', 't', 'normals', 'Area', 'Center', 'Indicator', 'tissue', 'cond', 'enclosingTissueIdx', 'condin', 'condout', 'contrast', 'eps0', 'mu0');

    %%   Add accurate integration for electric field/electric potential on neighbor facets
    %   Indices into neighbor triangles
    ineighborE      = knnsearch(Center, Center, 'k', RnumberE);   % [1:N, 1:RnumberE]
    ineighborP      = knnsearch(Center, Center, 'k', RnumberP);   % [1:N, 1:RnumberP]
    ineighborE      = ineighborE';          %   do transpose  
    ineighborP      = ineighborP';          %   do transpose  

    parpool(numThreads);
    EC                  = meshneighborints_En(P, t, normals, Area, Center, RnumberE, ineighborE);
    [PC, integralpd]    = meshneighborints_P(P, t, normals, Area, Center, RnumberP, ineighborP);
    delete(gcp('nocreate'));
    
    %%   Normalize sparse matrix EC by variable contrast (for speed up)
    N   = size(Center, 1);
    ii  = ineighborE;
    jj  = repmat(1:N, RnumberE, 1); 
    CO  = sparse(ii, jj, contrast(ineighborE));
    EC  = CO.*EC;

    save(filename_outputP, 'tneighbor', 'EC', 'PC', '-v7.3');
    
    %%  Electrode preconditioner M (left). Electrodes may be assigned to different tissues
    tic
    load electrode_data;
    ElectrodeIndexes = cell(max(IndicatorElectrodes), 1);
    for j = 1:max(IndicatorElectrodes)
        ElectrodeIndexes{j} = find(IndicatorElectrodes==j);
    end
    indexe = transpose(vertcat(ElectrodeIndexes{:}));   %   this index is not contiguous

    Ne          = length(indexe); 
    tempC       = Center(indexe, :); 
    tempA       = Area(indexe); 
    A           = repmat(tempA, 1, length(tempA));
    M           = (1/(4*pi))*1./dist(tempC').*A';       %   base preconditioner matrix
    for m = 1:Ne                                        %   base matrix with zero elements
        M(m, m) = 0;
    end
    for m = 1:Ne                                                    %   put in neighbor integrals
        indexf                  = indexe(m);                        %   global facet number on electrodes
        mneighbors              = ineighborP(:, indexf);            %   global neighbor numbers of facet indexf
        temp                    = intersect(indexe, mneighbors);    %   global numbers for neighbors of facet indexf within electrodes
        for n = 1:length(temp)
            indexintoM  = find(temp(n)==indexe);                    %   local number for the neighbor temp(n) in indexe - into matrix
            indexintoPD = find(mneighbors==temp(n));                %   local index for the neighbor temp(n) in integralpd - into integralpd
            M(indexintoM, m)  = M(indexintoM, m) + integralpd(indexintoPD, indexf)/(4*pi);
        end
    end
    M = inv(M);                                        %   direct inversion - replace
    disp([newline 'Preconditioner computed in ' num2str(toc) ' s']);
    
    %% Remove added paths
    warning off; rmpath(genpath(pwd)); warning on;

end