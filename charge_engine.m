function [c, resvec, electrodeCurrents, En_loc] = ...
    charge_engine(normals, Area, Center, condin, contrast, EC, PC, M, ElectrodeIndexes_global, indexe, ElectrodeIndexes_local, V, iter, relres, prec, weight)
%   Imitates commands executed in "bem2_charge_engine"
%
%   "filename_model" is the mesh data saved in "preprocess_model" under the
%   filename "filename_output" (not under "filename_outputP"!)
%
%   "filename_modelP" is the mesh data saved in "preprocess_model" under the
%   filename "filename_outputP"
%
%   "iter" is maximum possible number of iterations in the solution
%   (is 25 in original script)
%
%   "relres" is minimum acceptable relative residual
%   (is 1e-12 in original script)
%
%   precis is the precision for the FMM computations
%   (is 1e-2 in original script)
%
%   "weight" is weight of the charge conservation law
%   to be added (empirically found)
%   (is 1/2 in original script)
%
%   Variables "strdipolePplus", "strdipolePminus", "strdipolesig",
%   "strdipoleCurrent", "Ctr" come from function "setup_dipoles.m"
%
%%   Original Documentation:
%
%   This script computes the induced surface charge density for an
%   inhomogeneous multi-tissue object given the primary electric field, with
%   accurate neighbor integration
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
%   Model data is loaded from file by function and not in workspace
%   Some Variables from Script are now Parameters

    %% Add paths
    if ~isunix
        addpath(strcat(pwd, '\Engine'));
    else
        addpath(strcat(pwd, '/Engine'));
    end

    %% Iterative solution of linear equation
    % Right-hand side b of the matrix equation Zc = b
    % Surface charge density is normalized by eps0: real charge density is eps0*c
    b           = zeros(size(normals, 1), 1);   % Right-hand side of the matrix equation
    b(indexe)                   = M*V;  % Electrodes held at constant voltage
    %  GMRES iterative solution     
    MATVEC                      = @(c) bemf4_surface_field_lhs_v(c, Center, Area, contrast, normals, M, EC, PC, indexe, weight, condin, prec);
    [c, flag, rres, its, resvec]= gmres(MATVEC, b, [], relres, iter, [], [], b);

    %% Compute total current trough electrodes from charge distribution solution c
    % Normal electric field just inside
    % Only computed at triangle centers that belong to electrodes
    En_loc = bemf4_surface_field_electric_accurate(c, Center, Area, normals, EC, prec, indexe);
    electrodeCurrents = zeros(length(ElectrodeIndexes_global), 1);
    for j = 1:length(ElectrodeIndexes_global)
        % Indices of triangles belonging to electrode j in whole mesh
        index_glob           = ElectrodeIndexes_global{j};
        % Indices of triangles belonging to electrode j in reduced mesh consisting of only electrode triangles
        index_loc            = ElectrodeIndexes_local{j};
        % Compute Currents
        electrodeCurrents(j) = -sum(En_loc(index_loc).*Area(index_glob).*condin(index_glob));
    end

    %% Surface electric potential everywhere
    %Ptot = bemf4_surface_field_potential_accurate(c, Center, Area, PC);

    %% Remove added paths
    warning off; rmpath(genpath(pwd)); warning on;

end