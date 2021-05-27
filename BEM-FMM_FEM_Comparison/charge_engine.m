function [Einc, Pinc, c, resvec, conservation_law_error] = ...
    charge_engine(filename_model, iter, relres, weight, R, strdipolePplus, strdipolePminus, strdipolesig, strdipoleCurrent, Ctr, filename_solution)
%   Imitates commands executed in "bem2_charge_engine"
%
%   "filename_model" is the mesh data saved in "preprocess_model" under the
%   filename "filename_output" (not under "filename_outputP"!)
%
%   "iter" is maximum possible number of iterations in the solution
%   (is 25 in original script)
%
%   "relres" is minimum acceptable relative residual
%   (is 1e-12 in original script)
%
%   "weight" is weight of the charge conservation law
%   to be added (empirically found)
%   (is 1/2 in original script)
%
%   "R" is "radius of the enclosing sphere in m"??????????????????????
%   (is 0.01 or 0.005 in original script)
%
%   Variables "strdipolePplus", "strdipolePminus", "strdipolesig",
%   "strdipoleCurrent", "Ctr" come from function "setup_dipoles.m"
%
%   Set file to save to as "filename_solution"
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

    %% Load model
    load(filename_model, 'P', 't', 'normals', 'Area', 'Center', 'Indicator', 'tissue', 'cond', 'enclosingTissueIdx', 'condin', 'condout', 'contrast', 'eps0', 'mu0');

    %%  Right-hand side b of the matrix equation Zc = b
    %   Surface charge density is normalized by eps0: real charge density is eps0*c
    %   Gaussian integration is used here
    %[Einc, Pinc] = bemf3_inc_field_electric_gauss(strdipolePplus, strdipolePminus, strdipolesig, strdipoleCurrent, P, t);

    gaussRadius = 6 * R;
    [Einc, Pinc] = bemf3_inc_field_electric_gauss_selective(strdipolePplus, strdipolePminus, strdipolesig, strdipoleCurrent, P, t, Center, Ctr, gaussRadius);
    disp([newline 'Incident field calculated in ' num2str(toc) ' s']);

    b        = 2*(contrast.*sum(normals.*Einc, 2));                         %  Right-hand side of the matrix equation

    %%  GMRES iterative solution (native MATLAB GMRES is used)
    %   MATVEC is the user-defined function of c equal to the left-hand side of the matrix equation LHS(c) = b
    MATVEC = @(c) bemf4_surface_field_lhs(c, Center, Area, contrast, normals, weight, EC);     
    [c, flag, rres, its, resvec] = gmres(MATVEC, b, [], relres, iter, [], [], b);
    
    %%  Check charge conservation law (optional)
    conservation_law_error = sum(c.*Area)/sum(abs(c).*Area)

    %%   Topological low-pass solution filtering (repeat if necessary)
    % c = (c.*Area + sum(c(tneighbor).*Area(tneighbor), 2))./(Area + sum(Area(tneighbor), 2));

    %%  Save solution data (surface charge density, principal value of surface field)
    save(filename_solution, 'c', 'resvec', 'conservation_law_error');

%%   Not clear yet!! Do I need to compute and save this?? What about Einc, Pinc??
%%   Find and save surface electric potential
Padd = bemf4_surface_field_potential_accurate(c, Center, Area, PC);
%Padd = bemf4_surface_field_potential_subdiv(c, P, t, Area, 'barycentric', 3);
Ptot = Pinc + Padd;     %   Continuous total electric potential at interfaces

    %% Remove added paths
    warning off; rmpath(genpath(pwd)); warning on;

end