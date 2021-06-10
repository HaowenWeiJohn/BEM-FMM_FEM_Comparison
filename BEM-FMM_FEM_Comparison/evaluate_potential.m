function [Ptot, Pinc, Padd] = evaluate_potential(electrodes, filename_model, filename_modelP, R, strdipolePplus, strdipolePminus, strdipolesig, strdipoleCurrent, Ctr)
%   Computes electric potential for given dipole and surface charge
%   solution
%
%   "electrodes" is a mx3 matrix of rows of evaluation points
%
%   "filename_model" and "filename_modelP" is the mesh data
%   saved in "preprocess_model"
%
%   "R" is used for Gaussian integration
%   (suggested for now: R=1)
%
%   Variables "strdipolePplus", "strdipolePminus", "strdipolesig",
%   "strdipoleCurrent", "Ctr" come from function "setup_dipoles.m"

    %% Add paths
    if ~isunix
        addpath(strcat(pwd, '\Engine'));
    else
        addpath(strcat(pwd, '/Engine'));
    end
    
    %%  Load model
    load(filename_model, 'P', 't', 'normals', 'Area', 'Center', 'Indicator', 'tissue', 'cond', 'enclosingTissueIdx', 'condin', 'condout', 'contrast', 'eps0', 'mu0');

    %%  Right-hand side b of the matrix equation Zc = b
    %   Surface charge density is normalized by eps0: real charge density is eps0*c
    %   Gaussian integration is used here
    %[Einc, Pinc] = bemf3_inc_field_electric_gauss(strdipolePplus, strdipolePminus, strdipolesig, strdipoleCurrent, P, t);

    gaussRadius = 6 * R;
    [Einc, Pinc] = bemf3_inc_field_electric_gauss_selective(strdipolePplus, strdipolePminus, strdipolesig, strdipoleCurrent, P, t, Center, Ctr, gaussRadius);
    
    %%  Find surface electric potential
    Padd = bemf4_surface_field_potential_accurate(c, Center, Area, PC);
    %%Padd = bemf4_surface_field_potential_subdiv(c, P, t, Area, 'barycentric', 3);
    Ptot = Pinc + Padd;     %   Continuous total electric potential at interfaces
    
    %% Remove added paths
    warning off; rmpath(genpath(pwd)); warning on;

end