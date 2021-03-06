function [dipole_ctr, dipole_moment, dipole_n, VoltageDifference] = ...
                compute_dipole_potential(filename_dipoles, c, P, t, Center, Area, normals, electrodeCurrents, R, prec)
%   Imitates commands executed in "bem1_setup_dipole.m",
%   "bem1_setup_dipoles.m" and "bem1_precompute_dipoles.m" as well as
%   "bem6_comparison_analytical_all_dipoles_hor.m", 
%   "bem6_comparison_analytical_all_dipoles_ver.m" and
%   "bem6_comparison_analytical_single_dipole.m"
%
%   Attention: Current implementation only works for exactly 2 electrodes!
%
%   Unifies above scripts into one function and reads dipoles from ascii file
%   Dipole locations must be in mm and moments in Amm!
%   See "read_dipoles.m" for specifics
%
%   Does not compare with analytical solution!
%
%   "c" is solution charge distribution computed in "charge_engine.m"
%   "electrodeCurrents" as well
%
%   "P, t, Center, Area, normals" are usual model data computed in
%   "preprocess_model.m"
%
%   "R" controls distance in which more exact integration of potential is
%   used instead of FMM (suggestion for now: R=0)
%
%   "prec" is precision of FMM (suggestion for now: prec=1e-2)
%
%   "dipole_ctr", "dipole_moment", "dipole_n" are dipole data that has been
%   read from file
%
%   "VoltageDifference" is the voltage difference between the second and
%   the first electrode that were passed
%
%   Modifications by Paul Lunkenheimer
%
%%   Original Documentation "bem1_setup_dipole.m":
%
%   This script creates the base dipole (single)
%
%   Copyright SNM/WAW 2018-2020
%
%%   Original Documentation "bem1_setup_dipoles.m":
%
%   This script creates the cortical dipole layer (multiple dipoles)
%
%   Copyright SNM/WAW 2018-2020
%
%%   Original Documentation "bem1_precompute_dipoles.m":
%
%   This script compares analytical and numerical solutions for the
%   potential.
%   An infinitesimally short dipole within a four-layer sphere is considered
%
%   Copyright SNM 2018-2020
%
%%   Original Documentation "bem6_comparison_analytical_all_dipoles_hor.m":
%
%   This script compares analytical and numerical solutions for the
%   potential.
%   An infinitesimally short dipole within a four-layer sphere is considered
%
%   Copyright SNM 2018-2020
%
%%   Original Documentation "bem6_comparison_analytical_all_dipoles_ver.m":
%
%   This script compares analytical and numerical solutions for the
%   potential.
%   An infinitesimally short dipole within a four-layer sphere is considered
%
%   Copyright SNM 2018-2021
%
%%   Original Documentation "bem6_comparison_analytical_single_dipole.m":
%
%   This script compares analytical and numerical solutions for the
%   potential.
%   An infinitesimally short dipole within a four-layer sphere is considered
%
%   Copyright SNM 2018-2020
%
%%   Modifications:
%
%   Only set up location and moment of dipoles!
%   Execute as function not as script
%   Returns all computed data
%   No GUI output
%   Modified verbosity
%   Different timers
%   Some Variables from Script are now Parameters

    %% Add paths
    if ~isunix
        addpath(strcat(pwd, '\io'));
        addpath(strcat(pwd, '\Engine'));
    else
        addpath(strcat(pwd, '/io'));
        addpath(strcat(pwd, '/Engine'));
    end
    
    %% Load dipole(s) from file
    [dipole_ctr, dipole_moment, dipole_n] = read_dipoles(filename_dipoles);
    dipole_ctr                            = dipole_ctr*1e-03;       % only if the original data were in mm!
    dipole_moment                         = dipole_moment*1e-03;    % only if the original data were in Amm!
    
    %%  Compute electric potential for dipoles using reciprocity
    E = bemf5_volume_field_electric(dipole_ctr, c, P, t, Center, Area, normals, R, prec);

    VoltageDifference = - dot(E, dipole_moment, 2)/abs(electrodeCurrents(1));  % Reciprocity

    %% Remove added paths
    warning off; rmpath(genpath(pwd)); warning on;

end