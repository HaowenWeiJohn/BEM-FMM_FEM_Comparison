function [Ctr, NoDipoles, I0, strdipolePplus, strdipolePminus, dlength, strdipolesig, strdipoleCurrent, strdipolemvector, strdipolemcenter, strdipolemstrength] = ...
                setup_dipoles(filename_dipoles, dipole_length, dipole_tissue_index, cond, D)
%   Imitates commands executed in "bem1_setup_dipole.m",
%   "bem1_setup_dipoles.m" and "bem1_precompute_dipoles.m"
%
%   Unifies both scripts into one function and reads dipoles from file
%   Dipole locations must be in mm and moments in Amm!
%
%   Parameter "dipole_length" is the requested length of the dipoles loaded
%   from "filename_dipoles"
%   It must be given in m and may be a scalar as well as a row-vector
%   containing one length per dipole
%
%   Parameter "dipole_tissue_index" gives the tissue index the dipole is
%   located in
%   Please "preprocess_model" or "load_model" to check the index of
%   a tissue in the variable "tissue"
%   "dipole_tissue_index" can be one scalar or a row-vector
%   containing one tissue index per dipole
%
%   "D" is number of smaller subdipoles for magnetic dipole
%   subdivision (is 1 or 10 in original script)
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
%%   Modifications:
%
%   Execute as function not as script
%   Returns all computed data
%   No GUI output
%   Modified verbosity
%   Different timers
%   Some Variables from Script are now Parameters

    %% Add paths
    if ~isunix
        addpath(strcat(pwd, '\io'));
    else
        addpath(strcat(pwd, '/io'));
    end

    %% Load dipole(s) from file
    [Ctr, dipole_moment, NoDipoles]    = read_dipoles(filename_dipoles);
    Ctr = Ctr*1e-03;                                   % only if the original data were in mm!
    dipole_magnitude        = sqrt(sum(dipole_moment.*dipole_moment, 2));
    dipole_orientation      = dipole_moment./dipole_magnitude;
    dipole_magnitude        = dipole_magnitude*1e-03;  % only if the original data were in Amm!
    I0                      = dipole_magnitude./dipole_length;
    dlength                 = dipole_length;
    dipole_length           = dipole_length./2;
    strdipolePplus          = Ctr + dipole_length.*dipole_orientation;
    strdipolePminus         = Ctr - dipole_length.*dipole_orientation;
    if (length(dipole_tissue_index)==1)
        strdipolesig            = repmat(cond(dipole_tissue_index), 1, 2);
    elseif (all(length(dipole_tissue_index)==[NoDipoles 1]))
        strdipolesig            = repmat(cond(dipole_tissue_index), NoDipoles, 2);
    else
        error('Variable "dipole_tissue_index" must be scalar or a row-vector containing one tissue index per dipole');
    end
    strdipoleCurrent       = [+I0; -I0];
    
    %% Remove added paths
    warning off; rmpath(genpath(pwd)); warning on;

end