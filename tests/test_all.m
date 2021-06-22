%% This script tests all functions
% The function calls are executed in the necessary order and results
% compared against precomputed, saved results

%% Change to parent folder
cd ..

%% Check OS
if ~isunix
    slash = "\";
else
    slash = "/";
end

%% Specifiy all parameters and call "execute_all.m"

% Model
filename_mesh = "/home/paul/input/meshes/four_layer_surf_from_tets_11.mat";
%filename_mesh = "tests" + slash + "four_layer_surf_from_tets_3.mat";
filename_electrodes = "tests" + slash + "electrodes_five.txt";
filename_tissue = "tests" + slash + "fls_tissue.tiss";
filename_cond = "tests" + slash + "fls_conductivities.cond";

% Model preprocessing
numThreads = 4;     % Might not be valid on your machine!
% Current suggestions
TnumberE = 10;
GnumberE = 0;
RnumberP = 4;

% Iterative solver for charge distribution
% Current suggestions
iter   = 25;    % Number of iterations for iterative solver
relres = 1e-12; % Minimum acceptable relative residual for iterative solver
prec_charge   = 1e-2;  % Precision for FMM
weight = 1/2;   % Weight of the charge conservation law to be added

% Specific dipoles and their potential at electrodes
filename_dipoles = "tests" + slash + "dipoles_ten.txt";
% Current suggestions
R    = 0;       % Controls distance in which more exact integration of potential is used
prec_potential = 1e-2;    % Precision of FMM

% Results
filename_results = "tests" + slash + "result.mat";

% Start computation
execute_all(filename_mesh, filename_electrodes, filename_tissue, filename_cond, numThreads, TnumberE, GnumberE, RnumberP, iter, relres, prec_charge, weight, filename_dipoles, R, prec_potential, filename_results);

load(filename_results, 'cond', 'strge', 'dipole_ctr', 'dipole_moment', 'dipole_n', 'VoltageDifference', 'time_setup_electrodes', 'time_preprocess_model', 'time_solve_forward_problem_total', 'time_charge_engine', 'time_compute_dipole_potential');
VoltageDifferenceNum = horzcat(VoltageDifference{:});
dipole_ctr = dipole_ctr{2};
dipole_moment = dipole_moment{2};
disp(['Setup electrodes time: ' num2str(time_setup_electrodes)]);
disp(['Preprocess model time: ' num2str(time_preprocess_model)]);
disp(['Solve forward problem total time: ' num2str(time_solve_forward_problem_total)]);
disp(['Individual per electrode charge engine times: ']);
disp(vertcat(time_charge_engine{2:end}));
disp(['Individual per electrode compute dipole potential times: ']);
disp(vertcat(time_compute_dipole_potential{2:end}));

%% Compare with analytical solution

disp("Compute analytical solution and compare");
disp("The errors are computed per dipole over all of the electrodes");

dipole_length      = 0.0001;
dipole_magnitude   = sqrt(dot(dipole_moment, dipole_moment, 2));
dipole_orientation = dipole_moment./dipole_magnitude;
dipole_plus        = dipole_ctr + dipole_length/2.*dipole_orientation;
dipole_minus       = dipole_ctr - dipole_length/2.*dipole_orientation;
I0                 = dipole_magnitude./dipole_length;
% This is specific to the model used for this test
% Radius of the 4 tissue layers relative to a sphere with 100 mm radius
radfactor          = [0.92 0.86 0.80 0.78];

PotAnl = zeros(size(dipole_ctr, 1), strge.NumberOfElectrodes);

% Add paths
if ~isunix
    addpath(strcat(pwd, '\Engine'));
else
    addpath(strcat(pwd, '/Engine'));
end

No         = 1000;   % Length of computed series expansion of analytical solution
numThreads = 4;     % Might not be valid on your machine!
%parpool(numThreads);
for i = 1:size(dipole_ctr, 1)
    [PotAnl(i, :), ~] = a_p_4layer_infinite(No, I0(i, :), dipole_plus(i, :), dipole_minus(i, :), radfactor, cond, strge.PositionOfElectrodes);
    PotAnl(i, :)      = real(PotAnl(i, :));
end
%delete(gcp('nocreate'));

% Remove added paths
warning off; rmpath(genpath(pwd)); warning on;

VoltageDifferenceAnl = repmat(PotAnl(:, 1), 1, strge.NumberOfElectrodes) - PotAnl;
difference           = VoltageDifferenceAnl - VoltageDifferenceNum;
normalization_factor = sqrt(dot(VoltageDifferenceAnl, VoltageDifferenceAnl, 2));
Error2Norm           = sqrt(dot(difference, difference, 2))./normalization_factor
ErrorRDM             = VoltageDifferenceAnl./normalization_factor;
normalization_factor = sqrt(dot(VoltageDifferenceNum, VoltageDifferenceNum, 2));
ErrorRDM             = ErrorRDM - VoltageDifferenceNum./normalization_factor;
ErrorRDM             = 0.5*sqrt(dot(ErrorRDM, ErrorRDM, 2))

%% Change back to tests
cd tests
        
