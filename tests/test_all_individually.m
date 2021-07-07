%% This script tests all functions
% The function calls are executed in the necessary order and the result is
% compared against the analytical solution as well as a precomputed
% numerical solution.

machine = "laptop";

%% Change to parent folder
cd ..

%% Check OS
if ~isunix
    slash = "\";
else
    slash = "/";
end

%% Set up electrodes and refine mesh using "setup_electrodes.m"
% Loads mesh and electrode positions.
% In a radius of 1.5*AverageEdgeLength all triangles around an electrode
% position are marked as belonging to this electrode patch.
% The mesh is refined in this region.

%filename_mesh = "C:\Users\Paul\Documents\WWU\Masterarbeit\Meshing\multilayer_sphere_model\4_Layer_Sphere_Meshes\Surface_Meshes\Surface_meshes_coupled_with_volume_meshes\MAT\four_layer_surf_from_tets_7.mat";
filename_mesh = "tests" + slash + "four_layer_surf_from_tets_3.mat";
filename_electrodes = "tests" + slash + "electrodes_two.txt";

disp("Imprint the electrodes on the mesh");
[P, t, normals, Indicator, IndicatorElectrodes, strge] = ...
    setup_electrodes(filename_mesh, filename_electrodes);

% Scale from mm to m
strge.PositionOfElectrodes = strge.PositionOfElectrodes*1e-3;
strge.RadiusOfElectrodes   = strge.RadiusOfElectrodes*1e-3;

save("find_Sky_error_save1" + machine + ".mat", 'P', 't', 'normals', 'Indicator', 'IndicatorElectrodes', 'strge', '-v7.3');

%% Preprocess model using "preprocess_model.m"
% The mesh, tissue names and tissue conductivities are loaded.
% Additional data like the area per triangle, the triangle
% centers and conducivity (in, out and contrast) per triangle is computed.
% Correction matrices EC and PC for electric field and potential are
% computed. These correct the FMM computations by more exact integration of
% electric activity induced by near neighbors.
% The matrix M is a preconditioner for the iterative solver used in the
% electrode regions.

filename_tissue = "tests" + slash + "fls_tissue.tiss";
filename_cond = "tests" + slash + "fls_conductivities.cond";

numThreads = 4;     % Might not be valid on your machine!
% Current suggestions
TnumberE = 10;
GnumberE = 0;
RnumberP = 4;

disp("Process the mesh and model data");
[P, t, normals, Area, Center, Indicator, tissue, enclosingTissueIdx, cond, condin, condout, contrast, eps0, mu0, EC, PC, M, integralpd, ineighborE, ineighborP, ElectrodeIndexes_global, ElectrodeIndexes_local, V] = ...
            preprocess_model(P, t, normals, Indicator, IndicatorElectrodes, filename_cond, filename_tissue, numThreads, TnumberE, GnumberE, RnumberP);
        
save("find_Sky_error_save2" + machine + ".mat", 'TnumberE', 'GnumberE', 'RnumberP', 'P', 't', 'normals', 'Area', 'Center', 'Indicator', 'tissue', 'enclosingTissueIdx', 'cond', 'condin', 'condout', 'contrast', 'eps0', 'mu0', 'EC', 'PC', 'M', 'integralpd', 'ineighborE', 'ineighborP', 'ElectrodeIndexes_global', 'ElectrodeIndexes_local', 'V', '-v7.3');        
        
%% Compute solution
% This solution is only linked to the electrodes, but not yet to specific
% dipoles.

% Current suggestions
iter   = 25;    % Number of iterations for iterative solver
relres = 1e-12; % Minimum acceptable relative residual for iterative solver
prec   = 1e-2;  % Precision for FMM
weight = 1/2;   % Weight of the charge conservation law to be added

disp("Compute iterative solution for charge distribution");
[c, resvec, electrodeCurrents, En] = ...
    charge_engine(normals, Area, Center, condin, contrast, EC, PC, M, ElectrodeIndexes_global, ElectrodeIndexes_local, V, iter, relres, prec, weight);

save("find_Sky_error_save3" + machine + ".mat", 'iter', 'relres', 'prec', 'weight', 'c', 'resvec', 'electrodeCurrents', 'En', '-v7.3');        

%% Compute potential at electrodes generated by dipoles using reciprocity

filename_dipoles = "tests" + slash + "dipoles_ten.txt";

% Current suggestions
R    = 0;       % Controls distance in which more exact integration of potential is used
prec = 1e-2;    % Precision of FMM

disp("Compute potential evoked by given dipoles");
[dipole_ctr, dipole_moment, dipole_n, VoltageDifferenceNum] = ...
    compute_dipole_potential(filename_dipoles, c, P, t, Center, Area, normals, electrodeCurrents, R, prec);

save("find_Sky_error_save4" + machine + ".mat", 'R', 'prec', 'dipole_ctr', 'dipole_moment', 'dipole_n', 'VoltageDifferenceNum', '-v7.3');
            
%% Compare with analytical solution

disp("Compute analytical solution and compare");
disp("The errors are computed over all of the dipoles");

dipole_length      = 0.0001;
dipole_magnitude   = sqrt(dot(dipole_moment, dipole_moment, 2));
dipole_orientation = dipole_moment./dipole_magnitude;
dipole_plus        = dipole_ctr + dipole_length/2.*dipole_orientation;
dipole_minus       = dipole_ctr - dipole_length/2.*dipole_orientation;
I0                 = dipole_magnitude./dipole_length;
% This is specific to the model used for this test
% Radius of the 4 tissue layers relative to a sphere with 100 mm radius
radfactor          = [0.92 0.86 0.80 0.78];

PotAnl = zeros(size(dipole_ctr, 1), 2);

% Add paths
if ~isunix
    addpath(strcat(pwd, '\Engine'));
else
    addpath(strcat(pwd, '/Engine'));
end

No         = 100;   % Length of computed series expansion of analytical solution
numThreads = 4;     % Might not be valid on your machine!
%parpool(numThreads);
for i = 1:size(dipole_ctr, 1)
    [PotAnl(i, :), ~] = a_p_4layer_infinite(No, I0(i, :), dipole_plus(i, :), dipole_minus(i, :), radfactor, cond, strge.PositionOfElectrodes);
    PotAnl(i, :)      = real(PotAnl(i, :));
end
%delete(gcp('nocreate'));

% Remove added paths
warning off; rmpath(genpath(pwd)); warning on;

VoltageDifferenceAnl      = PotAnl(:, 1) - PotAnl(:, 2);
Error2Norm                = norm(VoltageDifferenceAnl - VoltageDifferenceNum)/norm(VoltageDifferenceAnl)
ErrorRDM                  = 0.5*norm(VoltageDifferenceAnl/norm(VoltageDifferenceAnl) - VoltageDifferenceNum/norm(VoltageDifferenceNum))

%% Change back to tests

cd tests

%% Compare with precomputed results

disp("Compare with precomputed results. Compute Differences:");
precomputed = load("result_test_all_individually.mat", 'VoltageDifferenceNum', 'Error2Norm', 'ErrorRDM');
disp("Potential: " + num2str(norm(precomputed.VoltageDifferenceNum - VoltageDifferenceNum)));
disp("Error 2-norm: " + num2str(precomputed.Error2Norm - Error2Norm));
disp("Error RDM: " + num2str(precomputed.ErrorRDM - ErrorRDM));