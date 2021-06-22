function execute_all(filename_mesh, filename_electrodes, filename_tissue, filename_cond, numThreads, TnumberE, GnumberE, RnumberP, iter, relres, prec_charge, weight, filename_dipoles, R, prec_potential, filename_results)
% This function calls all the functions necessary to:
%   setup the electrodes on the mesh
%   preprocess the model
%   solve the linear equation for the charge distribution
%   compute the potentials at the electrodes resulting from given dipoles
%
% As we assume that multiple electrodes are loaded, the first one is chosen
% as reference electrode and we iterate over all the others to compute the
% potentials evoked by the dipoles at different electrodes.

    %% Check OS
    if ~isunix
        slash = "\";
    else
        slash = "/";
    end
    
    %% Set up electrodes and refine mesh using "setup_electrodes.m"
    % Takes mesh and electrode positions.
    % In a radius of 1.5*AverageEdgeLength all triangles around an electrode
    % are marked as belonging to this electrode patch.
    % The mesh is refined in this region around each electrode.
    timer_setup_electrodes = tic;
    
    [P, t, normals, Indicator, IndicatorElectrodes, strge] = ...
        setup_electrodes(filename_mesh, filename_electrodes);

    % Scale from mm to m
    strge.PositionOfElectrodes = strge.PositionOfElectrodes*1e-3;
    strge.RadiusOfElectrodes   = strge.RadiusOfElectrodes*1e-3;
    
    time_setup_electrodes = toc(timer_setup_electrodes);

    %% Preprocess model using "preprocess_model.m"
    % The mesh, tissue names and tissue conductivities are loaded.
    % Additional data like the area per triangle, the triangle
    % centers and conducivity (in, out and contrast) per triangle is computed.
    % Correction matrices EC and PC for electric field and potential are
    % computed. These correct the FMM computations by more exact integration of
    % electric activity induced by near neighbors.
    % The matrix M is a preconditioner for the iterative solver used in the
    % electrode regions.
    timer_preprocess_model = tic;

    [P, t, normals, Area, Center, Indicator, tissue, enclosingTissueIdx, cond, condin, condout, contrast, eps0, mu0, EC, PC, M_total, integralpd, ineighborE, ineighborP, ElectrodeIndexes_global_total, ElectrodeIndexes_local_total, V_total] = ...
                preprocess_model(P, t, normals, Indicator, IndicatorElectrodes, filename_cond, filename_tissue, numThreads, TnumberE, GnumberE, RnumberP);
            
    time_preprocess_model = toc(timer_preprocess_model);
    
    %% Iterate over all electrodes to compute potentials evoked by the dipoles
    time_charge_engine                = cell(strge.NumberOfElectrodes, 1);
    time_compute_dipole_potential     = cell(strge.NumberOfElectrodes, 1);
    timer_solve_forward_problem_total = tic;
    
    % Input
    % Extract like this in order to optimize parfor loop
    ElectrodeIndexes_global_total1 = ElectrodeIndexes_global_total{1};
    ElectrodeIndexes_local_total1  = ElectrodeIndexes_local_total{1};
    V1                             = V_total(1:length(ElectrodeIndexes_global_total{1}), :);
    
    % Containers for output
    c                    = cell(strge.NumberOfElectrodes, 1);
    resvec               = cell(strge.NumberOfElectrodes, 1);
    electrodeCurrents    = cell(strge.NumberOfElectrodes, 1);
    En                   = cell(strge.NumberOfElectrodes, 1);
    dipole_ctr           = cell(strge.NumberOfElectrodes, 1);
    dipole_moment        = cell(strge.NumberOfElectrodes, 1);
    dipole_n             = cell(strge.NumberOfElectrodes, 1);
    VoltageDifference    = cell(strge.NumberOfElectrodes, 1);
    
    parpool(numThreads);
    parfor i=2:strge.NumberOfElectrodes
        %% Adjust variables to current electrode pair
        ElectrodeIndexes_global    = cell(2, 1);
        ElectrodeIndexes_global{1} = ElectrodeIndexes_global_total1;
        ElectrodeIndexes_local     = cell(2, 1);
        ElectrodeIndexes_local{1}  = ElectrodeIndexes_local_total1
        % Indices of triangles belonging to electrode i in whole mesh
        ElectrodeIndexes_global{2} = ElectrodeIndexes_global_total{i};
        % Indices of triangles belonging to electrode i in reduced mesh consisting of only electrode triangles
        ElectrodeIndexes_local{2} = ElectrodeIndexes_local_total{i};
        % Set voltages
        V = [V1; V_total(ElectrodeIndexes_local_total{i}, :)];
        % Indices of triangles belonging to electrode 1 and i in reduced mesh consisting of only electrode triangles
        index_loc = [1 : length(ElectrodeIndexes_global{1}), ElectrodeIndexes_local_total{i}];
        % Project preconditioner matrix M
        M = M_total(index_loc, index_loc);
        
        %% Compute solution
        % This solution is only linked to the electrodes, but not yet to specific
        % dipoles.
        timer_charge_engine = tic;

        [c{i}, resvec{i}, electrodeCurrents{i}, En{i}] = ...
            charge_engine(normals, Area, Center, condin, contrast, EC, PC, M, ElectrodeIndexes_global, ElectrodeIndexes_local, V, iter, relres, prec_charge, weight);
    
        time_charge_engine{i} = toc(timer_charge_engine);
    
        %% Compute potential at electrodes generated by dipoles using reciprocity
        timer_compute_dipole_potential = tic;

        [dipole_ctr{i}, dipole_moment{i}, dipole_n{i}, VoltageDifference{i}] = ...
            compute_dipole_potential(filename_dipoles, c{i}, P, t, Center, Area, normals, electrodeCurrents{i}, R, prec_potential);
    
        time_compute_dipole_potential{i} = toc(timer_compute_dipole_potential);
        disp(['Electrode ' num2str(i) ' computed']);
    end
    delete(gcp('nocreate'));
    
    VoltageDifference{1} = zeros(dipole_n{2}, 1);

    time_solve_forward_problem_total = toc(timer_solve_forward_problem_total);
    
    %% Save results
    save(filename_results, ...
            'IndicatorElectrodes', 'strge', ...
            'P', 't', 'normals', 'Area', 'Center', 'Indicator', 'tissue', 'enclosingTissueIdx', 'cond', 'condin', 'condout', 'contrast', 'eps0', 'mu0', 'EC', 'PC', 'M_total', 'integralpd', 'ineighborE', 'ineighborP', 'ElectrodeIndexes_global_total', 'ElectrodeIndexes_local_total', 'V_total', ...
            'c', 'resvec', 'electrodeCurrents', 'En', ...
            'dipole_ctr', 'dipole_moment', 'dipole_n', 'VoltageDifference', ...
            'time_setup_electrodes', 'time_preprocess_model', 'time_solve_forward_problem_total', 'time_charge_engine', 'time_compute_dipole_potential', ...
            '-v7.3');
 
end