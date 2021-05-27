%% Move this test to parent folder to run it
%   Does not varify "EC" and "PC"
%   (but code for computing these variables is exactly the same as in
%   original)

%% Check OS
if ~isunix
    slash = "\";
else
    slash = "/";
end

%% Preprocess a mesh using "preprocess_model.m"
filename_mesh = "tests" + slash + "four_layer_surf_from_tets_3.mat";
filename_tissue = "tests" + slash + "fls_tissue.tiss";
filename_cond = "tests" + slash + "fls_conductivities.cond";

filename_output = "CombinedMesh_test.mat";
filename_outputP = "CombinedMeshP_test.mat";

numThreads = 4;     % Might not be valid on your machine!
RnumberE = 128;
RnumberP = 128;

p = preprocess_model(filename_mesh, filename_cond, filename_tissue, filename_output, filename_outputP, numThreads, RnumberE, RnumberP);

%% And compare with results created with original "Model/model01_main_script"
clearvars -except p slash;

filename_original = "tests" + slash + "CombinedMesh_3.mat";

s = load(filename_original);
if isfield(s, 'name')
    s = rmfield(s, 'name');
end
if isfield(p, 'eps0')
    p = rmfield(p, 'eps0');
end
if isfield(p, 'mu0')
    p = rmfield(p, 'mu0');
end

names = fieldnames(s);
for i=1:length(names)
    name = names{i}
    assert(isfield(p, name));
    if (length(name)==6 && all(name=='tissue'))
       disp('Please compare "by hand":');
       disp(s.tissue);
       disp(p.tissue);
       continue;
    end
    if ~(all(s.(name)==p.(name), 'all'))
       disp("Unequal: " + name); 
    end
end

%% Output for user
disp(' ');
disp(' ');
disp('Please check by hand if tissue names are the same!');
disp('All tests ran succesfully');

