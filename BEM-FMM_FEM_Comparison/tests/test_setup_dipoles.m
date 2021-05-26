%% Move this test to parent folder to run it
%   Also move "bem1_setup_dipole.m" and "bem1_setup_dipoles.m"
%   These will be called in order to check modified setup of dipole(s)
%   against original scripts
%
%   First need to run "test_preprocess_model.m"

%% Check OS
if ~isunix
    slash = "\";
else
    slash = "/";
end

%% Set up 1 dipole using "setup_dipoles.m"
% Load model data
filename_model = "CombinedMesh_test.mat";
filename_modelP = "CombinedMeshP_test.mat";
p = {};
[p.P, p.t, p.normals, p.Area, p.Center, p.Indicator, p.tissue, p.cond, p.enclosingTissueIdx, p.condin, p.condout, p.contrast, p.eps0, p.mu0, p.tneighbor, p.EC, p.PC] = ...
            load_model(filename_model, filename_modelP);
        
filename_dipoles = "tests" + slash + "dipole.txt";
dipole_length = 4*1e-05;
dipole_tissue_index = 4;

p.D = 10;

[p.Ctr, p.M, p.I0, p.strdipolePplus, p.strdipolePminus, p.dlength, p.strdipolesig, p.strdipoleCurrent, p.strdipolemvector, p.strdipolemcenter, p.strdipolemstrength] = ...
                setup_dipoles(filename_dipoles, dipole_length, dipole_tissue_index, p.cond, p.D);
            


%% Set up same dipole using "bem1_setup_dipole.m" 
clearvars -except p slash;
bem1_setup_dipole;
close all;

%% Compare results        
s = whos;
%s = rmfield(s, 'p');
%s = rmfield(s, 'slash');
for k = 1:length(s)
    name = s(k).name;
    if (length(name)==1 && all(name=='p'))
        continue;
    end
    if (length(name)==5 && all(name=='slash'))
        continue;
    end
    if (length(name)==2 && all(name=='GM'))
        continue;
    end
    if (length(name)==1 && all(name=='R'))
        continue;
    end
    if (length(name)==2 && all(name=='WM'))
        continue;
    end
    if (length(name)==3 && all(name=='arg'))
        continue;
    end
    if (length(name)==1 && all(name=='d'))
        continue;
    end
    if (length(name)==7 && all(name=='dlength'))
        assert(abs(eval(name) -  p.(name)) < 1e4*eps);
        continue;
    end
    if (length(name)==2 && all(name=='f1'))
        continue;
    end
    if (length(name)==2 && all(name=='f2'))
        continue;
    end
    if (length(name)==6 && all(name=='indexg'))
        continue;
    end
    if (length(name)==7 && all(name=='indexg1'))
        continue;
    end
    if (length(name)==7 && all(name=='indexg2'))
        continue;
    end
    if (length(name)==7 && all(name=='indexg3'))
        continue;
    end
    if (length(name)==6 && all(name=='indexw'))
        continue;
    end
    if (length(name)==7 && all(name=='indexw1'))
        continue;
    end
    if (length(name)==7 && all(name=='indexw2'))
        continue;
    end
    if (length(name)==7 && all(name=='indexw3'))
        continue;
    end
    if (length(name)==1 && all(name=='m'))
        continue;
    end
    if (length(name)==3 && all(name=='str'))
        continue;
    end
    if (length(name)==2 && all(name=='t0'))
        continue;
    end
    if (length(name)==4 && all(name=='temp'))
        continue;
    end
    if (length(name)==5 && all(name=='theta'))
        continue;
    end
    if (length(name)==6 && all(name=='tissue'))
        eval(name)
        p.(name)
        continue;
    end
    if (length(name)==14 && all(name=='tissue_to_plot'))
        continue;
    end
    disp(name);
    assert(all(all(eval(name) ==  p.(name))));
end



