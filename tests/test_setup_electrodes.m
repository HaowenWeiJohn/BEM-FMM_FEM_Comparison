%% Move this test to parent folder to run it

%% Check OS
if ~isunix
    slash = "\";
else
    slash = "/";
end

%% Set up electrodes and refine mesh using "setup_electrodes.m"
filename_mesh = "tests" + slash + "four_layer_surf_from_tets_3.mat";
filename_electrodes = "tests" + slash + "electrodes.txt";

RadE = 5;

filename_output_mesh = "four_layer_surf_from_tets_3_refined.mat";
filename_output_electrodes = "four_layer_surf_from_tets_3_electrodes.mat";

setup_electrodes(filename_mesh, filename_electrodes, RadE, filename_output_mesh, filename_output_electrodes);

% %% And compare with results created with original "Model/model01_main_script"
% clearvars -except p slash;
% 
% filename_original = "tests" + slash + "CombinedMesh_3.mat";
% 
% s = load(filename_original);
% if isfield(s, 'name')
%     s = rmfield(s, 'name');
% end
% if isfield(p, 'eps0')
%     p = rmfield(p, 'eps0');
% end
% if isfield(p, 'mu0')
%     p = rmfield(p, 'mu0');
% end
% 
% names = fieldnames(s);
% for i=1:length(names)
%     name = names{i}
%     assert(isfield(p, name));
%     if (length(name)==6 && all(name=='tissue'))
%        disp('Please compare "by hand":');
%        disp(s.tissue);
%        disp(p.tissue);
%        continue;
%     end
%     if (length(name)==1 && all(name=='P'))
%        continue;
%     end
%     if (length(name)==length('enclosingTissueIdx') && all(name=='enclosingTissueIdx'))
%         assert((all(s.(name)(:, 1)==p.(name), 'all')));   
%         continue;
%     end
%     if (length(name)==1 && all(name=='t'))
%         for k=unique(p.Indicator)'
%             assert(sum(s.Indicator==k, 1)==sum(p.Indicator==k, 1));
%             assert(all(abs(p.P(p.t(p.Indicator==k, 1), :) - s.P(s.t(s.Indicator==k, 1), :)) < 1e4*eps, 'all'));
%             assert(all(abs(p.P(p.t(p.Indicator==k, 2), :) - s.P(s.t(s.Indicator==k, 2), :)) < 1e4*eps, 'all'));
%             assert(all(abs(p.P(p.t(p.Indicator==k, 3), :) - s.P(s.t(s.Indicator==k, 3), :)) < 1e4*eps, 'all'));
%         end
%         continue;
%     end
%     if (size(p.(name), 1)==size(p.t, 1))
%         for k=unique(p.Indicator)'
%             assert(sum(s.Indicator==k, 1)==sum(p.Indicator==k, 1));
%             assert((all(s.(name)(s.Indicator==k)==p.(name)(p.Indicator==k), 'all')));
%         end
%         continue;
%     end
%     assert((all(s.(name)==p.(name), 'all')));
% end
% 
% %% Output for user
% disp(' ');
% disp(' ');
% disp('Please check by hand if tissue names are the same!');
% disp('All tests ran succesfully');
% 
