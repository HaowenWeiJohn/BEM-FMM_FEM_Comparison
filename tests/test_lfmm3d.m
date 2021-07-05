%% Script for running a simple lfmm3d setup to check validity of results
% This may be used to check if you version of lfmm3d runs correctly

%% Change to parent folder
cd ..

%% Add paths
if ~isunix
    addpath(strcat(pwd, '\Engine'));
    addpath(strcat(pwd, '\io'));
    slash = "\";
else
    addpath(strcat(pwd, '/Engine'));
    addpath(strcat(pwd, '/io'));
    slash = "/";
end

%% Load data
filename_charges = "tests" + slash + "charges.txt";
[P, C, n_of_charges] = read_charges(filename_charges);

%% Run FMM

disp("Compute electric potential and field using FMM");
pg              = 2;                         % potential and field are evaluated
srcinfo.sources = P';                        % source/target points
srcinfo.charges = C';                        % charges
prec            = 1e-2;                      % FMM precision
U               = lfmm3d(prec, srcinfo, pg); % FMM

%% Compute electric field and potential directly

disp("Compute electric potential and field directly");
diff_x = P(:, 1) - P(:, 1)';                    % r_x-r'_x
diff_y = P(:, 2) - P(:, 2)';                    % r_y-r'_y
diff_z = P(:, 3) - P(:, 3)';                    % r_z-r'_z
dist = sqrt(diff_x.^2 + diff_y.^2 + diff_z.^2); % |r-r'|
ker = 1./dist;                                  % 1/|r-r'|
ker(1:size(P, 1)+1:end) = 0;                    % Set r'=r terms to 0
pot = (ker*C);                                  % Compute sum_i c_i/|r-r_i'|
ker_x = diff_x.*ker.^3;                         % r_x-r'_x/|r-r'|^3
ker_x(1:size(P, 1)+1:end) = 0;                  % Set r'=r terms to 0
ker_y = diff_y.*ker.^3;                         % r_y-r'_y/|r-r'|^3
ker_y(1:size(P, 1)+1:end) = 0;                  % Set r'=r terms to 0
ker_z = diff_z.*ker.^3;                         % r_z-r'_z/|r-r'|^3
ker_z(1:size(P, 1)+1:end) = 0;                  % Set r'=r terms to 0
field = zeros(size(P, 1), 3);
field(:, 1) = ker_x*C;                          % Compute sum_i c_i (r_x-r'_x)/|r-r_i'|^3
field(:, 2) = ker_y*C;                          % Compute sum_i c_i (r_y-r'_y)/|r-r_i'|^3
field(:, 3) = ker_z*C;                          % Compute sum_i c_i (r_z-r'_z)/|r-r_i'|^3

%% Compare results

disp("Compare the results");
diff_pot_2norm = norm(U.pot' - pot);
C_2norm = norm(C);
diff_field_2norm = sqrt(sum(sum((-U.grad' - field).^2)));
disp("2-norm of potential difference: " + num2str(diff_pot_2norm));
disp("...normalized by 2-norm of charges: " + num2str(diff_pot_2norm./C_2norm));
disp("2-norm of field difference: " + num2str(diff_field_2norm));
disp("...normalized by 2-norm of charges: " + num2str(diff_field_2norm./C_2norm));

%% Remove added paths
warning off; rmpath(genpath(pwd)); warning on;

%% Change back to tests
cd tests