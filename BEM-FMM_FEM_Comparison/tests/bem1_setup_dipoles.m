%   This script creates the cortical dipole layer (multiple dipoles)
%
%   Copyright SNM/WAW 2018-2020

%% Add paths
if ~isunix
    s = pwd;
    addpath(strcat(s, '\Engine'));
    addpath(strcat(s, '\tests'));
else
    s = pwd;
    addpath(strcat(s, '/Engine'));
    addpath(strcat(s, '/tests'));
end
clear s;

%% Load model
filename_model = "CombinedMesh_test.mat";
filename_modelP = "CombinedMeshP_test.mat";
load(filename_model, 'P', 't', 'normals', 'Area', 'Center', 'Indicator', 'tissue', 'cond', 'enclosingTissueIdx', 'condin', 'condout', 'contrast', 'eps0', 'mu0');
load(filename_modelP, 'tneighbor', 'EC', 'PC');
clear filename_model filename_modelP;

%% Hardcode Centers and Orientation for Comparison
Ctr = zeros(20, 3);
temp = zeros(20, 3);
% 10 radial dipoles
Ctr(1:10, :) = ...
    [35.2763689016264	-2.57993017776734	66.7021120909542
    26.0812979186662	58.4800984517253	40.0011747821773
    -58.7655385107809	-9.21086533117258	46.4975423370841
    -38.2683022799001	-10.9568285380806	64.1539940222064
    16.5989779718782	72.1707297385933	-14.7006700285406
    70.0458965823608	27.5813683657766	5.7524334889289
    65.1570434813793	33.3425409954718	-18.5225442294227
    -58.2394061805888	28.9710146456914	38.3288648166304
    66.1886830208771	-26.3340970749591	-25.0164660015855
    37.0670111623556	36.6497632362896	-54.6175936691989];
temp(1:10, :) = Ctr(1:10, :)./sqrt(sum(Ctr(1:10, :).^2, 2));

% 10 tangential dipoles
Ctr(11:20, :) = ...
    [-66.1336335024612	31.7471519034479	-17.8524750603816	
    -35.5535841948427	21.3939624580307	-63.0752805879182	
    -20.2283885308647	2.4855783549136	-72.6972090089156	
    39.8886919045368	-48.4317883256639	-41.9940964627532	
    -16.4165015469782	66.6115449153564	-31.5222232868737	
    24.6279720878085	12.2932275059592	-70.3035528855324	
    -67.7439864360699	-9.87414180462654	31.8355716985215	
    -42.4434771666759	2.74557694535596	62.3799892051849	
    -71.927969473119	-22.8813841323455	-1.74913340317495	
    -50.5343613870906	-56.0828575331148	-1.11418585490723];	
temp(11:20, :) = ...    
    [-0.0728714189476909	0.369195222386579	0.926490498638531
    -0.879859380554735	-0.0821854408215915	0.468073737531597
    -0.828470750463066	0.502275911521457	0.247699665588094
    -0.688850674379014	0.0591124888789937	-0.722489074011579
    -0.763760871231798	-0.420136365814883	-0.490055880176035
    -0.163374737894742	-0.960520313498607	-0.225187527128751
    0.248606487521355	0.639600248033905	0.727396959766303
    -0.791252455075087	0.266994615988862	-0.550121284237038
    -0.0829278743513777	0.332499878849458	-0.939450263835536
    -0.446522903971554	0.418063274454369	-0.79110074881857];
temp(11:20, :) = temp(11:20, :)./sqrt(sum(temp(11:20, :).^2, 2));
Ctr = Ctr*1e-03;

%%   Multiple dipole example
I0 = 1*1e-3/(4*1e-05);  %   source current, A
d  = 4*1e-05;           %   finite-dipole length, m
s  = 2.5e-3;            %   spacing from GM, m       
R  = 0.005;             %   radius of the enclosing sphere in m   
%Ctr= [0 0 75.5e-3];     %   position of enclosing sphere in m
GM              = load('meshsphere6_78.mat');   %   in mm
GM.P            = 1e-3*GM.P;                    %   now in m
GM.Center       = meshtricenter(GM.P, GM.t);    %   base for the dipole layer
%   Indexes into GM triangles strictly within the sphere
indexg1         = find( (GM.P(GM.t(:, 1), 1)-Ctr(1)).^2 + (GM.P(GM.t(:, 1), 2)-Ctr(2)).^2 + (GM.P(GM.t(:, 1), 3)-Ctr(3)).^2 < R^2);
indexg2         = find( (GM.P(GM.t(:, 2), 1)-Ctr(1)).^2 + (GM.P(GM.t(:, 2), 2)-Ctr(2)).^2 + (GM.P(GM.t(:, 2), 3)-Ctr(3)).^2 < R^2);
indexg3         = find( (GM.P(GM.t(:, 3), 1)-Ctr(1)).^2 + (GM.P(GM.t(:, 3), 2)-Ctr(2)).^2 + (GM.P(GM.t(:, 3), 3)-Ctr(3)).^2 < R^2);
indexg          = intersect(intersect(indexg1, indexg2), indexg3); 
M               = length(Ctr);
strdipolePplus              = Ctr + 0.00002*temp;
strdipolePminus             = Ctr - 0.00002*temp;
clear temp;
strdipolesig                = repmat(cond(4), 2*M, 1);
clear strdipoleCurrent; 
strdipoleCurrent(1:M, 1)    = +I0;
strdipoleCurrent(M+1:2*M, 1)= -I0;

NoDipoles = size(strdipolePplus, 1)

%%   Magnetic dipole subdivision (optional)
D = 1;                        %   number of smaller subdipoles
strdipolemvector   = zeros(D*M, 3);
strdipolemcenter   = zeros(D*M, 3);
strdipolemstrength = zeros(D*M, 1);
for m = 1:M
    temp = (1/D)*(strdipolePplus(m, :) - strdipolePminus(m, :));
    for d = 1:D 
        arg = d+D*(m-1);
        strdipolemvector(arg, :)     = temp;
        strdipolemcenter(arg, :)     = strdipolePminus(m, :) + (d-1/2)*temp;
        strdipolemstrength(arg, :)   = strdipoleCurrent(m);                  
    end
end

%%  Plot and check correct position
WM          = load('meshsphere6_73.mat');    %   in mm
WM.P        = 1e-3*WM.P;                    %   now in m
WM.Center   = meshtricenter(WM.P, WM.t);

indexw1 = find( (WM.P(WM.t(:, 1), 1)-Ctr(1)).^2 + (WM.P(WM.t(:, 1), 2)-Ctr(2)).^2 + (WM.P(WM.t(:, 1), 3)-Ctr(3)).^2 < R^2);
indexw2 = find( (WM.P(WM.t(:, 2), 1)-Ctr(1)).^2 + (WM.P(WM.t(:, 2), 2)-Ctr(2)).^2 + (WM.P(WM.t(:, 2), 3)-Ctr(3)).^2 < R^2);
indexw3 = find( (WM.P(WM.t(:, 3), 1)-Ctr(1)).^2 + (WM.P(WM.t(:, 3), 2)-Ctr(2)).^2 + (WM.P(WM.t(:, 3), 3)-Ctr(3)).^2 < R^2);
indexw  = intersect(intersect(indexw1, indexw2), indexw3); 

% Plot dipole(s) between WM and GM
f1 = figure;
str.EdgeColor = 'k'; str.FaceColor = [0 1 1]; str.FaceAlpha = 1.0; 
bemf2_graphics_base(WM.P, WM.t(indexw, :), str);
str.EdgeColor = 'k'; str.FaceColor = [0.5 0.5 0.5]; str.FaceAlpha = 1.0; 
bemf2_graphics_base(GM.P, GM.t(indexg, :), str);
bemf1_graphics_dipole(strdipolePplus, strdipolePminus, strdipoleCurrent, 4) 
axis 'equal';  axis 'tight';   
daspect([1 1 1]);
set(gcf,'Color','White');
camlight; lighting phong;
xlabel('x'); ylabel('y'); zlabel('z');
view(0, 0);

% Plot dipole(s) above WM
f2 = figure;
tissue_to_plot = 'WM';
t0 = t(Indicator==find(strcmp(tissue, tissue_to_plot)), :);    % (change indicator if necessary: 1-skin, 2-skull, etc.)
str.EdgeColor = 'none'; str.FaceColor = [0 1 1]; str.FaceAlpha = 1.0; 
bemf2_graphics_base(P, t0, str);
bemf1_graphics_dipole(strdipolePplus, strdipolePminus, strdipoleCurrent, 4) 
axis 'equal';  axis 'tight';   
daspect([1 1 1]);
set(gcf,'Color','White');
camlight; lighting phong;
xlabel('x'); ylabel('y'); zlabel('z');
view(0, 0);