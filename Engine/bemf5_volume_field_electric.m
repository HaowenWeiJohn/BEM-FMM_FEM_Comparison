function E = bemf5_volume_field_electric(Points, c, P, t, Center, Area, normals, R, prec, planeABCD)
%   Computes electric field for an array Points anywhere in space (line,
%   surface, volume). This field is due to surface charges at triangular
%   facets only. Includes accurate neighbor triangle integrals for
%   points located close to a charged surface.   
%   R is the dimensionless radius of the precise-integration sphere
%
%   Copyright SNM 2017-2020
%   R = is the local radius of precise integration in terms of average triangle size
%
%% Additional documentation by Paul Lunkenheimer
%
% Exact integration will be used for all neighbors in a radius of "R" times
% average triangle area.
% If "planeABCD" is given, only the neighbors in this plane are used for
% exact integration.
%
    

    if(nargin < 10)
        planeABCD = [];
    end
    
    %% Electric field at "Points" computed using FMM
    %   FMM 2019    
    srcinfo.sources = Center';                      %   source points
    targ            = Points';                      %   target points   
    pg      = 0;                                    %   nothing is evaluated at sources
    pgt     = 2;                                    %   field and potential are evaluated at target points
    srcinfo.charges = c.'.*Area';                   %   charges
    U               = lfmm3d(prec, srcinfo, pg, targ, pgt);
    E               = -U.gradtarg'/(4*pi);  

    %% Choose a set of triangles that are possible near neighbors for exact integration
    Size  = mean(sqrt(Area));
    if(isempty(planeABCD))  % Either all triangles are possible near neighbors
        eligibleTriangles = 1:size(t, 1);
    else    % Or only triangles that are close to the hyperplane "planeABCD"
        % If hyperplane is given by $\{x \in R^3 | <a,x> = z \}$,
        % it is defined by the affine linear map $\psi(x) = <a,x> - z$,
        % we allow triangles that suffice $|\psi(r_n)| \leq R*Size$.
        % Here $r_n$ is the center of triangle $t_n$
        % and $Size$ is the average triangle area.
        d1 = abs(planeABCD(1)*Center(:,1) + planeABCD(2)*Center(:,2) + planeABCD(3)*Center(:,3) + planeABCD(4));
        d2 = norm(planeABCD(1:3));
        d = d1./d2;
        eligibleTriangles = find(d <= R*Size);
    end
    
    % For every "eligibleTriangle" $t_n$, we compute the neighbors $t_m$
    % that $t_n$ will contribute exact integration to instead of FMM
    ineighborlocal   = rangesearch(Points, Center(eligibleTriangles, :), R*Size, 'NSMethod', 'kdtree');  
    
    %% Loop over $t_n$ and compute contributions to their neighbors $t_m$
    const = 4*pi;
    for j = 1:length(eligibleTriangles)
        n = eligibleTriangles(j);   % The triangle $t_n$
        index = ineighborlocal{j};  % and the triangles $t_m$ that $t_n$ will contribute to
        if ~isempty(index)
            %% Center-point electric potential integrals (FMM) to be subtracted
            % Per neighbors $t_m$ that $t_n$ will contribute to, we compute $r_m - r_n$
            % where $r_m, r_n$ are the centers of $t_m, t_n$
            temp = Points(index, :) - repmat(Center(n, :), length(index), 1);
            % $|r_m - r_n|$
            DIST = sqrt(dot(temp, temp, 2));
            % $A_n \nabla_r \frac{1}{|r_m - r_n|} = A_n \frac{r_m - r_n}{|r_m - r_n|^3}$
            I = Area(n)*temp./repmat(DIST.^3, 1, 3);
            % Subtract FMM result $c_n A_n \frac{r_m - r_n}{|r_m - r_n|^3}$
            E(index, :) = E(index, :) - c(n)*I/const;
            
            %% Exact integration
            % of neighbor integral $\int_{t_n} \nabla_r \frac{1}{\abs{r_m-r'}} \, dr'$,
            % where $r_m$ is center of $t_m$.
            % Vertices of the neighbors $t_n$
            r1      = P(t(n, 1), :);
            r2      = P(t(n, 2), :);
            r3      = P(t(n, 3), :);
            % Exact integration
            I       = potint2(r1, r2, r3, normals(n, :), Points(index, :));
            % Add exact result $\int_{t_n} - \nabla_r \frac{1}{4\pi} \frac{c_n}{\abs{r_m-r'}} \, dr'$
            E(index, :)= E(index, :) + (- c(n)*I/const);           
        end    
    end
end