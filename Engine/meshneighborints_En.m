function [EC] = meshneighborints_En(P, t, normals, Area, Center, RnumberE, ineighborE, TnumberE)
%   Accurate integration for electric field on neighbor facets using the solid angle approach
%   Copyright WAW/SNM 2020-2021
%
%% Additional documentation by Paul Lunkenheimer
%
% "ineighborE" is to be read as a row vector.
% In row $n$, this vector stores the indices $m$ of all triangles $t_m$ to which $t_n$ will contribute
% exact integration instead of FMM.
%
% "RnumberE" is the number of such neighbors of $t_n$
%
    %% Prepare constants and container for results
    N = size(t, 1);
    % Center-point normal electric field integrals (FMM) to be subtracted
    integralxc      = zeros(RnumberE, N);   % x-coordinate of E-field
    integralyc      = zeros(RnumberE, N);   % y-coordinate of E-field
    integralzc      = zeros(RnumberE, N);   % z-coordinate of E-field
    % Exact normal electric field to be added
    integrale = zeros(N, RnumberE);
    
    % When computing the contribution \int_{t_n} <\eta_m, \int_{t_m} \nabla_r \frac{1}{\abs{r-r'}} \, dr > \, dr'
    % to the normal electric field, we approximate the outer integral using
    % a quadrature rule.
    % Thus we use multiple points $r^n_p \in t_n$ and compute a weighted sum
    % of these evaluations.
    % For this we use the following coefficients and weights:
    gauss       = 25;   % Number of integration points
                        % Numbers 1, 4, 7, 13, 25 are permitted 
    if gauss == 1;  [coeffS, weightsS, IndexS]  = tri(1, 1); end
    if gauss == 4;  [coeffS, weightsS, IndexS]  = tri(4, 3); end
    if gauss == 7;  [coeffS, weightsS, IndexS]  = tri(7, 5); end
    if gauss == 13; [coeffS, weightsS, IndexS]  = tri(13, 7); end
    if gauss == 25; [coeffS, weightsS, IndexS]  = tri(25, 10); end
    if gauss == 0;  [coeffS, weightsS, IndexS]  = tri(50);      end
    

    %% Main loop for analytical double integrals (parallel)
    % This is the loop over columns of the system matrix
    tic
    if(isempty(gcp('nocreate')))
        error('A parallel pool must be initialized prior to running meshneighborints_En');
    end
    
    % Loop over $t_n$ and compute contributions to their neighbors $t_m$
    parfor n = 1:N
        %% Setup necessary values for double-integral
        % Evaluation points $r^m_p$ on $t_m$
        ObsPoints = zeros(IndexS, 3);
        for p = 1:IndexS
            ObsPoints(p, :)  = coeffS(1, p)*P(t(n, 1), :) +  coeffS(2, p)*P(t(n, 2), :) +  coeffS(3, p)*P(t(n, 3), :);
        end
        % Vertices of the neighbors $t_m$ that $t_n$ will contribute to
        index = ineighborE(:,n);
        r1 = P(t(index, 1), :);
        r2 = P(t(index, 2), :);
        r3 = P(t(index, 3), :);
        
        %% Exact integration
        % of <\eta_m, \int_{t_m} \nabla_r \frac{1}{\abs{r-r'^n_p}} \, dr >
        Int_temp = potint4b(r1, r2, r3, ObsPoints);
        Int_temp(:,1) = 0;  % Kill self-term
        % Now compute weighted sum over all $r^n_p$ on $t_n$
        % This quadrature rule is still normalized to $1$ and needs to be
        % multiplied by $A_n$
        Int = weightsS*Int_temp;
        if TnumberE>0
            Int(1) = mean(Int(2:TnumberE));   % Self integral, SNM, 06/09/21
        end
        integrale(n, :) = Int;       
           
        %% Center-point normal electric field integrals (FMM) to be subtracted
        % Per neighbors $t_m$ that $t_n$ will contribute to, we compute $r_m - r_n$
        % where $r_m, r_n$ are the centers of $t_m, t_n$
        temp    =  Center(index, :) - repmat(Center(n, :), RnumberE, 1);
        % $|r_m - r_n|$
        DIST    = sqrt(dot(temp, temp, 2));
        % $A_m \nabla_r \frac{1}{|r_m - r_n|} = A_m \frac{r_m - r_n}{|r_m - r_n|^3}$
        I       = Area(n)*temp./repmat(DIST.^3, 1, 3);
        % Set self integral to zero (as does FMM)
        I(1, :) = 0;
        integralxc(:, n) = I(:, 1);
        integralyc(:, n) = I(:, 2);
        integralzc(:, n) = I(:, 3);
        
    end
    disp([newline 'Integral evaluation time = ' num2str(toc) ' s']);    
    
    tic;
    %% Weight quadrature rule integral with $A_n$ and divide by $A_m$
    % $\frac{1}{A_m} A_n \sum_p w_p <\eta_m, \int_{t_m} \nabla_r \frac{1}{\abs{r-r'^n_p}} \, dr>$
    area_neighbor = Area(transpose(ineighborE));
    area_self = repmat(Area, 1, RnumberE);
    integrale = integrale .* area_self ./ area_neighbor;
    
    %% Compute NORMAL electric field from center-point electric field
    % $<\eta_m, A_m \nabla_r \frac{1}{|r_m - r_n|}>$
    integralc       = zeros(RnumberE, N);    
    for n = 1:N            
        index = ineighborE(:, n);                               
        integralc(:, n)  =       +(integralxc(:, n).*normals(index, 1) + ...
                                   integralyc(:, n).*normals(index, 2) + ...
                                   integralzc(:, n).*normals(index, 3));
    end
    
    %% Define useful sparse matrices EC, PC (for GMRES speed up)
    % Scale by constant $\frac{1}/{4\pi}$
    const           = 1/(4*pi);
    
    ii  = ineighborE;
    jj  = repmat([1:N], RnumberE, 1);
    EC  = sparse(ii, jj, const*(-integralc + transpose(integrale)));    % almost symmetric
    
    disp([newline 'Correction matrix construction time = ' num2str(toc) ' s']);
    
    %Uncomment the following line to save debugging information
    %save('integrals_test_TES', 'integrale', 'integralc', 'EC');
end