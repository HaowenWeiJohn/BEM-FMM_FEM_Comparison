function [PC, integralpd] = meshneighborints_P(P, t, normals, Area, Center, RnumberP, ineighborP)
%   Accurate integration for electric potential on neighbor facets
%   Copyright SNM 2018-2021
%
%% Additional documentation by Paul Lunkenheimer
%
% "ineighborP" is to be read as a row vector.
% In row $n$, this vector stores the indices $m$ of all triangles $t_m$ to which $t_n$ will contribute
% exact integration instead of FMM.
%
% "RnumberP" is the number of such neighbors of $t_n$
%

    %% Prepare constants and container for results
    tic 
    N = size(t, 1);
    integralpe      = zeros(RnumberP, N);    % Exact potential integrals to be added
    integralpc      = zeros(RnumberP, N);    % Center-point potential integrals (FMM) to be subtracted

    % When computing the contribution \int_{t_n} \frac{1}{|r-r'|} \, dr'
    % to the potential at $r \in t_m$, which point $r$ should we use? We
    % want one constant result for the whole of $t_m$.
    % Thus we use multiple points $r^m_p \in t_m$ and then average the
    % evaluations. For this we use the following coefficients and weights:
    gauss       = 25;   % Number of evaluation points for averaging
                        % Numbers 1, 4, 7, 13, 25 are permitted 
    if gauss == 1;  [coeffS, weightsS, IndexS]  = tri(1, 1);    end
    if gauss == 4;  [coeffS, weightsS, IndexS]  = tri(4, 3);    end
    if gauss == 7;  [coeffS, weightsS, IndexS]  = tri(7, 5);    end
    if gauss == 13; [coeffS, weightsS, IndexS]  = tri(13, 7);   end
    if gauss == 25; [coeffS, weightsS, IndexS]  = tri(25, 10);  end
    W           = weightsS';

    %% Main loop for analytical double integrals (parallel)
    % This is the loop over columns of the system matrix
    tic
    if(isempty(gcp('nocreate')))
        error('A parallel pool must be initialized prior to running meshneighborints_P');
    end
    
    % Loop over $t_n$ and compute contributions to their neighbors $t_m$
    parfor n = 1:N
        %% Setup necessary values for integral over $t_m$
        % Vertices of $t_n$
        r1      = P(t(n, 1), :);    %   [1x3]
        r2      = P(t(n, 2), :);    %   [1x3]
        r3      = P(t(n, 3), :);    %   [1x3]
        % Neighbors $t_m$ for which we will compute evaluation
        index       = ineighborP(:, n);
        % Container for resulting potentials at neighbors $t_m$
        IP          = zeros(RnumberP, 1);
        % Evaluation points $r^m_p$ on $t_m$
        ObsPoints   = zeros(RnumberP*IndexS, 3);
        for q = 1:RnumberP  % Iterate over number of $t_m$
            num = index(q);
            for p = 1:IndexS    % Iterate over number of $r^m_p$ on one $t_m$
                ObsPoints((q-1)*IndexS + p, :)  = coeffS(1, p)*P(t(num, 1), :) +  coeffS(2, p)*P(t(num, 2), :) +  coeffS(3, p)*P(t(num, 3), :);
            end
        end

        %% Exact integration
        % of \int_{t_n} \frac{1}{|r^m_p-r'|} \, dr'
        [JP, ~] = potint(r1, r2, r3, normals(n, :), ObsPoints);
        % Now take average over all $r^m_p$ per neighbor $t_m$
        for q = 1:RnumberP      
            IP(q)   = sum(W.*JP((q-1)*IndexS + [1:IndexS]), 1);
        end
        integralpe(:, n) = IP;

        %% Center-point electric potential integrals (FMM) to be subtracted
        % Per neighbors $t_m$ that $t_n$ will contribute to, we compute $r_m - r_n$
        % where $r_m, r_n$ are the centers of $t_m, t_n$
        temp    = Center(index, :) - repmat(Center(n, :), RnumberP, 1);
        % $|r_m - r_n|$
        DIST    = sqrt(dot(temp, temp, 2));
        % $A_n \frac{1}{|r_m - r_n|}$
        IPC     = Area(n)./DIST;
        % Set self integral to zero (as does FMM)
        IPC(1)  = 0;
        integralpc(:, n) = IPC;
    end    
    disp([newline 'Integral evaluation time = ' num2str(toc) ' s']);
    
    %%  Define useful sparse matrices EC, PC (for GMRES speed up)
    % Scale by constant $\frac{1}/{4\pi}$
    tic;
    N               = size(t, 1);
    const           = 1/(4*pi);  
    
    ii  = ineighborP;
    jj  = repmat([1:N], RnumberP, 1);
    PC  = sparse(ii, jj, const*(integralpe - integralpc));  % almost symmetric

    integralpd = integralpe - integralpc;
    
    disp([newline 'Correction matrix construction time = ' num2str(toc) ' s']);
end