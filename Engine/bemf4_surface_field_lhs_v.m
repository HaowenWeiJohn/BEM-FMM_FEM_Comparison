function LHS = bemf4_surface_field_lhs_v(c, Center, Area, contrast, normals, M, EC, PC, indexe, weight, condin, prec)   
%   Computes the left hand side of the charge equation for surface charges
%
%   Copyright SNM 2017-2020

%   LHS is the user-defined function of c equal to c - Z_times_c which is
%   exactly the left-hand side of the matrix equation Zc = b
    %timer_bemf4_surface_field_lhs_v = tic;
    disp("Computation of left-hand side of linear equation");

    [P0, E0]    = bemf4_surface_field_electric_plain(c, Center, Area, prec);%   Plain FMM result    
    correction  = contrast.*(EC*c);                                          %   Correction of plain FMM result
    
    LHS         = +c - 2*correction ...                                     %   This is the dominant (exact) matrix part and the "undo" terms for center-point FMM
                     - 2*(contrast.*sum(normals.*E0, 2));                   %   This is the full center-point FMM part
                     
                           
    correctionP  = PC*c;                                                    %   Correction of plain FMM result for potential
    P            = P0 + correctionP;                                        %   Exact results for potential
    LHS(indexe)  = M*P(indexe);                                             %   LHS for potential with preconditioner
    % M should be (using the full preconditioner) exactly the inverse of the linear operator applied to c(indexe)
    % So this should be (ignoring precision errors) simply: c(indexe)
    
%   Normal field just inside
    En          = -c(indexe)/2 + correction(indexe) ...                     %   This is the dominant (exact) matrix part and the "undo" terms for center-point FMM
                       + sum(normals(indexe, :).*E0(indexe, :), 2);         %   This is the full center-point FMM part      

%   Total current (normalized)
    I = sum(En.*Area(indexe).*condin(indexe))/sum(Area(indexe).*condin(indexe));    
    LHS         = LHS + weight*I;                                               %   Adding current conservation law
    
    %toc(timer_bemf4_surface_field_lhs_v);
end
