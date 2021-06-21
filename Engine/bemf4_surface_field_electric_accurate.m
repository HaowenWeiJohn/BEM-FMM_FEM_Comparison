function [En] = bemf4_surface_field_electric_accurate(c, Center, Area, normals, EC, prec, targetIndex)
%   Computes the normal electric field just inside
%
%   Copyright SNM 2017-2020
%
%% Additional documentation by Paul Lunkenheimer
%
% Allows evaluation of normal electric field at specific points
%
% (more efficient than global evaluation at all triangle centers).
% Pass complete "c", "Center", "Area", "normals", "EC".
% "prec" is requested precision of FMM.
% If argument "targetIndex" is passed, normal electric field will only be
% computed at the triangle centers specified as boolean vector in
% "targetIndex".
% If argument is not passed, normal electric field will be computed at all
% triangle centers "Center".

    %% Electric field computed using FMM
    % If targets are given, we only compute field at targets
    srcinfo.sources = Center';      % source points
    srcinfo.charges = (c.*Area)';   % real charges (as charge at center stretches over whole triangle)
    if(nargin == 6)
        targ = Center(targetIndex)';  % points of evaluation
        pg   = 0;        % nothing is evaluated at sources
        pgt  = 2;        % field and potential are evaluated at target points
        U    = lfmm3d(prec, srcinfo, pg, targ, pgt);
        E    = -U.gradtarg'/(4*pi);
    else
        pg   = 2;        % potential and field are evaluated at targets = source points
        U    = lfmm3d(prec, srcinfo, pg);
        E    = -U.grad'/(4*pi);
    end

    %% Normal electric field with correction by more exact integrals for near neighbors
    if(nargin == 6)
        correction = EC*c;                                           % Correction of plain FMM result
        En         = -c(targetIndex)/2 + correction(targetIndex) ... % This is the dominant (exact) matrix part and the "undo" terms for center-point FMM 
                        + sum(normals(targetIndex).*E, 2);
        % Paul:
        % I timed indexing first and then multiplication:
        %           EC(targetIndex, :)*c
        % versus normal multiplication and then indexing like is done
        % above.
        % The second, currently used option seems to be faster.
    else
        correction = EC*c;                  % Correction of plain FMM result
        En         = -c/2 + correction ...  % This is the dominant (exact) matrix part and the "undo" terms for center-point FMM
                        + sum(normals.*E, 2);
    end
end
