function x = l21norm(E)
%   L21 norm for the vector field (not complex)
%   SNM 2021

    x = sum(sqrt(dot(E, E, 2)));
end
