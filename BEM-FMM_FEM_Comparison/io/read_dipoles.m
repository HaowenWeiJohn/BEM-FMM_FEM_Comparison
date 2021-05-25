function [P, M, n_of_dip] = read_dipoles(filename)
%   Reads file containing dipoles
%
%   Dipoles must be seperated only by single line breaks
%   Each line contains 6 entries seperated by single tabs
%   First 3 entries are x, y and z coordinate of dipole location
%   Last 3 entries are x, y and z coordinate of dipole moment
%
%   By Paul Lunkenheimer

    FID = fopen(filename, 'r');
    n_of_dip = 0;
    while true
        line = fgetl(FID);
        if ~ischar(line); break; end
        n_of_dip = n_of_dip + 1;
        data = split(line, "	");
        P(n_of_dip, 1) = str2double(data{1});
        P(n_of_dip, 2) = str2double(data{2});
        P(n_of_dip, 3) = str2double(data{3});
        M(n_of_dip, 1) = str2double(data{4});
        M(n_of_dip, 2) = str2double(data{5});
        M(n_of_dip, 3) = str2double(data{6});
    end
    
end