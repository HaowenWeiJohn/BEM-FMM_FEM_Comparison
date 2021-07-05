function [P, C, n_of_charges] = read_charges(filename)
%   Reads ascii file containing charge distribution
%
%   Charges must be seperated only by single line breaks
%   Each line contains 4 entries seperated by single tabs
%   First 3 entries are x, y and z coordinate of charge location
%   Last entry is charge
%
%   By Paul Lunkenheimer

    FID = fopen(filename, 'r');
    n_of_charges = 0;
    while true
        line = fgetl(FID);
        if ~ischar(line); break; end
        n_of_charges = n_of_charges + 1;
        data = split(line, "	");
        P(n_of_charges, 1) = str2double(data{1});
        P(n_of_charges, 2) = str2double(data{2});
        P(n_of_charges, 3) = str2double(data{3});
        C(n_of_charges, 1) = str2double(data{4});
    end
    
end