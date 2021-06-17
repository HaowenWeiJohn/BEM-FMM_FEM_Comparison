function [P] = read_electrodes(filename)
%   Reads file containing electrodes
%
%   Electrodes must be seperated only by single line breaks
%   Each line contains 3 entries seperated by single spaces
%   Entries are x, y and z coordinate of electrode location
%
%   By Paul Lunkenheimer

    FID = fopen(filename, 'r');
    n_of_el = 0;
    while true
        line = fgetl(FID);
        if ~ischar(line); break; end
        n_of_el = n_of_el + 1;
        data = split(line, " ");
        P(n_of_el, 1) = str2double(data{1});
        P(n_of_el, 2) = str2double(data{2});
        P(n_of_el, 3) = str2double(data{3});
    end
    
end