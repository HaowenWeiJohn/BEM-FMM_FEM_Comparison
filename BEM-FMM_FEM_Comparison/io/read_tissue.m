function tissue = read_tissue(filename)
%   Reads file containing names of different tissues in headmodel
%
%   Names must be seperated only by single line breaks
%   and must be ordered according to the indices of the
%   corresponding mesh simplices
%
%   By Paul Lunkenheimer

    tissue = {};

    FID = fopen(filename, 'r');
    counter = 0;
    while true
        line = fgetl(FID);
        if ~ischar(line); break; end
        counter = counter + 1;
        tissue{counter} = line;
    end

end