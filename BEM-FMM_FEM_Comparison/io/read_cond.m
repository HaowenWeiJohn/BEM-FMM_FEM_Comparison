function cond = read_cond(filename)
%   Reads file containing conductivites of different tissues in mesh of
%   head
%
%   Conductivities must be seperated only by single line breaks
%   and must be ordered according to the indices of the
%   corresponding mesh simplices
%
%   By Paul Lunkenheimer

    cond = double.empty;

    FID = fopen(filename, 'r');
    while true
        line = fgetl(FID);
        if ~ischar(line); break; end
        cond = [cond str2double(line)];
    end

end