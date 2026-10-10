function T = read_csv(file)
%READ_CSV  Read a CSV file with a header row into a struct of columns.
%   Numeric columns become double vectors, text columns cell arrays of char.
%   Works in MATLAB and GNU Octave (no toolbox needed); quoted fields with commas are handled.
txt = fileread(file);
txt = strrep(txt, sprintf('\r'), '');
lines = regexp(txt, '\n', 'split');
lines = lines(~cellfun('isempty', lines));
hdr = splitline(lines{1});
n = numel(lines) - 1;
C = cell(n, numel(hdr));
for i = 1:n
    f = splitline(lines{i + 1});
    C(i, 1:numel(f)) = f;
end
T = struct();
for j = 1:numel(hdr)
    name = regexprep(hdr{j}, '[^A-Za-z0-9_]', '_');
    if isempty(name) || ~isletter(name(1)), name = ['c_' name]; end
    col = C(:, j);
    col(cellfun('isempty', col)) = {''};
    v = str2double(col);
    if all(~isnan(v) | strcmp(col, '') | strcmpi(col, 'nan'))
        T.(name) = v;
    else
        T.(name) = col;
    end
end
end

function f = splitline(s)
if any(s == '"')
    f = regexp(s, ',(?=(?:[^"]*"[^"]*")*[^"]*$)', 'split');
    f = regexprep(f, '^"|"$', '');
else
    f = strsplit(s, ',');
end
end
