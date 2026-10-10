function M = load_odb_json(file)
%LOAD_ODB_JSON  Mesh and fields exported by abaqus/export_figure_odb.py.
%   M.V (n x 2) node coordinates, M.U (n x 2) displacements, M.F (m x 3) corner-node indices of
%   the triangles, M.mat (m x 1) 0 = material A / 1 = B, M.S (m x 3) element-average stresses
%   S11, S22, S12 (global axes), M.damage (k x 3) [x y CSDMG] on the cohesive interface (CZM).
J = jsondecode(fileread(file));
N = J.nodes; if iscell(N), N = cell2mat(cellfun(@(r) r(:)', N, 'UniformOutput', false)); end
E = J.elements; if iscell(E), E = cell2mat(cellfun(@(r) r(:)', E, 'UniformOutput', false)); end
idx = zeros(max(N(:, 1)), 1); idx(N(:, 1)) = 1:size(N, 1);
M.V = N(:, 2:3); M.U = N(:, 4:5);
M.F = idx(E(:, 3:5)); M.mat = E(:, 2);
M.time = J.frame_time;
M.S = nan(size(E, 1), 3);
if isfield(J, 'stress_nodal')
    SN = J.stress_nodal; if iscell(SN), SN = cell2mat(cellfun(@(r) r(:)', SN, 'UniformOutput', false)); end
    for mi = 0:1
        r = SN(:, 1) == mi;
        nodal = nan(max(N(:, 1)), 3);
        nodal(SN(r, 2), :) = SN(r, 3:5);
        e = M.mat == mi;
        lab = E(e, 3:5);                                % corner-node labels of these elements
        for c = 1:3
            M.S(e, c) = mean(reshape(nodal(lab, c), [], 3), 2, 'omitnan');
        end
    end
end
M.damage = zeros(0, 3);
if isfield(J, 'damage') && ~isempty(J.damage)
    D = J.damage; if iscell(D), D = cell2mat(cellfun(@(r) r(:)', D, 'UniformOutput', false)); end
    M.damage = D;
end
end
