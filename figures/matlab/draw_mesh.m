function draw_mesh(ax, M, mode, scale, clim_)
%DRAW_MESH  Abaqus mesh (mode 'mesh': materials in colour, element edges) or element-average
%   S22 contour (mode 'S22', flat-shaded triangles, vector-friendly). scale = deformation scale.
S = ffm_style();
if nargin < 4, scale = 0; end
V = M.V + scale * M.U;
switch mode
    case 'mesh'
        C = repmat(S.matA, size(M.F, 1), 1); C(M.mat == 1, :) = repmat(S.matB, sum(M.mat == 1), 1);
        patch(ax, 'Faces', M.F, 'Vertices', V, 'FaceVertexCData', C, 'FaceColor', 'flat', ...
              'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 0.15);
    case 'S22'
        patch(ax, 'Faces', M.F, 'Vertices', V, 'FaceVertexCData', M.S(:, 2), 'FaceColor', 'flat', ...
              'EdgeColor', 'none');
        if nargin >= 5, caxis(ax, clim_); end
end
end
