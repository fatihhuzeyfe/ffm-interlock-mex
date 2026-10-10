function fill_roads(ax, xi, yi, W, ylo, yhi, S, roads)
%FILL_ROADS  Material A below / B above the interface polyline (xi, yi) and the printed roads
%   (raster 0 deg: lines parallel to the interface direction, spacing 0.2 w).
if nargin < 8, roads = true; end
xi = xi(:)'; yi = yi(:)';
PA = [xi W 0; yi ylo ylo]'; PB = [xi W 0; yi yhi yhi]';
patch(ax, PA(:, 1), PA(:, 2), S.matA, 'EdgeColor', 'none');
patch(ax, PB(:, 1), PB(:, 2), S.matB, 'EdgeColor', 'none');
if roads
    xs = linspace(0, W, 600);
    for yr = ylo + 0.1:0.2:yhi
        for k = 1:2
            if k == 1, P = PA; c = S.roadA; else, P = PB; c = S.roadB; end
            in = inpolygon(xs, yr * ones(size(xs)), P(:, 1), P(:, 2));
            yy = yr * ones(size(xs)); yy(~in) = NaN;
            plot(ax, xs, yy, '-', 'Color', c, 'LineWidth', 0.35);
        end
    end
end
plot(ax, xi, yi, 'k-', 'LineWidth', 1.0);
end
