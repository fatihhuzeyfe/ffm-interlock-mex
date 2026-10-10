function arrowhead(ax, p, d, len, col)
%ARROWHEAD  Filled arrow head at point p pointing along direction d (data units, equal axes).
d = d / norm(d); n = [-d(2) d(1)];
w = 0.45 * len;
P = [p; p - len * d + w * n; p - len * d - w * n];
patch(ax, P(:, 1), P(:, 2), col, 'EdgeColor', 'none');
end
