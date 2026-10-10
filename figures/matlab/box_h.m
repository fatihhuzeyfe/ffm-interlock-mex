function box_h(ax, y, v, col)
%BOX_H  Horizontal box plot of the values v at height y (no toolbox needed).
q = quant(v, [0.25 0.5 0.75]); iqr = q(3) - q(1);
lo = min(v(v >= q(1) - 1.5 * iqr)); hi = max(v(v <= q(3) + 1.5 * iqr));
h = 0.28;
patch(ax, [q(1) q(3) q(3) q(1)], [y - h y - h y + h y + h], col, 'EdgeColor', 'k', 'LineWidth', 0.5);
plot(ax, [q(2) q(2)], [y - h y + h], 'k-', 'LineWidth', 0.9);
plot(ax, [lo q(1)], [y y], 'k-', 'LineWidth', 0.5); plot(ax, [q(3) hi], [y y], 'k-', 'LineWidth', 0.5);
out = v(v < lo | v > hi);
plot(ax, out, y * ones(size(out)), '.', 'Color', [0.55 0.55 0.55], 'MarkerSize', 4);
end
