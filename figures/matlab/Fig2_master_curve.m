%FIG2_MASTER_CURVE  Onset-load gain of the interlocked interface over the flat one versus the
%   interface length scale lch/w. (a) tension, (b) shear: median of every lch/w bin (all
%   geometries, interface qualities, sizes and rasters) per material pair, whiskers = 10-90 %
%   range; (c) median gain of the geometry families in the two asymptotic regimes.
S = ffm_style();
B = read_csv(fullfile(S.data, 'fig2a_master_curve_binned.csv'));
F = read_csv(fullfile(S.data, 'fig2c_family_regimes.csv'));

fig = new_figure(S.W, 8.4);
XL = [0.04 600]; YL = [0.8 2.6];
loads = {'tension', 'shear'};
desc = {'Tension: gain vs. interface length scale', 'Shear: gain vs. interface length scale'};
pos = {[1.45 2.5 5.5 5.4], [7.75 2.5 5.5 5.4]};
off = [0.94 1.0 1.06];                          % small x offsets so the whiskers do not overlap
for L = 1:2
    ax = add_axes(fig, pos{L});
    patch(ax, [XL(1) 0.5 0.5 XL(1)], YL([1 1 2 2]), [0.94 0.94 0.94], 'EdgeColor', 'none');
    patch(ax, [2 XL(2) XL(2) 2], YL([1 1 2 2]), [0.94 0.94 0.94], 'EdgeColor', 'none');
    plot(ax, XL, [1 1], 'k-', 'LineWidth', 0.6);
    for p = 1:3
        b = sel(B.load, loads{L}) & sel(B.pair, S.pairs{p});
        [x, o] = sort(B.lch_over_w(b)); m = B.gain_median(b); lo = B.gain_p10(b); hi = B.gain_p90(b);
        m = m(o); lo = lo(o); hi = hi(o); xo = x * off(p);
        for i = 1:numel(xo)
            plot(ax, [xo(i) xo(i)], [lo(i) hi(i)], '-', 'Color', S.pairCol(p, :) * 0.5 + 0.5, 'LineWidth', 0.6);
        end
        plot(ax, xo, m, '-o', 'Color', S.pairCol(p, :), 'LineWidth', 1.1, 'MarkerSize', 3.2, 'MarkerFaceColor', 'w');
    end
    plat = median(B.gain_median(sel(B.load, loads{L}) & B.lch_over_w > 5));
    plot(ax, [2 XL(2)], [plat plat], 'k--', 'LineWidth', 0.6);
    txt(ax, 450, plat - 0.04, sprintf('plateau %.2f', plat), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    txt(ax, 0.14, 2.56, {'strength-', 'controlled'}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    txt(ax, 1.0, 2.56, 'transition', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'Rotation', 90);
    txt(ax, 35, 2.56, 'toughness-controlled', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    set(ax, 'XScale', 'log', 'XLim', XL, 'YLim', YL, 'XTick', [0.1 1 10 100], 'XTickLabel', {'0.1', '1', '10', '100'}, ...
        'YTick', 0.8:0.2:2.6);
    xlab(ax, 'interface length scale  \itl\rm_{ch} / \itw\rm');
    if L == 1
        ylab(ax, 'onset-load gain  \sigma_c^{interlock} / \sigma_c^{flat}');
        key(ax, {{'-o', 'Color', S.pairCol(1, :), 'MarkerFaceColor', 'w', 'MarkerSize', 3.2}, 'PLA/PETG'; ...
                 {'-o', 'Color', S.pairCol(2, :), 'MarkerFaceColor', 'w', 'MarkerSize', 3.2}, 'PLA/ABS'; ...
                 {'-o', 'Color', S.pairCol(3, :), 'MarkerFaceColor', 'w', 'MarkerSize', 3.2}, 'PETG/ABS'; ...
                 {'-', 'Color', [0.6 0.6 0.6]}, '10-90 % range'}, 'Location', 'west');
    end
    panel_caption(ax, char('a' + L - 1), desc{L});
end

% ---------------------------------------------------------------- (c) families
ax = add_axes(fig, [14.85 2.5 3.9 5.4]);
fams = {'sine', 'U', 'dovetail', 'T'};
mk = {'o', 's'};
for f = 1:4
    for L = 1:2
        y = f + (L - 1.5) * 0.32;
        r = sel(F.shape, fams{f}) & sel(F.load, loads{L});
        g1 = F.gain_median(r & sel(F.regime, 'strength'));
        g2 = F.gain_median(r & sel(F.regime, 'toughness'));
        c = S.shapeCol.(fams{f});
        plot(ax, [g1 g2], [y y], '-', 'Color', c * 0.5 + 0.5, 'LineWidth', 1.0);
        plot(ax, g1, y, mk{L}, 'MarkerSize', 4, 'MarkerEdgeColor', c, 'MarkerFaceColor', 'w', 'LineWidth', 0.8);
        plot(ax, g2, y, mk{L}, 'MarkerSize', 4, 'MarkerEdgeColor', c, 'MarkerFaceColor', c, 'LineWidth', 0.8);
    end
end
plot(ax, [1 1], [0.3 6.2], 'k-', 'LineWidth', 0.6);
set(ax, 'YDir', 'reverse', 'YLim', [0.3 6.2], 'YTick', 1:4, 'YTickLabel', fams, 'XLim', [0.85 2.2], ...
        'XTick', [1 1.5 2], 'TickLength', [0.015 0]);
xlab(ax, 'median gain');
key(ax, {{'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 4}, 'strength-controlled'; ...
         {'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 4}, 'toughness-controlled'; ...
         {'o', 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerSize', 4}, 'tension'; ...
         {'s', 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerSize', 4}, 'shear'}, 'Location', 'southeast');
panel_caption(ax, 'c', 'Geometry families in the two regimes');

save_figure(fig, 'Figure_2');
