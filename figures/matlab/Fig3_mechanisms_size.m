%FIG3_MECHANISMS_SIZE  Where the crack starts and which geometry gains most.
%   (a, b) share of onset mechanisms vs lch/w (interlocked geometries), tension / shear;
%   (c) size effect; (d) raster direction; (e) gain of every geometry at the reference interface
%   (k = 1, w = 5 mm, raster 0 deg): one column per geometry (icon and d/w, alpha below), dots
%   = material pairs, grey bar = range over the pairs; tension above, shear below.
S = ffm_style();
A  = read_csv(fullfile(S.data, 'fig2b_all_points.csv'));
Z  = read_csv(fullfile(S.data, 'fig3b_size_effect.csv'));
T  = read_csv(fullfile(S.data, 'fig3c_gain_table.csv'));
R  = read_csv(fullfile(S.data, 'fig3d_raster_ratio.csv'));
GE = read_csv(fullfile(S.data, 'fig_all_geometries.csv'));

fig = new_figure(S.W, 15.9);
loads = {'tension', 'shear'};

% ---------------------------------------------------------------- (a, b) mechanisms
il = (sel(A.shape, 'T') | sel(A.shape, 'U') | sel(A.shape, 'dovetail')) & A.GIIc_over_GIc == 2;
edges = logspace(log10(0.04), log10(600), 22); mid = sqrt(edges(1:end - 1) .* edges(2:end));
colM = [0.72 0.80 0.88; S.blue; S.red];
for L = 1:2
    ax = add_axes(fig, [1.3 + (L - 1) * 4.6, 11.4, 3.5, 3.6]);
    m = il & sel(A.load, loads{L});
    Y = zeros(numel(mid), 3); n = zeros(numel(mid), 1);
    for b = 1:numel(mid)
        in = m & A.lch_over_w >= edges(b) & A.lch_over_w < edges(b + 1);
        n(b) = sum(in);
        if n(b)
            Y(b, :) = [sum(in & sel(A.mechanism, 'interface_base')), sum(in & sel(A.mechanism, 'interface_feature')), ...
                       sum(in & (sel(A.mechanism, 'bulk_A') | sel(A.mechanism, 'bulk_B')))] / n(b);
        end
    end
    ok = n > 20;
    h = area(ax, mid(ok), Y(ok, :), 'LineStyle', 'none');
    for k = 1:3, set(h(k), 'FaceColor', colM(k, :)); end
    set(ax, 'XScale', 'log', 'XLim', [min(mid(ok)) max(mid(ok))], 'YLim', [0 1], 'XTick', [0.1 1 10 100], ...
        'XTickLabel', {'0.1', '1', '10', '100'});
    ylab(ax, 'share of onsets'); xlab(ax, '\itl\rm_{ch} / \itw\rm');
    if L == 1
        lg = legend(ax, h, {'base interface', 'on the feature', 'bulk fracture'}, 'Orientation', 'horizontal', ...
                    'FontName', S.font, 'FontSize', S.fsSmall, 'TextColor', 'k');
        set(lg, 'Box', 'off', 'Units', 'centimeters', 'Position', [1.2 15.15 8.0 0.4]);
    end
    panel_caption(ax, char('a' + L - 1), sprintf('Onset location, %s', loads{L}));
end

% ---------------------------------------------------------------- (c) size effect
ax = add_axes(fig, [10.75 11.4 3.5 3.6]);
mk = {'o', 's'}; ls = {'-', '--'};
for L = 1:2
    for p = 1:3
        r = sel(Z.load, loads{L}) & sel(Z.pair, S.pairs{p});
        w = unique(Z.width_mm(r)); y = zeros(size(w));
        for i = 1:numel(w), y(i) = median(Z.gain_median(r & Z.width_mm == w(i))); end
        plot(ax, w, y, [ls{L} mk{L}], 'Color', S.pairCol(p, :), 'MarkerFaceColor', 'w', 'MarkerSize', 3.2, 'LineWidth', 1.0);
    end
end
plot(ax, [1 11], [1 1], 'k-', 'LineWidth', 0.6);
set(ax, 'XLim', [1.4 10.6], 'XTick', [2 5 10], 'YLim', [0.95 2.1]);
xlab(ax, 'feature width  w  [mm]'); ylab(ax, 'median gain');
key(ax, {{'-o', 'Color', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 3}, 'tension'; ...
         {'--s', 'Color', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 3}, 'shear'}, 'Location', 'northeast');
panel_caption(ax, 'c', 'Size effect (k = 1)');

% ---------------------------------------------------------------- (d) raster direction
ax = add_axes(fig, [15.3 11.4 3.5 3.6]);
x = 100 * (R.r90_over_r0 - 1);
e = -8:0.25:8; c = histc(x, e); c = c(1:end - 1);
for i = 1:numel(c)
    if c(i) > 0
        patch(ax, [e(i) e(i + 1) e(i + 1) e(i)], [3 3 c(i) c(i)], [0.45 0.45 0.45], 'EdgeColor', 'none');
    end
end
q = quant(x, [0.05 0.5 0.95]);
plot(ax, [0 0], [3 3e5], 'k-', 'LineWidth', 0.6);
set(ax, 'YScale', 'log', 'YLim', [3 3e5], 'XLim', [-8 8], 'YTick', [10 1e2 1e3 1e4]);
xlab(ax, 'change of onset load  [%]'); ylab(ax, 'configurations');
txt(ax, -7.6, 1.6e5, {sprintf('median %+.1f %%', q(2)), sprintf('90 %%: %+.1f ... %+.1f %%', q(1), q(3))}, 'VerticalAlignment', 'top');
panel_caption(ax, 'd', 'Raster 90\circ vs. 0\circ');

% ---------------------------------------------------------------- (e) ranking of all geometries
order = {'sine_ar0.5', 'sine_ar1', 'sine_ar1.5', 'sine_ar2', 'U_ar0.5', 'U_ar1', 'U_ar1.5', 'U_ar2', ...
    'dovetail_ar0.5_a10', 'dovetail_ar0.5_a20', 'dovetail_ar0.5_a30', 'dovetail_ar1_a10', 'dovetail_ar1_a20', ...
    'dovetail_ar1_a30', 'dovetail_ar1.5_a10', 'dovetail_ar1.5_a20', 'dovetail_ar1.5_a30', 'dovetail_ar2_a10', ...
    'dovetail_ar2_a20', 'dovetail_ar2_a30', 'T_ar0.5', 'T_ar1', 'T_ar1.5', 'T_ar2'};
fam = [repmat({'sine'}, 1, 4), repmat({'U'}, 1, 4), repmat({'dovetail'}, 1, 12), repmat({'T'}, 1, 4)];
xs = zeros(1, numel(order)); xp = 0;
for i = 1:numel(order)
    if i > 1 && ~strcmp(fam{i}, fam{i - 1}), xp = xp + 0.7; end
    xs(i) = xp; xp = xp + 1;
end
XL = [-0.7, xs(end) + 0.7]; X0 = 1.7; WID = 16.9;
xcm = @(xd) X0 + (xd - XL(1)) / (XL(2) - XL(1)) * WID;
for L = 1:2
    ax = add_axes(fig, [X0, 6.15 - (L - 1) * 2.95, WID, 2.7]);
    patch(ax, [XL(1) XL(2) XL(2) XL(1)], [0.8 0.8 1 1], [0.96 0.89 0.89], 'EdgeColor', 'none');
    plot(ax, XL, [1 1], 'k-', 'LineWidth', 0.6);
    for i = 1:numel(order)
        v = zeros(1, 3);
        for p = 1:3
            r = sel(T.geometry, order{i}) & sel(T.pair, S.pairs{p}) & sel(T.load, loads{L});
            v(p) = T.gain(find(r, 1));
        end
        plot(ax, [xs(i) xs(i)], [min(v) max(v)], '-', 'Color', [0.73 0.73 0.73], 'LineWidth', 0.9);
        for p = 1:3
            plot(ax, xs(i), v(p), 'o', 'MarkerSize', 3.6, 'MarkerFaceColor', S.pairCol(p, :), 'MarkerEdgeColor', 'w', 'LineWidth', 0.4);
        end
    end
    set(ax, 'XLim', XL, 'YLim', [0.8 2.5], 'XTick', [], 'YTick', [1 1.5 2 2.5]);
    ylab(ax, 'gain vs. flat');
    txt(ax, XL(1) + 0.2, 2.42, loads{L}, 'FontWeight', 'bold', 'FontSize', S.fs, 'VerticalAlignment', 'top');
    txt(ax, XL(1) + 0.2, 0.83, 'below flat (< 1)', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
    if L == 1
        key(ax, {{'o', 'MarkerFaceColor', S.pairCol(1, :), 'MarkerEdgeColor', 'w', 'MarkerSize', 4}, 'PLA/PETG'; ...
                 {'o', 'MarkerFaceColor', S.pairCol(2, :), 'MarkerEdgeColor', 'w', 'MarkerSize', 4}, 'PLA/ABS'; ...
                 {'o', 'MarkerFaceColor', S.pairCol(3, :), 'MarkerEdgeColor', 'w', 'MarkerSize', 4}, 'PETG/ABS'; ...
                 {'-', 'Color', [0.73 0.73 0.73], 'LineWidth', 0.9}, 'range over pairs'}, 'Orientation', 'horizontal');
        lg = findobj(fig, 'Type', 'legend'); set(lg(1), 'Units', 'centimeters', 'Position', [11.0 8.35 7.4 0.4]);
    end
end
% geometry icons with d/w (and alpha) under each column
for i = 1:numel(order)
    g = sel(GE.geometry, order{i}); W = GE.W(find(g, 1)); ar = GE.ar(find(g, 1)); al = GE.alpha_deg(find(g, 1));
    ia = add_axes(fig, [xcm(xs(i)) - 0.3, 2.45, 0.6, 0.55]);
    fill_roads(ia, GE.x(g), GE.y(g), W, -1.4, 2.3, S, false);
    set(ia, 'XLim', [0 W], 'YLim', [-1.4 2.3], 'Visible', 'off');
    if strcmp(fam{i}, 'dovetail'), s = {sprintf('%g', ar), sprintf('%g\\circ', al)}; else, s = sprintf('%g', ar); end
    fig_text(fig, xcm(xs(i)), 2.38, s, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 5.2);
end
fig_text(fig, X0 - 0.15, 2.38, {'d/w', '\alpha'}, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 5.2);
famName = {'sinusoidal', 'U', 'dovetail', 'T'}; famIdx = [1 4; 5 8; 9 20; 21 24];
for f = 1:4
    x1 = xcm(xs(famIdx(f, 1)) - 0.4); x2 = xcm(xs(famIdx(f, 2)) + 0.4);
    fig_line(fig, [x1 x2], [1.55 1.55]);
    fig_text(fig, (x1 + x2) / 2, 1.48, famName{f}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
             'FontWeight', 'bold', 'FontSize', 6.5);
end
fig_text(fig, X0 + WID / 2, 0.98, 'Gain of every geometry at the reference interface (k = 1, w = 5 mm, raster 0\circ)', ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
fig_text(fig, X0 + WID / 2, 0.55, '(e)', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
         'FontSize', S.fsPanel, 'FontWeight', 'bold');

save_figure(fig, 'Figure_3');
