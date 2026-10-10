%FIG5_CZM  Cohesive-zone counterpart of the FFM predictions.
%   (a, b) nominal stress vs displacement / FFM onset displacement, PLA/PETG, k = 1, w = 5 mm,
%   dots = FFM onset; (c) CZM onset / FFM onset vs lch/w; (d) gain vs lch/w: FFM master curve,
%   FFM and CZM for the dovetail sweep; (e) Abaqus CZM model of the dovetail in tension:
%   deformed mesh and interface damage before onset, at onset and after debonding.
S = ffm_style();
C  = read_csv(fullfile(S.data, 'fig5a_czm_curves.csv'));
Z  = read_csv(fullfile(S.data, 'fig5b_czm_summary.csv'));
B  = read_csv(fullfile(S.data, 'fig2a_master_curve_binned.csv'));

fig = new_figure(S.W, 14.4);
ref = Z.k_sigma == 1 & Z.k_G == 1 & Z.w_mm == 5;
geos  = {'flat', 'sine_ar1', 'U_ar1', 'T_ar1', 'dovetail_ar0.5_a10', 'dovetail_ar1_a20', 'dovetail_ar2_a30'};
gname = {'flat', 'sinusoidal', 'U', 'T', 'dovetail 0.5/10\circ', 'dovetail 1/20\circ', 'dovetail 2/30\circ'};
gcol  = [0.45 0.45 0.45; S.green; S.orange; S.purple; S.sky; S.blue; 0 0.25 0.45];
loads = {'tension', 'shear'};

% ---------------------------------------------------------------- (a, b) curves
for L = 1:2
    x0 = 1.4 + (L - 1) * 6.1;
    ax = add_axes(fig, [x0, 9.4, 4.6, 3.9]);
    for g = 1:numel(geos)
        r = ref & sel(Z.geometry, geos{g}) & sel(Z.pair, 'PLA/PETG') & sel(Z.load, loads{L});
        if ~any(r), continue; end
        c = sel(C.name, Z.name{find(r, 1)});
        plot(ax, C.u_over_u_ffm(c), C.sigma(c), '-', 'Color', gcol(g, :), 'LineWidth', 0.9);
        plot(ax, 1, Z.sigma_ffm(find(r, 1)), 'o', 'MarkerSize', 3.6, 'MarkerFaceColor', gcol(g, :), 'MarkerEdgeColor', 'k', 'LineWidth', 0.4);
    end
    set(ax, 'XLim', [0 5], 'YLim', [0 60], 'XTick', 0:5);
    xlab(ax, '\itu\rm / \itu\rm_{FFM}'); ylab(ax, '\sigma_{nom}  [MPa]');
    if L == 1                                   % geometry key across the top of the figure
        ents = cell(numel(geos) + 1, 2);
        for g = 1:numel(geos), ents(g, :) = {{'-', 'Color', gcol(g, :), 'LineWidth', 1.2}, gname{g}}; end
        ents(end, :) = {{'o', 'MarkerFaceColor', [0.6 0.6 0.6], 'MarkerEdgeColor', 'k', 'MarkerSize', 3.6}, 'FFM onset'};
        key(ax, ents, 'Orientation', 'horizontal');
        lg = findobj(fig, 'Type', 'legend'); set(lg(1), 'Units', 'centimeters', 'Position', [0.8 13.85 17.4 0.45]);
        set(ax, 'Units', 'centimeters', 'Position', [x0, 9.4, 4.6, 3.9]);
    end
    panel_caption(ax, char('a' + L - 1), sprintf('CZM response, %s (PLA/PETG)', loads{L}));
end

% ---------------------------------------------------------------- (c) CZM / FFM vs lch/w
ax = add_axes(fig, [13.6 9.4 4.6 3.9]);
patch(ax, [0.08 10 10 0.08], [0.9 0.9 1.1 1.1], [0.93 0.93 0.93], 'EdgeColor', 'none');
plot(ax, [0.08 10], [1 1], 'k-', 'LineWidth', 0.6);
mk = {'o', 's'};
for g = 1:numel(geos)
    for L = 1:2
        r = sel(Z.geometry, geos{g}) & sel(Z.load, loads{L});
        if ~any(r), continue; end
        fc = gcol(g, :); if L == 2, fc = 'w'; end
        plot(ax, Z.lch_over_w(r) .* (1 + 0.04 * (g - 4)), Z.czm_over_ffm(r), mk{L}, 'MarkerSize', 3.4, ...
             'MarkerEdgeColor', gcol(g, :), 'MarkerFaceColor', fc, 'LineWidth', 0.7);
    end
end
set(ax, 'XScale', 'log', 'XLim', [0.08 10], 'YLim', [0.45 1.15], 'XTick', [0.1 1 10], 'XTickLabel', {'0.1', '1', '10'});
xlab(ax, '\itl\rm_{ch} / \itw\rm'); ylab(ax, '\sigma_c^{CZM} / \sigma_c^{FFM}');
key(ax, {{'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 3.4}, 'tension'; ...
         {'s', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 3.4}, 'shear'}, 'Location', 'southwest');
panel_caption(ax, 'c', 'Onset: CZM relative to FFM');

% ---------------------------------------------------------------- (d) gain: FFM vs CZM
ax = add_axes(fig, [1.4 2.4 6.4 4.4]);
b = sel(B.load, 'tension') & sel(B.pair, 'PLA/PETG');
[x, o] = sort(B.lch_over_w(b)); p10 = B.gain_p10(b); p90 = B.gain_p90(b); m = B.gain_median(b);
for i = 1:numel(x), plot(ax, [x(i) x(i)], [p10(o(i)) p90(o(i))], '-', 'Color', S.blue * 0.4 + 0.6, 'LineWidth', 0.6); end
plot(ax, x, m(o), '-', 'Color', S.blue, 'LineWidth', 1.0);
sw = sel(Z.geometry, 'dovetail_ar1_a20') & sel(Z.pair, 'PLA/PETG') & sel(Z.load, 'tension') & ~isnan(Z.gain_czm);
[xs, os] = sort(Z.lch_over_w(sw)); gf = Z.gain_ffm(sw); gc = Z.gain_czm(sw);
plot(ax, xs, gf(os), '--o', 'Color', 'k', 'MarkerSize', 4.0, 'MarkerFaceColor', 'w', 'LineWidth', 0.7);
plot(ax, xs, gc(os), '-x', 'Color', S.red, 'MarkerSize', 6, 'LineWidth', 1.1);
plot(ax, [0.05 30], [1 1], 'k-', 'LineWidth', 0.6);
set(ax, 'XScale', 'log', 'XLim', [0.08 20], 'YLim', [0.8 2.2], 'XTick', [0.1 1 10], 'XTickLabel', {'0.1', '1', '10'});
xlab(ax, '\itl\rm_{ch} / \itw\rm'); ylab(ax, 'gain vs. flat');
key(ax, {{'-', 'Color', S.blue}, 'FFM, all geometries (median, 10-90 %)'; ...
         {'--o', 'Color', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 4}, 'FFM, dovetail 1/20\circ'; ...
         {'-x', 'Color', S.red, 'MarkerSize', 6, 'LineWidth', 1.1}, 'CZM, dovetail 1/20\circ'}, 'Location', 'northwest');
panel_caption(ax, 'd', 'Gain in tension (PLA/PETG): FFM and CZM');

% ---------------------------------------------------------------- (e) Abaqus CZM model, damage
tf = {'t012', 't018', 't030'}; tl = {'before onset', 'onset', 'after debonding'};
win = [4.5 15.5 -2.2 6.6]; iw = 3.05; ih = iw * (win(4) - win(3)) / (win(2) - win(1));
cmapD = interp1([0 0.35 0.7 1], [1 0.95 0.55; 0.96 0.56 0.12; 0.62 0.10 0.32; 0.05 0.02 0.08], linspace(0, 1, 128));
for k = 1:3
    pos = [9.1 + (k - 1) * 3.2, 3.75, iw, ih];
    ax = add_axes(fig, pos);
    f = fullfile(S.data, ['odb_czm_dovetail_ar1_a20_tension_' tf{k} '.json']);
    if exist(f, 'file') ~= 2
        txt(ax, 0.5, 0.5, 'export missing', 'Units', 'normalized', 'HorizontalAlignment', 'center');
        set(ax, 'XTick', [], 'YTick', []); continue
    end
    M = load_odb_json(f);
    sc = 3;
    draw_mesh(ax, M, 'mesh', sc);
    set(findobj(ax, 'Type', 'patch'), 'EdgeColor', 'none');
    D = M.damage;
    if ~isempty(D)
        [~, ia] = min(sum((M.V - permute(D(:, 1:2), [3 2 1])) .^ 2, 2), [], 1);
        P = M.V(ia(:), :) + sc * M.U(ia(:), :);
        scatter(ax, P(:, 1), P(:, 2), 5, D(:, 3), 'filled');
        colormap(ax, cmapD); caxis(ax, [0 1]);
    end
    axis(ax, 'equal'); set(ax, 'XLim', win(1:2), 'YLim', win(3:4), 'XTick', [], 'YTick', []);
    set(ax, 'Units', 'centimeters', 'Position', pos);
    title(ax, {tl{k}, sprintf('t = %.2f', M.time)}, 'FontName', S.font, 'FontSize', 6, 'FontWeight', 'normal', 'Color', 'k');
    if k == 2
        cb = colorbar(ax, 'Orientation', 'horizontal');
        set(cb, 'Units', 'centimeters', 'Position', [11.6 2.95 4.6 0.16], 'FontName', S.font, 'FontSize', 5.5, 'Color', 'k');
        xlabel(cb, 'interface damage', 'FontWeight', 'bold', 'FontSize', 6);
        set(ax, 'Units', 'centimeters', 'Position', pos);
    end
end
xe = 9.1 + (3 * 3.2 - 0.15) / 2;
fig_text(fig, xe, 1.04, 'Abaqus CZM model, dovetail in tension (deformation \times3)', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
fig_text(fig, xe, 0.62, '(e)', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', S.fsPanel, 'FontWeight', 'bold');

save_figure(fig, 'Figure_5');
