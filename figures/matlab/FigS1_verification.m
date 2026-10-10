%FIGS1_VERIFICATION  Verification. (a) Griffith benchmark, (b) FFM strength vs exact solution,
%   (c) onset stress vs element size, (d) change production -> fine mesh per setup,
%   (e, f) coarse and fine Abaqus meshes at the dovetail corners, (g) close-up of the fine mesh
%   at the re-entrant corner, (h) stress sigma_22 of the uncracked model under tension.
S = ffm_style();
BG = read_csv(fullfile(S.data, 'figS1a_benchmark_G.csv'));
BS = read_csv(fullfile(S.data, 'figS1b_benchmark_strength.csv'));
MC = read_csv(fullfile(S.data, 'figS1c_mesh_convergence.csv'));

fig = new_figure(S.W, 13.6);
y1 = 8.5; h1 = 3.8;

% ---------------------------------------------------------------- (a)
ax = add_axes(fig, [1.3 y1 3.4 h1]);
patch(ax, [0.007 3 3 0.007], [0.98 0.98 1.02 1.02], [0.93 0.93 0.93], 'EdgeColor', 'none');
plot(ax, [0.007 3], [1 1], 'k-', 'LineWidth', 0.6);
plot(ax, BG.da, BG.G_ratio, '-o', 'Color', S.blue, 'MarkerFaceColor', 'w', 'MarkerSize', 3);
plot(ax, BG.da, BG.closure_over_energy, '-s', 'Color', S.red, 'MarkerFaceColor', 'w', 'MarkerSize', 2.8);
set(ax, 'XScale', 'log', 'XLim', [0.007 3], 'YLim', [0.95 1.05], 'XTick', [0.01 0.1 1], 'XTickLabel', {'0.01', '0.1', '1'});
xlab(ax, '\Delta\ita\rm  [mm]'); ylab(ax, 'FE / exact');
key(ax, {{'-o', 'Color', S.blue, 'MarkerFaceColor', 'w', 'MarkerSize', 3}, 'energy release rate'; ...
         {'-s', 'Color', S.red, 'MarkerFaceColor', 'w', 'MarkerSize', 2.8}, 'closure / energy'}, 'Location', 'southwest');
panel_caption(ax, 'a', 'Griffith benchmark');

% ---------------------------------------------------------------- (b)
ax = add_axes(fig, [5.9 y1 3.1 h1]);
plot(ax, BS.a_over_lch, BS.sigma_exact, 'k-', 'LineWidth', 0.8);
plot(ax, BS.a_over_lch, BS.sigma_ffm_fe, 'o', 'MarkerEdgeColor', S.blue, 'MarkerFaceColor', 'w', 'MarkerSize', 4, 'LineWidth', 0.9);
for i = 1:numel(BS.a_over_lch)
    txt(ax, BS.a_over_lch(i) + 0.08, BS.sigma_ffm_fe(i) + 1.6, sprintf('%+.2f %%', 100 * (BS.sigma_ffm_fe(i) / BS.sigma_exact(i) - 1)), 'FontSize', 5.4);
end
set(ax, 'XLim', [0 1.7], 'YLim', [20 52]);
xlab(ax, '\ita\rm / \itl\rm_{ch}'); ylab(ax, 'onset stress  [MPa]');
key(ax, {{'-', 'Color', 'k'}, 'exact FFM'; {'o', 'MarkerEdgeColor', S.blue, 'MarkerFaceColor', 'w', 'MarkerSize', 4}, 'FE + FFM'}, 'Location', 'northeast');
panel_caption(ax, 'b', 'FFM strength vs. exact');

% ---------------------------------------------------------------- (c)
ax = add_axes(fig, [10.2 y1 3.1 h1]);
hh = [0.1 0.05 0.025];
setups = unique(MC.setup, 'stable');
patch(ax, [0.046 0.054 0.054 0.046], [0.99 0.99 1.05 1.05], [0.93 0.93 0.93], 'EdgeColor', 'none');
for i = 1:numel(setups)
    r = sel(MC.setup, setups{i}); ex = MC.sigma_extrapolated(r);
    y = [median(MC.sigma_coarse(r) ./ ex) median(MC.sigma_medium(r) ./ ex) median(MC.sigma_fine(r) ./ ex)];
    c = [0.55 0.55 0.55]; if ~isempty(strfind(setups{i}, 'sine')) && ~isempty(strfind(setups{i}, 'shear')), c = S.red; end
    plot(ax, hh, y, '-o', 'Color', c, 'MarkerFaceColor', 'w', 'MarkerSize', 2.8, 'LineWidth', 0.8);
end
plot(ax, [0.02 0.12], [1 1], 'k-', 'LineWidth', 0.6);
set(ax, 'XScale', 'log', 'XDir', 'reverse', 'XLim', [0.022 0.11], 'XTick', [0.025 0.05 0.1], ...
        'XTickLabel', {'0.025', '0.05', '0.1'}, 'XMinorTick', 'off', 'YLim', [0.995 1.04]);
xlab(ax, 'element size  \ith\rm / \itw\rm'); ylab(ax, '\sigma_c / \sigma_c^{extrapolated}');
key(ax, {{'-o', 'Color', [0.55 0.55 0.55], 'MarkerFaceColor', 'w', 'MarkerSize', 2.8}, 'setup (median)'; ...
         {'-o', 'Color', S.red, 'MarkerFaceColor', 'w', 'MarkerSize', 2.8}, 'with contact'}, 'Location', 'northwest');
panel_caption(ax, 'c', 'Mesh convergence');

% ---------------------------------------------------------------- (d)
ax = add_axes(fig, [16.0 y1 2.7 h1]);
ord = [setups(cellfun('isempty', strfind(setups, 'sine'))); setups(~cellfun('isempty', strfind(setups, 'sine')))];
for i = 1:numel(ord)
    r = sel(MC.setup, ord{i}); v = max(100 * abs(MC.change_medium_fine(r)), 1e-3);
    c = [0.82 0.82 0.82]; if ~isempty(strfind(ord{i}, 'sine')) && ~isempty(strfind(ord{i}, 'shear')), c = S.red; end
    box_h(ax, i, v, c);
    p = strsplit(ord{i}, '__');
    lab = sprintf('%s, %s, %s', strrep(strrep(p{1}, 'dovetail_', 'DT '), '_', ' '), strrep(p{2}, '-', '/'), upper(p{3}(1)));
    if strcmp(p{4}, 'r90'), lab = [lab ' r90']; end
    text(ax, 1.6e-3, i, [lab '  '], 'HorizontalAlignment', 'right', 'FontName', S.font, 'FontSize', 4.9, 'Color', 'k', 'Interpreter', 'none');
end
set(ax, 'XScale', 'log', 'XLim', [2e-3 10], 'XTick', [0.01 0.1 1 10], 'XTickLabel', {'0.01', '0.1', '1', '10'}, ...
        'YLim', [0.4 numel(ord) + 0.6], 'YTick', [], 'YDir', 'reverse');
xlab(ax, '|\Delta\sigma_c|  [%]');
title(ax, sprintf('median %.2f %%', 100 * median(abs(MC.change_medium_fine))), 'FontName', S.font, 'FontSize', 6, 'FontWeight', 'normal', 'Color', 'k');
panel_caption(ax, 'd', 'Production vs. fine mesh');

% ---------------------------------------------------------------- (e-h) Abaqus meshes and field
fc = fullfile(S.data, 'odb_mesh_coarse.json'); ff = fullfile(S.data, 'odb_mesh_fine.json');
win = [0.8 2.2 -0.35 1.25]; iw = 3.4; ih = iw * (win(4) - win(3)) / (win(2) - win(1));
fig_text(fig, 9.5, 5.75, 'Abaqus models: dovetail d/w = 1, \alpha = 20\circ, PLA/PETG', 'HorizontalAlignment', 'center', 'FontSize', S.fsSmall);
if exist(fc, 'file') == 2 && exist(ff, 'file') == 2
    Mc = load_odb_json(fc); Mf = load_odb_json(ff);
    yc = mean(reshape(Mf.V(Mf.F, 2), [], 3), 2);
    s22 = Mf.S(:, 2); Mf.S(:, 2) = s22 / median(s22(abs(yc) > 3 & ~isnan(s22)));
    P = {Mc, 'mesh', win, 'Coarse mesh (h/w = 0.1)', sprintf('%d elements', size(Mc.F, 1)); ...
         Mf, 'mesh', win, 'Fine mesh (h/w = 0.025)', sprintf('%d elements', size(Mf.F, 1)); ...
         Mf, 'mesh', [1.38 1.62 -0.1 0.1653], 'Fine mesh, corner close-up', ''; ...
         Mf, 'S22', win, 'Stress \sigma_{22}, tension', ''};
    for k = 1:4
        pos = [1.0 + (k - 1) * 4.35, 1.4, iw, ih];
        ax = add_axes(fig, pos);
        if strcmp(P{k, 2}, 'mesh')
            draw_mesh(ax, P{k, 1}, 'mesh');
            if k == 3, set(findobj(ax, 'Type', 'patch'), 'LineWidth', 0.3); end
        else
            draw_mesh(ax, P{k, 1}, 'S22', 0, [0.4 1.6]); colormap(ax, parula(256));
            cb = colorbar(ax, 'Units', 'centimeters', 'Position', [pos(1) + iw + 0.12, pos(2), 0.16, ih]);
            set(cb, 'FontName', S.font, 'FontSize', 5.5, 'Color', 'k'); ylabel(cb, '\sigma_{22} / \sigma_{nom}', 'FontWeight', 'bold');
        end
        w4 = P{k, 3};
        axis(ax, 'equal'); set(ax, 'XLim', w4(1:2), 'YLim', w4(3:4), 'XTick', [], 'YTick', []);
        set(ax, 'Units', 'centimeters', 'Position', pos);
        if ~isempty(P{k, 5})
            text(ax, 0.03, 0.97, P{k, 5}, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontName', S.font, ...
                 'FontSize', 5.6, 'Color', 'k', 'BackgroundColor', 'w', 'Margin', 1);
        end
        panel_caption(ax, char('e' + k - 1), P{k, 4}, 0.3);
    end
else
    fig_text(fig, 9.5, 3.5, 'Abaqus mesh exports missing: run  python scripts/export_figure_meshes.py', 'HorizontalAlignment', 'center');
end

save_figure(fig, 'Figure_S1_verification');
