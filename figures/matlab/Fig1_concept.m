%FIG1_CONCEPT  (a) interface families (periodic unit cells, printed roads), (b) the 24 geometries
%   of the study: aspect ratio d/w and dovetail flank angle alpha, (c) candidate crack paths of a
%   dovetail, (d) coupled stress-energy criterion of the governing path.
S  = ffm_style();
G  = read_csv(fullfile(S.data, 'fig1a_interfaces.csv'));
GE = read_csv(fullfile(S.data, 'fig_all_geometries.csv'));
P  = read_csv(fullfile(S.data, 'fig1b_candidate_paths.csv'));
I  = read_csv(fullfile(S.data, 'fig1b_interface.csv'));
L  = read_csv(fullfile(S.data, 'fig1b_crack_lengths.csv'));
CV = read_csv(fullfile(S.data, 'fig1c_ffm_curves.csv'));
ON = read_csv(fullfile(S.data, 'fig1c_ffm_onset.csv'));

fig = new_figure(S.W, 19.6);
YL = 1.7;

% ---------------------------------------------------------------- (a) families
names  = {'flat', 'sine_ar1', 'U_ar1', 'T_ar1', 'dovetail_ar1_a20'};
titles = {'flat', 'sinusoidal', 'U', 'T', 'dovetail'};
cw = 2.9; gap = 0.55; x0 = 0.9;
axa = cell(1, 5);
for i = 1:5
    ax = add_axes(fig, [x0 + (i - 1) * (cw + gap), 16.3, cw, cw * 2 * YL / 4]);
    axa{i} = ax;
    r = sel(G.geometry, names{i}); W = G.W(find(r, 1));
    fill_roads(ax, G.x(r), G.y(r), W, -YL, YL, S);
    plot(ax, [0 0], [-YL YL], 'k--', 'LineWidth', 0.5); plot(ax, [W W], [-YL YL], 'k--', 'LineWidth', 0.5);
    axis(ax, 'equal'); set(ax, 'XLim', [0 W], 'YLim', [-YL YL], 'Visible', 'off');
    txt(ax, W / 2, -YL - 0.2, titles{i}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', S.fs);
    if i == 1
        txt(ax, 0.2, -1.25, 'A', 'FontWeight', 'bold', 'FontSize', S.fs);
        txt(ax, 0.2, 1.25, 'B', 'FontWeight', 'bold', 'FontSize', S.fs);
    end
    if i == 5
        plot(ax, [1.5 2.5], [-0.35 -0.35], 'k-', 'LineWidth', 0.5);
        arrowhead(ax, [1.5 -0.35], [-1 0], 0.15, 'k'); arrowhead(ax, [2.5 -0.35], [1 0], 0.15, 'k');
        txt(ax, 2.0, -0.62, 'w', 'HorizontalAlignment', 'center', 'FontSize', S.fs, 'FontAngle', 'italic');
        plot(ax, [3.2 3.2], [0 1], 'k-', 'LineWidth', 0.5);
        arrowhead(ax, [3.2 0], [0 -1], 0.15, 'k'); arrowhead(ax, [3.2 1], [0 1], 0.15, 'k');
        txt(ax, 3.35, 0.5, 'd', 'FontSize', S.fs, 'FontAngle', 'italic');
        t = linspace(70, 90, 20) * pi / 180;
        plot(ax, 2.5 + 0.65 * cos(t), 0.65 * sin(t), 'k-', 'LineWidth', 0.5);
        plot(ax, [2.5 2.5], [0 0.8], 'k:', 'LineWidth', 0.5);
        txt(ax, 2.62, 0.95, '\alpha', 'FontSize', S.fs);
    end
end
% loading
fig_text(fig, 18.55, 19.35, 'loading', 'HorizontalAlignment', 'center', 'FontSize', S.fsSmall);
annotation(fig, 'arrow', [0.955 0.955], [0.89 0.955], 'HeadLength', 5, 'HeadWidth', 5, 'LineWidth', 0.8);
fig_text(fig, 18.55, 17.25, 'tension', 'HorizontalAlignment', 'center', 'FontSize', S.fsSmall);
annotation(fig, 'arrow', [0.925 0.99], [0.845 0.845], 'HeadLength', 5, 'HeadWidth', 5, 'LineWidth', 0.8);
fig_text(fig, 18.55, 16.25, 'shear', 'HorizontalAlignment', 'center', 'FontSize', S.fsSmall);
panel_caption(axa{3}, 'a', 'Interface families: periodic unit cells with printed roads (raster 0\circ)', 0.75);

% ---------------------------------------------------------------- (b) design space
rowsG = {{'sine_ar0.5', 'sine_ar1', 'sine_ar1.5', 'sine_ar2'}, {'U_ar0.5', 'U_ar1', 'U_ar1.5', 'U_ar2'}, ...
         {'T_ar0.5', 'T_ar1', 'T_ar1.5', 'T_ar2'}, ...
         {'dovetail_ar0.5_a10', 'dovetail_ar1_a10', 'dovetail_ar1.5_a10', 'dovetail_ar2_a10'}, ...
         {'dovetail_ar0.5_a20', 'dovetail_ar1_a20', 'dovetail_ar1.5_a20', 'dovetail_ar2_a20'}, ...
         {'dovetail_ar0.5_a30', 'dovetail_ar1_a30', 'dovetail_ar1.5_a30', 'dovetail_ar2_a30'}};
rowL = {'sinusoidal', 'U', 'T', 'dovetail, \alpha = 10\circ', 'dovetail, \alpha = 20\circ', 'dovetail, \alpha = 30\circ'};
gx0 = 3.3; gy0 = 12.9; vals = [0.5 1 1.5 2]; iw = 1.45; ih = 0.9; dx = 1.62; dyy = 1.02;
for c = 1:4
    fig_text(fig, gx0 + (c - 1) * dx + iw / 2, gy0 + 0.35, sprintf('d/w = %g', vals(c)), ...
             'HorizontalAlignment', 'center', 'FontSize', S.fsSmall);
end
for r = 1:6
    yb = gy0 - r * dyy;
    fig_text(fig, gx0 - 0.15, yb + ih / 2, rowL{r}, 'HorizontalAlignment', 'right', 'FontSize', S.fsSmall);
    for c = 1:4
        g = sel(GE.geometry, rowsG{r}{c});
        ia = add_axes(fig, [gx0 + (c - 1) * dx, yb, iw, ih]);
        W = GE.W(find(g, 1));
        fill_roads(ia, GE.x(g), GE.y(g), W, -1.4, 2.3, S, false);
        set(ia, 'XLim', [0 W], 'YLim', [-1.4 2.3], 'XTick', [], 'YTick', [], 'Box', 'on', 'LineWidth', 0.4);
    end
end
fig_text(fig, gx0 + 1.5 * dx + iw / 2, gy0 - 6 * dyy - 0.35, 'The 24 geometries of the study (plus flat reference)', ...
         'HorizontalAlignment', 'center', 'FontSize', S.fs);
fig_text(fig, gx0 + 1.5 * dx + iw / 2, gy0 - 6 * dyy - 0.8, '(b)', 'HorizontalAlignment', 'center', ...
         'FontSize', S.fsPanel, 'FontWeight', 'bold');

% ---------------------------------------------------------------- (c) candidate paths
ax = add_axes(fig, [11.2 7.9 4.6 4.2]);
W = max(I.x);
fill_roads(ax, I.x, I.y, W, -1.3, 2.0, struct('matA', S.matA * 0.5 + 0.5, 'matB', S.matB * 0.5 + 0.5), false);
pids = unique(P.pid, 'stable');
segcol = [S.blue; S.red; S.green; S.purple; S.orange; S.sky];
for k = 1:numel(pids)
    r = sel(P.pid, pids{k}); xy = [P.x(r) P.y(r)];
    if pids{k}(1) == 'I'
        d = diff(xy); nrm = [-d(:, 2) d(:, 1)] ./ sqrt(sum(d .^ 2, 2)); nrm = [nrm; nrm(end, :)];
        sgn = 1 - 2 * (pids{k}(end) == 'b');
        xo = xy + 0.06 * sgn * nrm;
        c = segcol(mod(str2double(pids{k}(2:end - 1)) - 1, 6) + 1, :);
        plot(ax, xo(:, 1), xo(:, 2), '-', 'Color', c, 'LineWidth', 1.3);
        arrowhead(ax, xo(end, :), xo(end, :) - xo(end - 1, :), 0.14, c);
    else
        plot(ax, xy(:, 1), xy(:, 2), '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.5);
        arrowhead(ax, xy(end, :), xy(end, :) - xy(1, :), 0.10, [0.45 0.45 0.45]);
    end
end
corners = [1.5 0; 2.5 0; 1.13603 1; 2.86397 1];
plot(ax, corners(:, 1), corners(:, 2), 'o', 'MarkerSize', 4.5, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
gov = ON.pid{1};
r = sel(P.pid, gov); gp = [P.x(r) P.y(r)];
s = [0; cumsum(sqrt(sum(diff(gp) .^ 2, 2)))];
for a = L.da(sel(L.pid, gov))'
    q = [interp1(s, gp(:, 1), a) interp1(s, gp(:, 2), a)];
    plot(ax, [q(1) q(1)], [q(2) - 0.05 q(2) - 0.2], 'k-', 'LineWidth', 0.45);
end
axis(ax, 'equal'); set(ax, 'XLim', [0 W], 'YLim', [-1.3 2.0], 'XTick', 0:4, 'YTick', -1:2);
xlab(ax, 'x / w'); ylab(ax, 'y / w');
key(ax, {{'-', 'Color', S.blue, 'LineWidth', 1.3}, 'interface paths (both directions)'; ...
         {'--', 'Color', [0.45 0.45 0.45]}, 'bulk paths (fans from corners)'; ...
         {'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 4.5}, 'singular points'; ...
         {'-', 'Color', 'k', 'LineWidth', 0.45}, 'tick marks: 12 crack lengths'}, 'Location', 'northeast');
lg = findobj(fig, 'Type', 'legend'); set(lg(1), 'Units', 'centimeters', 'Position', [15.95 10.2 3.0 1.8]);
panel_caption(ax, 'c', 'Candidate crack paths (dovetail, d/w = 1, \alpha = 20\circ)');

% ---------------------------------------------------------------- (d) coupled criterion
ax = add_axes(fig, [2.2 2.4 9.6 3.0]);
r = sel(CV.pid, gov); a = CV.da_mm(r); ss = CV.sigma_stress(r); sg = CV.sigma_energy(r);
env = max(ss, sg); [smin, j] = min(env); ytop = 6 * smin;
patch(ax, [a; flipud(a)], [env; ytop * ones(size(a))], [0.92 0.92 0.92], 'EdgeColor', 'none');
r2 = sel(CV.pid, ON.pid{2});
plot(ax, CV.da_mm(r2), max(CV.sigma_stress(r2), CV.sigma_energy(r2)), '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.7);
plot(ax, a, ss, '-', 'Color', S.blue, 'LineWidth', 1.3);
plot(ax, a, sg, '-', 'Color', S.red, 'LineWidth', 1.3);
plot(ax, a(j), smin, 'p', 'MarkerSize', 10, 'MarkerFaceColor', S.yellow, 'MarkerEdgeColor', 'k', 'LineWidth', 0.6);
set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', [min(a) max(a)], 'YLim', [0.5 * smin ytop]);
xlab(ax, 'crack increment  \Delta\ita\rm  [mm]'); ylab(ax, '\sigma  [MPa]');
key(ax, {{'-', 'Color', S.blue, 'LineWidth', 1.3}, 'stress condition \sigma_s(\Delta\ita\rm)'; ...
         {'-', 'Color', S.red, 'LineWidth', 1.3}, 'energy condition \sigma_G(\Delta\ita\rm)'; ...
         {'-', 'Color', [0.6 0.6 0.6]}, sprintf('envelope of path %s', ON.pid{2}); ...
         {'s', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0.92 0.92 0.92], 'MarkerSize', 7}, 'both conditions met'; ...
         {'p', 'MarkerFaceColor', S.yellow, 'MarkerEdgeColor', 'k', 'MarkerSize', 8}, ...
         sprintf('onset: %.1f MPa, \\Delta\\ita\\rm_c = %.2f mm', smin, a(j))}, 'Location', 'eastoutside');
set(ax, 'Units', 'centimeters', 'Position', [2.2 2.4 9.6 3.0]);
lg = findobj(fig, 'Type', 'legend'); set(lg(1), 'Units', 'centimeters', 'Position', [12.2 2.45 6.2 2.9]);
panel_caption(ax, 'd', sprintf('Coupled stress-energy criterion of path %s (w = %g mm, k_\\sigma = %g, k_G = %g)', ...
    gov, ON.w_mm(1), ON.k_sigma(1), ON.k_G(1)), 1.15);

save_figure(fig, 'Figure_1');
