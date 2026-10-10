%FIG4_SENSITIVITY_DESIGN  (a, b) Sobol indices of ln(gain) in tension and shear, (c) surrogate
%   prediction for geometries left out of the training, (d-g) Gaussian-process design maps of
%   the dovetail (PLA/PETG) in the strength- and toughness-controlled regimes.
S = ffm_style();
SB = read_csv(fullfile(S.data, 'fig4a_sobol.csv'));
DM = read_csv(fullfile(S.data, 'fig4b_design_map.csv'));
PA = read_csv(fullfile(S.data, 'fig4c_surrogate_parity.csv'));
CV = read_csv(fullfile(S.data, 'fig4_cv_scores.csv'));

fig = new_figure(S.W, 15.4);
loads = {'tension', 'shear'};
order = {'k_sigma', 'geometry', 'k_G', 'width_mm', 'pair', 'GIIc_over_GIc', 'raster'};
names = {'strength k_\sigma', 'geometry', 'toughness k_G', 'size w', 'material pair', 'G_{IIc}/G_{Ic}', 'raster'};

% ---------------------------------------------------------------- (a, b) Sobol
for L = 1:2
    ax = add_axes(fig, [2.7 + (L - 1) * 5.1, 9.6, 4.2, 4.6]);
    for k = 1:7
        r = sel(SB.load, loads{L}) & sel(SB.factor, order{k});
        st = SB.ST(r); s1 = SB.S1(r);
        patch(ax, [0 st st 0], [k - 0.32 k - 0.32 k + 0.32 k + 0.32], [0.82 0.82 0.82], 'EdgeColor', 'none');
        patch(ax, [0 s1 s1 0], [k - 0.32 k - 0.32 k + 0.32 k + 0.32], S.blue, 'EdgeColor', 'none');
        txt(ax, st + 0.015, k, sprintf('%.2f', st), 'FontSize', 5.3);
    end
    set(ax, 'YDir', 'reverse', 'YLim', [0.4 7.6], 'XLim', [0 0.75], 'YTick', 1:7, 'TickLength', [0.012 0]);
    if L == 1, set(ax, 'YTickLabel', names); else, set(ax, 'YTickLabel', {}); end
    xlab(ax, 'Sobol index of ln(gain)');
    if L == 1
        key(ax, {{'s', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', S.blue, 'MarkerSize', 6}, 'first order S_1'; ...
                 {'s', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0.82 0.82 0.82], 'MarkerSize', 6}, 'total S_T'}, ...
            'Location', 'southeast');
    end
    panel_caption(ax, char('a' + L - 1), sprintf('Sensitivity, %s', loads{L}));
end

% ---------------------------------------------------------------- (c) surrogate parity (subsample)
ax = add_axes(fig, [13.9 9.6 4.6 4.6]);
rng(1); idx = randperm(numel(PA.gain_ffm), min(4000, numel(PA.gain_ffm)));
tn = sel(PA.load(idx), 'tension');
plot(ax, PA.gain_ffm(idx(tn)), PA.gain_surrogate(idx(tn)), '.', 'Color', S.blue * 0.6 + 0.4, 'MarkerSize', 3);
plot(ax, PA.gain_ffm(idx(~tn)), PA.gain_surrogate(idx(~tn)), '.', 'Color', S.red * 0.6 + 0.4, 'MarkerSize', 3);
plot(ax, [0.75 2.9], [0.75 2.9], 'k-', 'LineWidth', 0.8);
axis(ax, 'square'); set(ax, 'XLim', [0.75 2.9], 'YLim', [0.75 2.9], 'XTick', 1:0.5:2.5, 'YTick', 1:0.5:2.5);
xlab(ax, 'gain, FFM'); ylab(ax, 'gain, surrogate');
r2 = CV.R2_leave_one_geometry_out(sel(CV.load, 'both'));
txt(ax, 0.85, 2.78, {sprintf('R^2 = %.2f', r2), '(4000 random points)'}, 'VerticalAlignment', 'top');
key(ax, {{'.', 'Color', S.blue, 'MarkerSize', 8}, 'tension'; {'.', 'Color', S.red, 'MarkerSize', 8}, 'shear'}, 'Location', 'southeast');
panel_caption(ax, 'c', 'Surrogate, unseen geometries');

% ---------------------------------------------------------------- (d-g) design maps
lws = [0.5 20]; reg = {'strength-controlled', 'toughness-controlled'};
k = 0; axes4 = cell(1, 4);
for L = 1:2
    for j = 1:2
        k = k + 1;
        ax = add_axes(fig, [1.5 + (k - 1) * 4.05, 2.4, 3.55, 3.55]); axes4{k} = ax;
        r = sel(DM.load, loads{L}) & abs(DM.lch_over_w - lws(j)) < 1e-9;
        ar = unique(DM.ar(r)); al = unique(DM.alpha_deg(r));
        Zg = reshape(DM.gain(r), numel(ar), numel(al))';
        contourf(ax, ar, al, Zg, 0.8:0.1:2.6, 'LineStyle', 'none');
        [Cc, hc] = contour(ax, ar, al, Zg, [1 1.25 1.5 1.75 2 2.25], 'LineColor', 'k', 'LineWidth', 0.5);
        clabel(Cc, hc, 'FontName', S.font, 'FontSize', 5, 'Color', 'k', 'LabelSpacing', 200);
        colormap(ax, parula(256)); caxis(ax, [0.8 2.6]);
        [gx, gy] = meshgrid([0.5 1 1.5 2], [10 20 30]);
        plot(ax, gx(:), gy(:), 'o', 'MarkerSize', 2.5, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 0.4, 'Clipping', 'off');
        set(ax, 'XLim', [0.5 2], 'YLim', [10 30], 'XTick', [0.5 1 1.5 2], 'YTick', 10:5:30);
        xlab(ax, 'aspect ratio d/w');
        if k == 1, ylab(ax, 'flank angle \alpha [\circ]'); else, set(ax, 'YTickLabel', {}); end
        title(ax, sprintf('\\itl\\rm_{ch}/\\itw\\rm = %g', lws(j)), 'FontName', S.font, 'FontSize', 6, 'FontWeight', 'normal', 'Color', 'k');
        panel_caption(ax, char('d' + k - 1), sprintf('%s, %s', loads{L}, reg{j}));
    end
end
cb = colorbar(axes4{4}, 'Units', 'centimeters', 'Position', [17.9 2.4 0.22 3.55]);
set(cb, 'FontName', S.font, 'FontSize', 5.5, 'Color', 'k');
ylabel(cb, 'gain vs. flat', 'FontSize', 6, 'FontWeight', 'bold');
set(axes4{4}, 'Units', 'centimeters', 'Position', [1.5 + 3 * 4.05, 2.4, 3.55, 3.55]);
fig_text(fig, 9.5, 6.75, 'Design maps of the dovetail (PLA/PETG): white dots = FE geometries, lines = iso-gain', ...
         'HorizontalAlignment', 'center', 'FontSize', S.fsSmall);

save_figure(fig, 'Figure_4');
