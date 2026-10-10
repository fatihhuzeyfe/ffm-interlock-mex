function panel_caption(ax, letter, description, below)
%PANEL_CAPTION  Panel description (centred, normal) and panel letter (centred, bold) BELOW the
%   axes. below = distance from the bottom of the axes to the description [cm] (room for the
%   tick labels and the x label; default 1.15 cm).
S = ffm_style();
if nargin < 4, below = 1.15; end
u = get(ax, 'Units'); set(ax, 'Units', 'centimeters'); p = get(ax, 'Position'); set(ax, 'Units', u);
fig = get(ax, 'Parent');
w = max(p(3), 4);
if ~isempty(description)
    annotation(fig, 'textbox', 'Units', 'centimeters', 'Position', [p(1) + p(3) / 2 - w / 2, p(2) - below - 0.42, w, 0.42], ...
        'String', description, 'FontName', S.font, 'FontSize', S.fs, 'Color', 'k', 'LineStyle', 'none', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Margin', 0);
    yl = p(2) - below - 0.85;
else
    yl = p(2) - below - 0.42;
end
annotation(fig, 'textbox', 'Units', 'centimeters', 'Position', [p(1) + p(3) / 2 - 0.5, yl, 1.0, 0.42], ...
    'String', ['(' letter ')'], 'FontName', S.font, 'FontSize', S.fsPanel, 'FontWeight', 'bold', 'Color', 'k', ...
    'LineStyle', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Margin', 0);
end
