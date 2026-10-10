function panel_label(ax, txt, dx, dy)
%PANEL_LABEL  Bold panel letter, dx/dy = offset of the letter from the top-left axes corner [cm].
S = ffm_style();
if nargin < 3, dx = -0.9; end
if nargin < 4, dy = 0.15; end
u = get(ax, 'Units'); set(ax, 'Units', 'centimeters'); p = get(ax, 'Position'); set(ax, 'Units', u);
annotation(get(ax, 'Parent'), 'textbox', 'Units', 'centimeters', ...
    'Position', [p(1) + dx, p(2) + p(4) + dy, 0.6, 0.5], 'String', txt, 'FontName', S.font, ...
    'FontSize', S.fsPanel, 'FontWeight', 'bold', 'LineStyle', 'none', 'VerticalAlignment', 'bottom', ...
    'Margin', 0);
end
