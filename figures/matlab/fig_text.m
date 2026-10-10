function h = fig_text(fig, x, y, s, varargin)
%FIG_TEXT  Text at a position in centimetres on the figure (rotation allowed), black, figure font.
S = ffm_style();
ov = getappdata(fig, 'cmOverlay');
if isempty(ov) || ~isgraphics(ov)
    u = get(fig, 'Units'); set(fig, 'Units', 'centimeters'); p = get(fig, 'Position'); set(fig, 'Units', u);
    ov = axes('Parent', fig, 'Units', 'centimeters', 'Position', [0 0 p(3) p(4)], 'XLim', [0 p(3)], ...
              'YLim', [0 p(4)], 'Visible', 'off', 'HitTest', 'off');
    setappdata(fig, 'cmOverlay', ov);
end
h = text(ov, x, y, s, 'FontName', S.font, 'FontSize', S.fs, 'Color', 'k', 'Clipping', 'off', varargin{:});
end
