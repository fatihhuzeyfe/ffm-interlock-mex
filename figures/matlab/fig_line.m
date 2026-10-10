function fig_line(fig, x, y, varargin)
%FIG_LINE  Line in figure centimetres (uses the overlay axes of fig_text).
fig_text(fig, 0, 0, '');                       % make sure the overlay axes exists
ov = getappdata(fig, 'cmOverlay');
line(ov, x, y, 'Color', 'k', 'LineWidth', 0.6, 'Clipping', 'off', varargin{:});
end
