function ax = add_axes(fig, pos)
%ADD_AXES  Axes at an exact position [left bottom width height] in centimetres.
S = ffm_style();
ax = axes('Parent', fig, 'Units', 'centimeters', 'Position', pos, 'FontName', S.font, ...
          'FontSize', S.fs, 'LineWidth', S.lw, 'TickDir', 'in', 'Box', 'on', 'Layer', 'top', ...
          'TickLength', [0.018 0.018], 'XColor', 'k', 'YColor', 'k');
hold(ax, 'on');
end
