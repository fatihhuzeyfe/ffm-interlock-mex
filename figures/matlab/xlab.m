function xlab(ax, s)
%XLAB  Bold black x label.
S = ffm_style();
xlabel(ax, s, 'FontName', S.font, 'FontSize', S.fs, 'FontWeight', 'bold', 'Color', 'k');
end
