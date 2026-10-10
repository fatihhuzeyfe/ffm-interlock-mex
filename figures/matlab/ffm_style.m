function S = ffm_style()
%FFM_STYLE  Common settings of all figures. Edit here to restyle every figure at once.
%   All lengths in centimetres. Fonts: Helvetica (standard PostScript font, stays editable
%   text in CorelDRAW / Illustrator when the EPS is imported).
here     = fileparts(mfilename('fullpath'));
S.data   = fullfile(here, '..', 'data');
S.out    = fullfile(here, '..', 'output_matlab');
if ~exist(S.out, 'dir'), mkdir(S.out); end
S.font   = 'Helvetica';
S.fs     = 7;          % axis labels and ticks [pt]
S.fsSmall= 6;          % annotations
S.fsPanel= 9;          % panel letters
S.lw     = 0.6;        % axes line width
S.W      = 19.0;       % full page width [cm]
% Okabe-Ito colour-blind safe palette
S.blue   = [0 114 178]/255;   S.orange = [230 159 0]/255;  S.green = [0 158 115]/255;
S.red    = [213 94 0]/255;    S.purple = [204 121 167]/255; S.sky  = [86 180 233]/255;
S.yellow = [240 228 66]/255;  S.grey   = [0.50 0.50 0.50];
S.pairs   = {'PLA/PETG', 'PLA/ABS', 'PETG/ABS'};
S.pairCol = [S.blue; S.red; S.green];
S.matA   = [201 220 235]/255;  S.matB = [243 220 194]/255;      % material A (lower), B (upper)
S.roadA  = [157 181 203]/255;  S.roadB = [217 184 146]/255;     % printed roads
S.shapeCol = struct('T', S.purple, 'U', S.orange, 'dovetail', S.blue, 'sine', S.green, 'flat', S.grey);
end
