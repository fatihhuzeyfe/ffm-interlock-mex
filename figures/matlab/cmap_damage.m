function c = cmap_damage(n)
%CMAP_DAMAGE  Light yellow (undamaged) -> orange -> dark purple/black (fully damaged).
if nargin < 1, n = 256; end
a = [0.99 1.00 0.64; 0.98 0.60 0.03; 0.73 0.21 0.33; 0.34 0.06 0.43; 0.00 0.00 0.02];
c = interp1(linspace(0, 1, size(a, 1)), a, linspace(0, 1, n));
end
