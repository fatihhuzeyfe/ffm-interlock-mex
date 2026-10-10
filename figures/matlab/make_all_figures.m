%MAKE_ALL_FIGURES  Draw all figures (EPS + PDF + PNG in figures/output_matlab).
%   Data: figures/data (CSV from  python analysis/export_figure_data.py;  odb_*.json from
%   python scripts/export_figure_meshes.py). MATLAB R2016b or newer, no toolboxes.
close all;
Fig1_concept;
Fig2_master_curve;
Fig3_mechanisms_size;
Fig4_sensitivity_design;
FigS1_verification;
if exist(fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'fig5b_czm_summary.csv'), 'file')
    Fig5_czm;
end
