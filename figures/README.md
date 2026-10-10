# Figures

`matlab/` contains the scripts of all figures of the article (MATLAB R2016b or newer, no toolboxes).
They read the data written to `figures/data/` by `analysis/export_figure_data.py`,
`scripts/czm_postprocess.py` and `scripts/export_figure_meshes.py` (step 2 of `docs/REPRODUCE.md`)
and write EPS, PDF and PNG files to `figures/output_matlab/`:

    cd figures/matlab
    make_all_figures

Every panel is placed at a fixed position in centimetres (`add_axes(fig, [left bottom width height])`);
fonts, sizes and colours of all figures are set in `matlab/ffm_style.m`.
