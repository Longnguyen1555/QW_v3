% Demo: linewidth versus temperature, with digitized experimental references.
% Run this script from the folder containing moap2.m.
% CAUTION: gamma is a phenomenological regularization, not a fitted rate.
cfg = struct();
cfg.run.task = 'linewidth_sweep';
cfg.sweep.parameter = 'T_K';
cfg.sweep.values = [5 25 50 75 100 125 150 175 200 225 250 275 300];
cfg.fields.B_T = 10;
cfg.structure.Lz_nm = 10;
cfg.oap.broadening_meV = 0.5;
cfg.linewidth.level_fraction = 0.5;
cfg.linewidth.quarter_fraction = 0.25;
cfg.linewidth.show_experiment = true;
cfg.linewidth.plot_quarter = true;
r = moap2(cfg);
