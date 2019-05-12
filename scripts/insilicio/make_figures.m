close all; clear all; clc; pwd
addpath(genpath('../../'))

load('baseline.mat')
load('results_baseline_phantom.mat')
load('numericalbrainphantom.mat')

chi_ref = numericalbrainphantom.chimap_ppm;
B0dir = numericalbrainphantom.B0dir;

SNR = 100
range = max(numericalbrainphantom.RDF_ppm(:)) - min(numericalbrainphantom.RDF_ppm(:))
noise = range/SNR * randn(size(numericalbrainphantom.RDF_ppm));

spatial_res = numericalbrainphantom.voxelSize_mm;
N = size(numericalbrainphantom.chimap_ppm);
mask = numericalbrainphantom.mask;
phs_tissue = numericalbrainphantom.RDF_ppm + mask .* noise;


%%%%%%%%%%%%%%%%%%
% create figures %
%%%%%%%%%%%%%%%%%%
system('rm *.png *.csv *.aux *.log')

close all;

% chi ground truth
clim = [-0.05, 0.2]
iy = 84
plot_orientations(chi_ref, 'iy', iy, 'CLim', clim)
create_tikzPDF('chi_ref')

% noisy RDF
plot_orientations(phs_tissue, 'iy', iy, 'CLim', [-0.05, 0.05])
create_tikzPDF('RDF_noisy')

% tkd
imagine(results_tkd.chimaps_ppm{:})
h = plot_metrics([results_tkd.metrics{:}])
write_csv(h, 'tkd_hyperparam_selection')
create_tikzPDF('L2_hyperparam_selection')
for i = [2, 3, 4, 6]
    plot_orientations(results_tkd.chimaps_ppm{i}, 'CLim', clim)
    create_tikzPDF(['chi_tkd_lambda' num2str(i)])
end

% L2
imagine(results_L2.chimaps_ppm{:})
h = plot_metrics([results_L2.metrics{:}])
write_csv(h, 'L2_hyperparam_selection')
create_tikzPDF('L2_hyperparam_selection')
for i = [2, 3, 4, 6]
    plot_orientations(results_L2.chimaps_ppm{i}, 'CLim', clim)
    create_tikzPDF(['chi_L2_lambda' num2str(i)])
end


% fuse shells two L2s
chi1 = results_L2.chimaps_ppm{1};
chi2 = results_L2.chimaps_ppm{6};
metrics1 = results_L2.metrics{1};
metrics2 = results_L2.metrics{6};
chi12 = fuse_shell(chi1, chi2, metrics1, metrics2);
metrics12 = compute_metrics(chi_ref, chi12);

imagine(chi_ref, chi1, chi2, chi12)
h = plot_metrics([metrics1, metrics2, metrics12])
write_csv(h, 'fuse_L2_lambda1_lambda6')
create_tikzPDF('fuse_L2_lambda1_lambda6')
plot_orientations(chi12, 'CLim', clim)
create_tikzPDF('chi12_L2_lambda1_lambda6')


% fuse tkd and L2
chi1 = results_tkd.chimaps_ppm{2};
chi2 = results_L2.chimaps_ppm{4};
metrics1 = results_tkd.metrics{2};
metrics2 = results_L2.metrics{4};
chi12 = fuse_shell(chi1, chi2, metrics1, metrics2);
metrics12 = compute_metrics(chi_ref, chi12);

imagine(chi_ref, chi1, chi2, chi12)
h = plot_metrics([metrics1, metrics2, metrics12])
write_csv(h, 'fuse_tkd2_L24')
create_tikzPDF('fuse_tkd2_L24')
plot_orientations(chi12, 'CLim', clim)
create_tikzPDF('chi12_tkd2_L24')


% fuse tkd l2 on zcs
chi1 = results_tkd.chimaps_ppm{2};
chi2 = results_L2.chimaps_ppm{4};
metrics1 = results_tkd.metrics{2};
metrics2 = results_L2.metrics{4};
chi12 = fuse_zcd(B0dir, chi1, chi2, metrics1, metrics2);
metrics12 = compute_metrics(chi_ref, chi12);

imagine(chi_ref, chi1, chi2, chi12)
h = plot_metrics([metrics1, metrics2, metrics12])
write_csv(h, 'fuse_tkd2_L24_zcd')
create_tikzPDF('fuse_tkd2_L24_zcd')
plot_orientations(chi12, 'CLim', clim)
create_tikzPDF('chi12_tkd2_L24_zcd')
plot_orientations(chi12, 'CLim', clim, 'iz', 88)
create_tikzPDF('chi12_tkd2_L24_iz88_zcd')


system('rm *.aux *.log')

close all;
