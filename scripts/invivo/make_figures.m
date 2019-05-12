close all; clear all; clc;
addpath(genpath('../../../'))

load('baseline.mat')
load('results_baseline.mat')

DataParams.voxelSize_mm = spatial_res;
DataParams.B0dir = [0, 0, 1];

%%%%%%%%%%%%%%%%%%
% create figures %
%%%%%%%%%%%%%%%%%%
system('rm *.png *.csv *.aux *.log')

close all;

% chi ground truth
plot_orientations(chi_33, 'CLim', [-0.1, 0.15])
create_tikzPDF('chi_33')
plot_orientations(chi_cosmos, 'CLim', [-0.1, 0.15])
create_tikzPDF('chi_cosmos')


% tkd
imagine(results_tkd.chimaps_ppm{:})
h = plot_metrics([results_tkd.metrics_chi_33{:}])
write_csv(h, 'tkd_hyperparam_selection')
create_tikzPDF('tkd_hyperparam_selection')
for i = 1:length(results_tkd.thre_tkd)
    plot_orientations(results_tkd.chimaps_ppm{i}, 'CLim', [-0.1, 0.15])
    create_tikzPDF(['chi_tkd_lambda' num2str(i)])
end

% L2
imagine(results_L2.chimaps_ppm{:})
h = plot_metrics([results_L2.metrics_chi_33{:}])
write_csv(h, 'L2_hyperparam_selection')
create_tikzPDF('L2_hyperparam_selection')
for i = 1:length(results_L2.reg_param)
    plot_orientations(results_L2.chimaps_ppm{i}, 'CLim', [-0.1, 0.15])
    create_tikzPDF(['chi_L2_lambda' num2str(i)])
end


% fuse shells two L2s
chi1 = results_L2.chimaps_ppm{1};
chi2 = results_L2.chimaps_ppm{3};
metrics1 = results_L2.metrics_chi_33{1};
metrics2 = results_L2.metrics_chi_33{3};
chi12 = fuse_shell(chi1, chi2, metrics1, metrics2);
metrics12_chi_33 = compute_metrics(chi_33, chi12);

imagine(chi_33, chi1, chi2, chi12)
h = plot_metrics([metrics1, metrics2, metrics12_chi_33])
write_csv(h, 'fuse_L2_lambdas')
create_tikzPDF('fuse_L2_lambdas')
plot_orientations(chi12, 'CLim', [-0.1, 0.15])
create_tikzPDF('chi12_L2_lambdas')


% fuse tkd and L2
chi1 = results_tkd.chimaps_ppm{3};
chi2 = results_L2.chimaps_ppm{2};
metrics1 = results_tkd.metrics_chi_33{3};
metrics2 = results_L2.metrics_chi_33{2};
chi12 = fuse_shell(chi1, chi2, metrics1, metrics2);
metrics12_chi_33 = compute_metrics(chi_33, chi12);

imagine(chi_33, chi1, chi2, chi12)
h = plot_metrics([metrics1, metrics2, metrics12_chi_33])
write_csv(h, 'fuse_tkd_L2')
create_tikzPDF('fuse_tkd_L2')
plot_orientations(chi12, 'CLim', [-0.1, 0.15])
create_tikzPDF('chi12_tkd_L2')

% fuse tkd l2 on zcs
chi1 = results_tkd.chimaps_ppm{3}
chi2 = results_L2.chimaps_ppm{2};
metrics1 = results_tkd.metrics_chi_33{3};
metrics2 = results_L2.metrics_chi_33{2};
chi12 = fuse_zcd(DataParams.B0dir, chi1, chi2, metrics1, metrics2);
metrics12_chi_33 = compute_metrics(chi_33, chi12);

imagine(chi_33, chi1, chi2, chi12)
h = plot_metrics([metrics1, metrics2, metrics12_chi_33])
write_csv(h, 'fuse_tkd_L2_zcd')
create_tikzPDF('fuse_tkd_L2_zcd')
plot_orientations(chi12, 'CLim', [-0.1, 0.15])
create_tikzPDF('chi12_tkd_L2_zcd')


system('rm *.aux *.log')

close all;