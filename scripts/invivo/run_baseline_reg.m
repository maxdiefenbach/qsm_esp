close all; clear all; clc;
addpath(genpath('../../'))

load('baseline.mat')

N = size(phs_tissue)

%%-------------------------------------------------------------------------
%% create dipole kernel
%%-------------------------------------------------------------------------
[ky,kx,kz] = meshgrid(-N(1)/2:N(1)/2-1, -N(2)/2:N(2)/2-1, -N(3)/2:N(3)/2-1);

kx = (kx / max(abs(kx(:)))) / spatial_res(1);
ky = (ky / max(abs(ky(:)))) / spatial_res(2);
kz = (kz / max(abs(kz(:)))) / spatial_res(3);
k2 = kx.^2 + ky.^2 + kz.^2;
R_tot = eye(3);     % orientation matrix for transverse acquisition
kernel = fftshift( 1/3 - (kx * R_tot(3,1) + ky * R_tot(3,2) + kz * R_tot(3,3)).^2 ./ (k2 + eps) );
kernel_disp = fftshift(mean(abs(kernel), 4));


%%-------------------------------------------------------------------------
%% TKD recon
%%-------------------------------------------------------------------------
thresholds = [0.19, 0.21, 0.26];
for i = 1:length(thresholds)
    thre_tkd = thresholds(i)
    kernel_inv = zeros(N);
    kernel_inv( abs(kernel) > thre_tkd ) = 1 ./ kernel(abs(kernel) > thre_tkd);

    chi_tkd = real( ifftn( fftn(phs_tissue) .* kernel_inv ) ) .* mask;

    results_tkd.thre_tkd{i} = thre_tkd;
    results_tkd.chimaps_ppm{i} = chi_tkd;
    results_tkd.metrics_chi_33{i} = compute_metrics(chi_33, chi_tkd)
    results_tkd.metrics_chi_cosmos{i} = compute_metrics(chi_cosmos, chi_tkd)
end


%%-------------------------------------------------------------------------
%% closed-form L2 recon
%%-------------------------------------------------------------------------
[k1, k2, k3] = ndgrid(0:N(1)-1,0:N(2)-1,0:N(3)-1);

E1 = 1 - exp(2i .* pi .* k1 / N(1));
E2 = 1 - exp(2i .* pi .* k2 / N(2));
E3 = 1 - exp(2i .* pi .* k3 / N(3));

EtE = abs(E1).^2 + abs(E2).^2 + abs(E3).^2;
DtD = abs(kernel).^2;

regularizationParameters = [6e-2, 9e-2, 2e-1]
for i = 1:length(regularizationParameters)
    reg_param = regularizationParameters(i)

    chi_L2 = real( ifftn(conj(kernel) .* fftn(phs_tissue) ./ (DtD + reg_param * EtE))) .* mask;

    results_L2.reg_param{i} = reg_param;
    results_L2.chimaps_ppm{i} = chi_L2;
    results_L2.metrics_chi_33{i} = compute_metrics(chi_33, chi_L2)
    results_L2.metrics_chi_cosmos{i} = compute_metrics(chi_cosmos, chi_L2)
end

save('results_baseline.mat', 'results_tkd', 'results_L2')
