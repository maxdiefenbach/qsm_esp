%
%	qsm2016_script_recon_evaluation.m
%
%	Evaluation script for the QSM reconstruction challenge 2016
%
%	Find the results and report at http://qsm.neuroimaging.at
%
%	2016 BB, CL, FS


%%-------------------------------------------------------------------------
%% load data
%%-------------------------------------------------------------------------
% set(0,'DefaultFigureWindowStyle','docked')
addpath data/
addpath utils/

load phs_tissue;            % tissue phase from transversal orientation (in ppm, normalized by gyro*TE*B0)
load phs_wrap;              % raw phase from trans orient without processing (in radians)
load phs_unwrap;            % phase from trans orient processed with Laplacian unwrapping and BET masking (in radians)
load spatial_res;           % voxel size
load msk;                   % brain mask => obtained by eroding the BET mask by 5 voxels (by setting peel=5 in LBV)
load magn;                  % magnitude from transversal orientation
load magn_raw;              % raw magnitude without masking from transversal orientation
load mp_rage;               % mprage
load chi_33;                % chi_33 from STI solution (in ppm)
load chi_cosmos;            % COSMOS from 12 orientations (in ppm)
N = size(msk);

%%-------------------------------------------------------------------------
%% display data
%%-------------------------------------------------------------------------
imagesc3d2(chi_cosmos, N/2, 1, [90,90,-90], [-0.10,0.14], [], 'Cosmos')
imagesc3d2(chi_33, N/2, 2, [90,90,-90], [-0.10,0.14], [], '\chi_{33}')
imagesc3d2(msk, N/2, 3, [90,90,-90], [0,1], [], 'Mask')
imagesc3d2(phs_tissue, N/2, 4, [90,90,-90], [-0.05,0.05], [], 'Input Phase')
imagesc3d2(phs_wrap, N/2, 5, [90,90,-90], [-pi,pi], [], 'Raw Phase')
imagesc3d2(phs_unwrap, N/2, 6, [90,90,-90], [-pi,pi], [], 'Unwrapped Phase')
imagesc3d2(magn, N/2, 7, [90,90,-90], [0,0.5], [], 'Magnitude')
imagesc3d2(magn_raw, N/2, 8, [90,90,-90], [0,0.5], [], 'Raw Magnitude')
imagesc3d2(mp_rage, N/2, 9, [90,90,-90], [0,0.9], [], 'MPRAGE')

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
mosaic(squeeze(kernel_disp(1+end/2,:,:)), 1, 1, 11, '', [0,2/3])

%%-------------------------------------------------------------------------
%% TKD recon
%%-------------------------------------------------------------------------

thre_tkd = 0.19;      % TKD threshold parameter

kernel_inv = zeros(N);
kernel_inv( abs(kernel) > thre_tkd ) = 1 ./ kernel(abs(kernel) > thre_tkd);

chi_tkd = real( ifftn( fftn(phs_tissue) .* kernel_inv ) ) .* msk;
imagesc3d2(chi_tkd, N/2, 21, [90,90,-90], [-0.10,0.14], [], 'TKD')

%%-------------------------------------------------------------------------
%% closed-form L2 recon
%%-------------------------------------------------------------------------

[k1, k2, k3] = ndgrid(0:N(1)-1,0:N(2)-1,0:N(3)-1);

E1 = 1 - exp(2i .* pi .* k1 / N(1));
E2 = 1 - exp(2i .* pi .* k2 / N(2));
E3 = 1 - exp(2i .* pi .* k3 / N(3));

EtE = abs(E1).^2 + abs(E2).^2 + abs(E3).^2;
DtD = abs(kernel).^2;
reg_param = 9e-2;   % gradient regularization parameter
chi_L2 = real( ifftn(conj(kernel) .* fftn(phs_tissue) ./ (DtD + reg_param * EtE))) .* msk;
imagesc3d2(chi_L2, N/2, 22, [90,90,-90], [-0.10,0.14], [], 'CF L2')

%%-------------------------------------------------------------------------
%% load data
%%-------------------------------------------------------------------------
load chi_33;                % chi_33 from STI solution (in ppm)
load chi_cosmos;            % COSMOS from 12 orientations (in ppm)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% LOAD AND ADD YOUR QSM TO THE LISTS BELOW %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% create 4d volume of all results
qsm_results = cat(4, chi_tkd, chi_L2);
qsm_names = {'TKD'; 'CFL2'; };
disp([ 'totally ', num2str(size(qsm_results, 4)), ' algorithms loaded']);

%%-------------------------------------------------------------------------
%% compute and print performance metrics
%%-------------------------------------------------------------------------
disp('evaluation results');
load evaluation_mask;
wm_mask = (evaluation_mask>6);
wm_px = sum(wm_mask(:));
gm_mask = (evaluation_mask>0) & (evaluation_mask<=6);
gm_px = sum(gm_mask(:));

disp('   Algorithm: RSME  : HFEN  : SSIM   : WM ERROR : GM ERROR : W+GM ERROR : (all relative to chi33)')
for ccc = 1:(size(qsm_results, 4))
    gm = qsm_results(:,:,:,ccc);
    wm = qsm_results(:,:,:,ccc);
    gm_error(ccc) = sum(abs((gm(gm_mask) - chi_33(gm_mask)))) / gm_px ;
    wm_error(ccc) = sum(abs((wm(wm_mask) - chi_33(wm_mask)))) / wm_px ;
    rmse(ccc) = compute_rmse(qsm_results(:,:,:,ccc), chi_33);
    hfen(ccc) = compute_hfen(qsm_results(:,:,:,ccc), chi_33);
    ssim(ccc) = compute_ssim(qsm_results(:,:,:,ccc), chi_33);
    mean_error(ccc) = (gm_error(ccc) + wm_error(ccc))/2;

    disp([ sprintf('%15s ', char(qsm_names(ccc))), ...
        ' : ', num2str(rmse(ccc), '%-5.6f'), ...
        ' : ', num2str(hfen(ccc), '%5.6f'), ...
        ' : ', num2str(ssim(ccc), '%5.6f'), ...
        '  : ', num2str(wm_error(ccc), '%5.3f'), ...
        '    : ', num2str(gm_error(ccc), '%5.3f'), ...
        '    : ', num2str((gm_error(ccc) + wm_error(ccc))/2 , '%5.6f')]);
end




% mnd: save maps with correct orientations
phs_tissue = flip(permute(phs_tissue, [2, 1, 3]), 3);
mp_rage = flip(permute(mp_rage, [2, 1, 3]), 3);
mask = flip(permute(msk, [2, 1, 3]), 3);
evaluation_mask = flip(permute(evaluation_mask, [2, 1, 3]), 3);

chi_33 = flip(permute(chi_33, [2, 1, 3]), 3);
chi_cosmos = flip(permute(chi_cosmos, [2, 1, 3]), 3);
chi_tkd = flip(permute(chi_tkd, [2, 1, 3]), 3);
chi_L2 = flip(permute(chi_L2, [2, 1, 3]), 3);

save('baseline.mat', ...
     'phs_tissue', 'mp_rage', 'mask', 'evaluation_mask', 'spatial_res', ...
     'chi_cosmos', 'chi_33', 'chi_tkd', 'chi_L2')