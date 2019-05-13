close all; clear all; clc;

unzip('numericalbrainphantom.zip')

nii = load_nii('theo_susceptibility.nii.gz')

numericalbrainphantom.chimap_ppm = flip(permute(nii.img, [2, 1, 3]), 3);
numericalbrainphantom.voxelSize_mm = [nii.hdr.dime.pixdim(3), nii.hdr.dime.pixdim(2), nii.hdr.dime.pixdim(1)];
numericalbrainphantom.B0dir = [0, 0, 1];
numericalbrainphantom.fieldStrength_T = 3;
numericalbrainphantom.centerFreq_Hz = 42.58e6 * numericalbrainphantom.fieldStrength_T;
nii = load_nii('mask_data.nii.gz')
numericalbrainphantom.mask = flip(permute(nii.img, [2, 1, 3]), 3);
numericalbrainphantom.RDF_ppm = forwardSimulate_RDF_ppm(numericalbrainphantom);

numericalbrainphantom
imagine(numericalbrainphantom.chimap_ppm, ...
        numericalbrainphantom.mask, ...
        numericalbrainphantom.RDF_ppm)

save('numericalbrainphantom.mat', 'numericalbrainphantom')