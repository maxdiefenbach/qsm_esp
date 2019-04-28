%% Fourier transform
function [out] = ftn(in)
    out = fftshift(fftn(ifftshift(in)));
end