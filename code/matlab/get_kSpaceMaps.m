function [refF2, errormapF2] = get_kSpaceMaps(ref, img)

    refF = fft3c(ref);
    refF2 = abs(refF).^2;

    errormapF2 = abs(fft3c(img) - refF).^2;

end